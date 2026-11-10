"""Stage Numerics.dll and the pinned Roslyn assemblies and build the pynumerics wheel.

The wheel is pure Python plus ``Numerics.dll``, the two Roslyn assemblies
(``Microsoft.CodeAnalysis.dll`` and ``Microsoft.CodeAnalysis.CSharp.dll``, used by
``pynumerics.callbacks.compile_csharp``) and ``pynumerics.runtimeconfig.json``
under ``pynumerics/_bin``, so one ``py3-none-any`` wheel serves Linux and Windows; the
.NET 10 runtime is its only platform requirement. ``Numerics.dll`` is taken byte for byte
from ``lib/net10.0`` of the ``RMC.Numerics`` NuGet package whose version equals
``NUMERICS_VERSION`` in ``src/pynumerics/_version.py`` (or, with ``--dll``, from a locally
built assembly).

The identity of the staged assembly is written into the wheel's ``_version.py`` at build time:
``NUMERICS_VERSION`` (the package version, or the local build's informational version),
``NUMERICS_ASSEMBLY_VERSION`` (read from the assembly's version resource), ``NUMERICS_SHA256``
(of the file), and, for a local build, ``NUMERICS_SOURCE_COMMIT`` and ``NUMERICS_BUILD_COMMAND``
(empty for a package assembly). The version and assembly version must equal the values tracked in
``_version.py``, and a non-empty tracked ``NUMERICS_SHA256`` must equal the digest, so a release
cannot ship an assembly other than the one the source declares. The tracked file is never modified:
the wheel is built from a temporary copy of the package.

The package is found, in order: ``--nupkg PATH``; the NuGet global packages folder
(``NUGET_PACKAGES`` or ``~/.nuget/packages``), where ``dotnet restore`` of any project
referencing ``RMC.Numerics`` leaves it; with ``--download``, the flat container of ``--feed``
(nuget.org by default). With ``--dll PATH`` that file is staged instead; its source commit is taken
from the ``+<commit>`` suffix of its informational version unless ``--source-commit`` is given.
The staged assembly stays in ``_bin`` (gitignored) so a source checkout can run the tests against it.

The Roslyn assemblies come from ``lib/ROSLYN_TARGET_FRAMEWORK`` of the Microsoft.CodeAnalysis.Common
and Microsoft.CodeAnalysis.CSharp packages at ``ROSLYN_VERSION``, found in the NuGet global packages
folder or, with ``--download``, on nuget.org; each file's SHA-256 must equal its pin in ``_version.py``.
Their other dependencies (System.Collections.Immutable, System.Reflection.Metadata) are part of the
.NET 10 shared framework and are not bundled. ``--roslyn-only`` stages just the Roslyn assemblies.

Usage, from ``python/pynumerics``::

    python scripts/build_wheel.py [--nupkg PATH | --dll PATH] [--download] [--feed URL]
                                  [--stage-only | --roslyn-only] [--out DIR]

Needs the ``build`` package (``pip install build``) unless ``--stage-only``. Prints the
assembly identity, then the wheel's path and SHA-256 on its last lines.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import os
import pathlib
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import urllib.request
import xml.etree.ElementTree as ElementTree
import zipfile

PACKAGE = pathlib.Path(__file__).resolve().parents[1]
BIN = PACKAGE / "src" / "pynumerics" / "_bin"
VERSION_FILE = PACKAGE / "src" / "pynumerics" / "_version.py"
VERSION_ENTRY = "pynumerics/_version.py"
ASSEMBLY_ENTRY = "lib/net10.0/Numerics.dll"
RUNTIME_CONFIG = "pynumerics.runtimeconfig.json"
NUGET_FEED = "https://api.nuget.org/v3-flatcontainer"
PACKAGE_PATH = "{feed}/rmc.numerics/{v}/rmc.numerics.{v}.nupkg"
ROSLYN_CONTAINER = "https://api.nuget.org/v3-flatcontainer/{p}/{v}/{p}.{v}.nupkg"
ROSLYN_FILES = (
    ("microsoft.codeanalysis.common", "Microsoft.CodeAnalysis.dll", "ROSLYN_CODEANALYSIS_SHA256"),
    ("microsoft.codeanalysis.csharp", "Microsoft.CodeAnalysis.CSharp.dll", "ROSLYN_CSHARP_SHA256"),
)
"""(NuGet package id, assembly file, _version.py pin) of each bundled Roslyn assembly."""
ROSLYN_LICENSE = "Microsoft.CodeAnalysis.LICENSE.txt"
IDENTITY_NAMES = ("NUMERICS_VERSION", "NUMERICS_ASSEMBLY_VERSION", "NUMERICS_SHA256",
                  "NUMERICS_SOURCE_COMMIT", "NUMERICS_BUILD_COMMAND")
"""The _version.py constants written into the wheel from the staged assembly."""
DEFAULT_BUILD_COMMAND = "dotnet build Numerics/Numerics.csproj -c Release -f net10.0"
WHEEL_SOURCES = ("pyproject.toml", "README.md", "LICENSE", "src")
"""The files of the package a wheel is built from (copied to a temporary directory)."""


def _pin(name: str, text: str | None = None) -> str:
    """Read one string constant from _version.py (or the given source text) without importing it."""
    source = VERSION_FILE.read_text(encoding="utf-8") if text is None else text
    match = re.search(rf'^{name} = "([^"\r\n]*)"', source, re.MULTILINE)
    if match is None:
        raise SystemExit(f"{VERSION_FILE.name} lacks a one-line string {name}.")
    return match.group(1)


def _nuget_root() -> pathlib.Path:
    """Return the NuGet global packages folder."""
    return pathlib.Path(os.environ.get("NUGET_PACKAGES") or pathlib.Path.home() / ".nuget" / "packages")


def stage_roslyn(download: bool) -> list[str]:
    """Stage the pinned Roslyn assemblies into _bin from the NuGet cache (or nuget.org), verify them, return digests."""
    version, framework = _pin("ROSLYN_VERSION"), _pin("ROSLYN_TARGET_FRAMEWORK")
    root = _nuget_root()
    if not (BIN / ROSLYN_LICENSE).is_file():
        raise SystemExit(f"{BIN / ROSLYN_LICENSE} is missing; it is tracked in the repository.")
    digests = []
    for package, file, pin_name in ROSLYN_FILES:
        pinned = _pin(pin_name)
        cached = root / package / version / "lib" / framework / file
        if cached.is_file():
            assembly, source = cached.read_bytes(), str(cached)
        elif download:
            url = ROSLYN_CONTAINER.format(p=package, v=version)
            with urllib.request.urlopen(url, timeout=120) as response:  # noqa: S310 - fixed https nuget.org URL
                payload = response.read()
            with zipfile.ZipFile(io.BytesIO(payload)) as archive:
                assembly = archive.read(f"lib/{framework}/{file}")
            source = url
        else:
            raise SystemExit(f"{package} {version} is not in {root}; restore it or pass --download.")
        digest = hashlib.sha256(assembly).hexdigest()
        if digest != pinned:
            raise SystemExit(f"{file} from {source} has SHA-256 {digest}; _version.py pins {pinned}.")
        (BIN / file).write_bytes(assembly)
        print(f"staged {file} {digest} from {source}")
        digests.append(digest)
    return digests


def pep440(numerics_version: str) -> str:
    """Return the PEP 440 form of an RMC.Numerics version: ``2.3.0`` stays, ``2.3.0-dev.1`` is ``2.3.0.dev1``."""
    match = re.fullmatch(r"(\d+\.\d+\.\d+)(?:-dev\.(\d+))?", numerics_version)
    if match is None:
        raise SystemExit(f"Unsupported RMC.Numerics version {numerics_version!r} (expected X.Y.Z or X.Y.Z-dev.N).")
    return match.group(1) + (f".dev{match.group(2)}" if match.group(2) else "")


def version_resource_string(assembly: bytes, key: str) -> str:
    """Return one value of the string table in a PE file's version resource.

    The .NET SDK writes the assembly version (key ``Assembly Version``) and the informational
    version (key ``ProductVersion``) into the Win32 version resource of every assembly it builds.
    Each entry is a ``String`` structure: ``wLength``, ``wValueLength`` (in UTF-16 code units,
    including the terminating null) and ``wType`` (1 for text), the null-terminated UTF-16LE key,
    padding to a 32-bit boundary, then the null-terminated UTF-16LE value. Resource data is 32-bit
    aligned in the file, so the padding is computed from the file offset.
    """
    needle = (key + "\0").encode("utf-16-le")
    start = assembly.find(needle)
    while start >= 6:
        _, value_length, value_type = struct.unpack_from("<HHH", assembly, start - 6)
        value_start = (start + len(needle) + 3) & ~3
        if value_type == 1 and value_length > 0 and start % 2 == 0:
            value = assembly[value_start:value_start + 2 * value_length].decode("utf-16-le")
            return value.split("\0", 1)[0]
        start = assembly.find(needle, start + 2)
    raise SystemExit(f"The assembly has no {key!r} entry in its version resource.")


def _nuspec_version(archive: zipfile.ZipFile) -> str:
    """Return the <version> of the .nuspec at the root of a NuGet package."""
    names = [name for name in archive.namelist() if "/" not in name and name.endswith(".nuspec")]
    if len(names) != 1:
        raise SystemExit(f"Expected one .nuspec at the package root, found {names}.")
    root = ElementTree.fromstring(archive.read(names[0]))
    version = next((element.text for element in root.iter() if element.tag.rsplit("}", 1)[-1] == "version"), None)
    if not version:
        raise SystemExit(f"{names[0]} has no <version>.")
    return version.strip()


def _package_bytes(version: str, nupkg: pathlib.Path | None, download: bool, feed: str) -> tuple[bytes, str]:
    """Return the .nupkg bytes and a description of where they came from."""
    if nupkg is not None:
        return nupkg.read_bytes(), str(nupkg)
    root = _nuget_root()
    cached = root / "rmc.numerics" / version / f"rmc.numerics.{version}.nupkg"
    if cached.is_file():
        return cached.read_bytes(), str(cached)
    if not download:
        raise SystemExit(f"RMC.Numerics {version} is not in {root}; pass --nupkg, restore a project that references "
                         "it, or pass --download.")
    url = PACKAGE_PATH.format(feed=feed.rstrip("/"), v=version.lower())
    with urllib.request.urlopen(url, timeout=120) as response:  # noqa: S310 - https feed URL chosen by the caller
        return response.read(), url


def stage(nupkg: pathlib.Path | None, download: bool, dll: pathlib.Path | None = None, *, feed: str = NUGET_FEED,
          source_commit: str | None = None, build_command: str = DEFAULT_BUILD_COMMAND) -> dict[str, str]:
    """Stage Numerics.dll into _bin from a package (or a local build with ``dll``) and return its identity.

    The identity maps each name in ``IDENTITY_NAMES`` to the value written into the wheel's
    ``_version.py``. The version and assembly version must equal the tracked values, and a
    non-empty tracked ``NUMERICS_SHA256`` must equal the file's digest.
    """
    dist_version, numerics_version = _pin("__version__"), _pin("NUMERICS_VERSION")
    if dist_version.split("+")[0].split(".post")[0] != pep440(numerics_version):
        raise SystemExit(f"Distribution version {dist_version} does not match the Numerics version {numerics_version}.")
    if dll is not None:
        assembly, source = dll.read_bytes(), str(dll)
        informational = version_resource_string(assembly, "ProductVersion")
        version, _, metadata = informational.partition("+")
        commit = source_commit if source_commit is not None else metadata
        if not re.fullmatch(r"[0-9a-f]{40}", commit):
            raise SystemExit(f"{dll} names no source commit (informational version {informational!r}); "
                             "pass --source-commit.")
        built_from = (commit, build_command)
    else:
        payload, source = _package_bytes(numerics_version, nupkg, download, feed)
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            version = _nuspec_version(archive)
            assembly = archive.read(ASSEMBLY_ENTRY)
        built_from = ("", "")
    identity = {
        "NUMERICS_VERSION": version,
        "NUMERICS_ASSEMBLY_VERSION": version_resource_string(assembly, "Assembly Version"),
        "NUMERICS_SHA256": hashlib.sha256(assembly).hexdigest(),
        "NUMERICS_SOURCE_COMMIT": built_from[0],
        "NUMERICS_BUILD_COMMAND": built_from[1],
    }
    for name in ("NUMERICS_VERSION", "NUMERICS_ASSEMBLY_VERSION"):
        if identity[name] != _pin(name):
            raise SystemExit(f"Numerics.dll from {source} has {name} {identity[name]}; _version.py declares "
                             f"{_pin(name)}.")
    pinned = _pin("NUMERICS_SHA256")
    if pinned and identity["NUMERICS_SHA256"] != pinned:
        raise SystemExit(f"Numerics.dll from {source} has SHA-256 {identity['NUMERICS_SHA256']}; _version.py pins "
                         f"{pinned}.")
    if not (BIN / RUNTIME_CONFIG).is_file():
        raise SystemExit(f"{BIN / RUNTIME_CONFIG} is missing; it is tracked in the repository.")
    (BIN / "Numerics.dll").write_bytes(assembly)
    print(f"staged Numerics.dll {identity['NUMERICS_SHA256']} from {source}")
    return identity


def write_identity(text: str, identity: dict[str, str]) -> str:
    """Return _version.py source text with each identity constant set to its value (line endings kept)."""
    for name, value in identity.items():
        if re.search(r'["\\\r\n]', value):
            raise SystemExit(f"{name} value {value!r} cannot be written as a plain string literal.")
        text, count = re.subn(rf'^{name} = "[^"\r\n]*"', lambda _, n=name, v=value: f'{n} = "{v}"', text,
                              flags=re.MULTILINE)
        if count != 1:
            raise SystemExit(f"{VERSION_FILE.name} must define {name} once as a one-line string; found {count}.")
    return text


def build(out: pathlib.Path, identity: dict[str, str]) -> pathlib.Path:
    """Build the wheel with the identity written into its _version.py, verify its contents, and return its path."""
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pynumerics-wheel-") as temporary:
        source = pathlib.Path(temporary) / "pynumerics"
        source.mkdir()
        ignore = shutil.ignore_patterns("__pycache__", "*.egg-info", "*.pyc")
        for name in WHEEL_SOURCES:
            if (PACKAGE / name).is_dir():
                shutil.copytree(PACKAGE / name, source / name, ignore=ignore)
            else:
                shutil.copy2(PACKAGE / name, source / name)
        version_file = source / "src" / "pynumerics" / "_version.py"
        version_file.write_bytes(write_identity(VERSION_FILE.read_bytes().decode("utf-8"), identity).encode("utf-8"))
        subprocess.run([sys.executable, "-m", "build", "--wheel", "--outdir", str(out), str(source)], check=True)
    wheels = sorted(out.glob("pynumerics-*.whl"), key=lambda p: p.stat().st_mtime)
    if not wheels:
        raise SystemExit(f"no pynumerics wheel in {out}")
    wheel = wheels[-1]
    if not wheel.name.endswith("-py3-none-any.whl"):
        raise SystemExit(f"{wheel.name} is not a py3-none-any wheel.")
    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
        for name in ("Numerics.dll", RUNTIME_CONFIG, ROSLYN_LICENSE, *(file for _, file, _ in ROSLYN_FILES)):
            if f"pynumerics/_bin/{name}" not in names:
                raise SystemExit(f"{wheel.name} lacks pynumerics/_bin/{name}.")
        packed = hashlib.sha256(archive.read("pynumerics/_bin/Numerics.dll")).hexdigest()
        for _, file, pin_name in ROSLYN_FILES:
            roslyn = hashlib.sha256(archive.read(f"pynumerics/_bin/{file}")).hexdigest()
            if roslyn != _pin(pin_name):
                raise SystemExit(f"The wheel's {file} has SHA-256 {roslyn}; expected {_pin(pin_name)}.")
        written = archive.read(VERSION_ENTRY).decode("utf-8")
    if packed != identity["NUMERICS_SHA256"]:
        raise SystemExit(f"The wheel's Numerics.dll has SHA-256 {packed}; expected {identity['NUMERICS_SHA256']}.")
    for name, value in identity.items():
        if _pin(name, written) != value:
            raise SystemExit(f"The wheel's _version.py has {name} {_pin(name, written)!r}; expected {value!r}.")
    return wheel


def main(argv: list[str]) -> int:
    """Stage the assemblies, optionally build the wheel, and print the identity, the wheel path and its digest."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--nupkg", type=pathlib.Path, help="the RMC.Numerics .nupkg to take Numerics.dll from")
    parser.add_argument("--dll", type=pathlib.Path, help="a locally built Numerics.dll to stage instead of a package's")
    parser.add_argument("--download", action="store_true", help="download the packages if they are not cached")
    parser.add_argument("--feed", default=NUGET_FEED, help=f"flat container the RMC.Numerics package is downloaded "
                        f"from (default {NUGET_FEED})")
    parser.add_argument("--source-commit", help="the commit a --dll build was made from (default: from its "
                        "informational version)")
    parser.add_argument("--build-command", default=DEFAULT_BUILD_COMMAND,
                        help=f"how a --dll build was made (default {DEFAULT_BUILD_COMMAND!r})")
    parser.add_argument("--stage-only", action="store_true", help="stage the _bin assemblies and stop")
    parser.add_argument("--roslyn-only", action="store_true", help="stage only the Roslyn assemblies and stop")
    parser.add_argument("--out", type=pathlib.Path, default=PACKAGE / "dist", help="where the wheel is written")
    args = parser.parse_args(argv)
    if args.dll is not None and args.nupkg is not None:
        parser.error("--dll and --nupkg are exclusive")
    if args.dll is None and args.source_commit is not None:
        parser.error("--source-commit applies only to --dll")
    stage_roslyn(args.download)
    if args.roslyn_only:
        return 0
    identity = stage(args.nupkg, args.download, args.dll, feed=args.feed, source_commit=args.source_commit,
                     build_command=args.build_command)
    for name, value in identity.items():
        print(f"{name} = {value!r}")
    if args.stage_only:
        return 0
    wheel = build(args.out.resolve(), identity)
    print(wheel)
    print(f"{hashlib.sha256(wheel.read_bytes()).hexdigest()}  {wheel.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
