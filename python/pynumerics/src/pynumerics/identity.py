"""Report which ``Numerics`` assembly the process loaded.

One process loads one ``Numerics`` assembly. Dependent packages assert that the loaded
assembly is the one they were built against: its version, and the SHA-256 of the file
it was loaded from.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

from pynumerics._version import NUMERICS_ASSEMBLY_VERSION, NUMERICS_SHA256
from pynumerics.runtime import RuntimeLoadError, initialize


@dataclass(frozen=True)
class NumericsIdentity:
    """Identity of the loaded ``Numerics`` assembly.

    Parameters
    ----------
    version : str
        The assembly version, for example ``"2.3.0.0"``.
    informational_version : str
        The informational (product) version attribute, or an empty string.
    location : pathlib.Path
        The file the runtime loaded the assembly from.
    sha256 : str
        Lower-case hexadecimal SHA-256 of that file.
    """

    version: str
    informational_version: str
    location: Path
    sha256: str

    @property
    def is_packaged(self) -> bool:
        """Whether the loaded file is byte-identical to the pinned NuGet assembly."""
        return self.sha256 == NUMERICS_SHA256 and self.version == NUMERICS_ASSEMBLY_VERSION


def file_sha256(path: str | Path) -> str:
    """Return the lower-case hexadecimal SHA-256 of a file.

    Parameters
    ----------
    path : str or pathlib.Path
        The file to hash.

    Returns
    -------
    str
        The 64-character digest.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def numerics_identity() -> NumericsIdentity:
    """Return the identity of the ``Numerics`` assembly loaded in this process.

    Starts the runtime first if it is not running. The assembly is located through the
    runtime (``typeof(Numerics.Tools).Assembly``), not through the discovery rules, so
    the report describes what is actually loaded.

    Returns
    -------
    NumericsIdentity
        Version, informational version, file location and SHA-256.

    Raises
    ------
    pynumerics.RuntimeLoadError
        If the runtime cannot start.
    """
    initialize()
    import clr  # noqa: PLC0415
    import Numerics  # noqa: PLC0415
    import System  # noqa: PLC0415 - resolvable only after initialize
    from System.Reflection import AssemblyInformationalVersionAttribute  # noqa: PLC0415

    assembly = System.Reflection.Assembly.GetAssembly(clr.GetClrType(Numerics.Tools))
    name = assembly.GetName()
    info = ""
    attribute = System.Attribute.GetCustomAttribute(assembly, clr.GetClrType(AssemblyInformationalVersionAttribute))
    if attribute is not None:
        info = str(getattr(attribute, "__implementation__", attribute).InformationalVersion)
    location = Path(str(assembly.Location))
    return NumericsIdentity(
        version=str(name.Version),
        informational_version=info,
        location=location,
        sha256=file_sha256(location),
    )


def assert_numerics_identity(version: str | None = None, sha256: str | None = None) -> NumericsIdentity:
    """Assert that the loaded ``Numerics`` assembly has the given version and digest.

    Parameters
    ----------
    version : str, optional
        Expected assembly version, for example the ``Numerics`` reference recorded in
        ``Hydrologics.dll``. None, the default, uses the version this package pins.
    sha256 : str, optional
        Expected SHA-256 of the loaded file. None, the default, uses the digest this
        package pins, which the wheel build writes; a source checkout pins none.

    Returns
    -------
    NumericsIdentity
        The identity that matched.

    Raises
    ------
    pynumerics.RuntimeLoadError
        If the loaded assembly's version or digest differs from the expectation, or if no
        digest is given and this package pins none (a source checkout).
    """
    identity = numerics_identity()
    expected_version = NUMERICS_ASSEMBLY_VERSION if version is None else version
    expected_sha = NUMERICS_SHA256 if sha256 is None else sha256.lower()
    if not expected_sha:
        raise RuntimeLoadError(
            "This pynumerics source checkout pins no Numerics.dll digest: NUMERICS_SHA256 is written into "
            "the wheel when it is built. Pass sha256, or install the wheel."
        )
    if identity.version != expected_version:
        raise RuntimeLoadError(
            f"The loaded Numerics assembly is version {identity.version} ({identity.location}); "
            f"expected {expected_version}."
        )
    if identity.sha256 != expected_sha:
        raise RuntimeLoadError(
            f"The loaded Numerics assembly {identity.location} has SHA-256 {identity.sha256}; "
            f"expected {expected_sha}."
        )
    return identity
