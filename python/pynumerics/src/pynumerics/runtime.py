"""Host the .NET runtime that runs RMC.Numerics inside Python.

This module is the single owner of the .NET runtime in a Python process that uses
RMC.Numerics. It loads CoreCLR through pythonnet, references ``Numerics.dll`` once,
and lets dependent packages add their own assemblies (for example ``Hydrologics.dll``
or ``RMC.BestFit.dll``) to the same runtime afterwards, so one process holds one
CoreCLR and one ``Numerics`` assembly.

``pythonnet.load("coreclr", ...)`` must run before any ``import clr``: on Windows a
bare ``import clr`` starts the legacy .NET Framework, that choice cannot be undone for
the life of the process, and the net10.0 assemblies cannot load there. The loader
refuses to continue when it detects that state.

The assembly directory is resolved once per process, and the first match wins:

1. the ``dll_dir`` argument of :func:`initialize`;
2. the ``PYNUMERICS_DLL_DIR`` environment variable;
3. the ``_bin`` directory packaged beside this module.

A ``dll_dir`` or ``PYNUMERICS_DLL_DIR`` that is not a directory is an error rather
than a reason to keep searching. The runtime configuration is always the packaged
``_bin/pynumerics.runtimeconfig.json`` unless the chosen directory holds its own.

The packaged configuration turns off quick JIT for methods with loops
(``System.Runtime.TieredCompilation.QuickJitForLoops``) and profile-guided tiering
(``System.Runtime.TieredPGO``). Most Numerics work from Python is one call that loops over an
array, and a process that keeps compiling new methods delays the promotion of unoptimized
first-tier code for seconds, so the array members ran up to about three times slower than their
optimized code. Compiling those methods optimized at first use removes that, at a start-up
cost of about twenty milliseconds. The settings change only when machine code is optimized, never the
floating-point operations it performs, so every value is unchanged; the ``DOTNET_TC_QuickJitForLoops``
and ``DOTNET_TieredPGO`` environment variables override them.

The .NET installation is ``DOTNET_ROOT`` when that variable is set, otherwise the one
clr_loader discovers (the default installation directory on Windows and macOS, or the
parent of the ``dotnet`` executable on ``PATH``).
"""

from __future__ import annotations

import os
import sys
import threading
from pathlib import Path
from typing import Any

ENV_DLL_DIR = "PYNUMERICS_DLL_DIR"
"""Environment variable naming a directory that holds ``Numerics.dll``."""

ASSEMBLY_FILE = "Numerics.dll"
RUNTIME_CONFIG_FILE = "pynumerics.runtimeconfig.json"

_PACKAGED_BIN = Path(__file__).resolve().parent / "_bin"

_runtime_dir: Path | None = None
_extra_assemblies: dict[str, Path] = {}
_LOCK = threading.RLock()  # pythonnet.load is not safe to race; concurrent first calls start it once


class RuntimeLoadError(RuntimeError):
    """Error raised when the .NET runtime or an assembly cannot be loaded.

    :func:`initialize` raises it when a runtime other than CoreCLR, such as the legacy
    .NET Framework, already occupies the process; when ``dll_dir`` or
    ``PYNUMERICS_DLL_DIR`` is not a directory; when no directory holding
    ``Numerics.dll`` is found; when no runtime configuration is available; and when a
    ``Numerics`` assembly from another file is already loaded. :func:`add_assembly`
    raises it when the file does not exist, would load a second ``Numerics`` assembly,
    or names an assembly already loaded from another file. A runtime that already occupies the process cannot be
    replaced, so the first case requires a new Python process.
    """


def initialize(dll_dir: str | os.PathLike[str] | None = None) -> Path:
    """Start the .NET runtime and load ``Numerics.dll``.

    Loads CoreCLR with the runtime configuration; when pythonnet already hosts CoreCLR
    in the process, uses that runtime instead of loading another. The call appends the
    assembly directory to ``sys.path`` and references ``Numerics.dll`` there by full path,
    so a copy elsewhere on ``sys.path`` is never loaded instead, after which namespaces
    such as ``Numerics.Distributions`` can be imported. After a successful call, later calls return immediately and
    ignore their arguments. Every wrapper in this package calls ``initialize`` itself;
    call it directly to choose the assembly directory or to surface loading problems
    early.

    Parameters
    ----------
    dll_dir : str or os.PathLike, optional
        Directory holding ``Numerics.dll``. None, the default, tries the
        ``PYNUMERICS_DLL_DIR`` environment variable and then the packaged ``_bin``
        directory.

    Returns
    -------
    pathlib.Path
        The directory ``Numerics.dll`` was loaded from.

    Raises
    ------
    pynumerics.RuntimeLoadError
        If a runtime other than CoreCLR already occupies the process, for example after
        a bare ``import clr`` on Windows; if ``dll_dir`` or ``PYNUMERICS_DLL_DIR`` is
        not a directory; if no directory holding ``Numerics.dll`` is found; if no
        runtime configuration is available; if no .NET installation is found; or if a
        ``Numerics`` assembly from another file is already loaded in the process.

    Notes
    -----
    Errors raised by pythonnet while starting the runtime, for example when no
    compatible .NET runtime is installed, propagate unchanged. The call is
    thread-safe: threads that call it together start the runtime once.

    Examples
    --------
    >>> import pynumerics
    >>> pynumerics.initialize()  # doctest: +SKIP
    >>> pynumerics.numerics_identity().version  # doctest: +SKIP
    '2.3.0.0'
    """
    global _runtime_dir
    if _runtime_dir is not None:
        return _runtime_dir
    with _LOCK:
        if _runtime_dir is not None:
            return _runtime_dir

        _guard_against_foreign_runtime()
        directory = _resolve_dll_dir(dll_dir)
        runtime_config = _resolve_runtime_config(directory)

        import pythonnet

        if pythonnet.get_runtime_info() is None:
            pythonnet.load("coreclr", runtime_config=str(runtime_config), dotnet_root=str(_resolve_dotnet_root()))
        _assert_coreclr(pythonnet.get_runtime_info())

        import clr  # noqa: PLC0415 - legal only after pythonnet.load above

        if str(directory) not in sys.path:
            sys.path.append(str(directory))
        _reference(clr, directory / ASSEMBLY_FILE)
        _runtime_dir = directory
        return directory


def add_assembly(path: str | os.PathLike[str]) -> Any:
    """Load another assembly into the runtime that hosts RMC.Numerics.

    Starts the runtime first if it is not running. The assembly's directory is
    appended to ``sys.path`` and the assembly is referenced by its full path, so its namespaces
    can be imported afterwards. Adding the same file again is a no-op. Dependent
    packages, such as the Hydrologics bindings, load their engine through this
    function so the process keeps one ``Numerics`` assembly: references from the added
    assembly to ``Numerics`` bind to the one already loaded.

    Parameters
    ----------
    path : str or os.PathLike
        Path to the assembly file, for example ``.../Hydrologics.dll``.

    Returns
    -------
    System.Reflection.Assembly
        The loaded assembly.

    Raises
    ------
    pynumerics.RuntimeLoadError
        If the file does not exist, if it is a ``Numerics`` assembly (only the one this
        package loads is allowed), if an assembly of the same name from another file is
        already loaded, or if the runtime cannot start.
    """
    file = Path(path).resolve()
    if not file.is_file():
        raise RuntimeLoadError(f"The assembly '{file}' does not exist.")
    if file.stem.lower() == "numerics":
        raise RuntimeLoadError(
            f"Refusing to load '{file}': one process loads one Numerics assembly, and "
            "pynumerics already owns it."
        )
    initialize()
    with _LOCK:
        import clr  # noqa: PLC0415 - the runtime is started above

        directory = str(file.parent)
        if directory not in sys.path:
            sys.path.append(directory)
        assembly = _reference(clr, file)
        _extra_assemblies[file.stem] = file
        return assembly


def loaded_assemblies() -> dict[str, Path]:
    """Return the assemblies added through :func:`add_assembly`, by name.

    Returns
    -------
    dict of str to pathlib.Path
        Assembly name to the file it was requested from.
    """
    with _LOCK:
        return dict(_extra_assemblies)


def runtime_dir() -> Path:
    """Return the directory ``Numerics.dll`` was loaded from.

    Starts the .NET runtime first if it is not already running.

    Returns
    -------
    pathlib.Path
        The assembly directory resolved by :func:`initialize`.

    Raises
    ------
    pynumerics.RuntimeLoadError
        If the runtime is not running yet and cannot be started.
    """
    return initialize()


def is_initialized() -> bool:
    """Return whether :func:`initialize` has completed in this process.

    Returns
    -------
    bool
        True after a successful :func:`initialize`.
    """
    return _runtime_dir is not None


def get_runtime_info() -> Any:
    """Return pythonnet's description of the .NET runtime hosting RMC.Numerics.

    Starts the .NET runtime first if it is not already running.

    Returns
    -------
    clr_loader.types.RuntimeInfo
        The runtime description, with ``kind``, ``version``, ``initialized``,
        ``shutdown``, and ``properties`` fields; ``kind`` is ``"CoreCLR"``.

    Raises
    ------
    pynumerics.RuntimeLoadError
        If the runtime is not running yet and cannot be started.
    """
    initialize()
    import pythonnet

    return pythonnet.get_runtime_info()


def _reference(clr: Any, file: Path) -> Any:
    """Reference an assembly by its full path and confirm the runtime bound that file.

    A reference by simple name loads the first file of that name on ``sys.path``, so a copy of the
    assembly in an earlier ``sys.path`` entry (the working directory, another package's binaries)
    would be loaded in place of the resolved one. The full path binds the resolved file; references
    by name made afterwards bind to the assembly already loaded.

    Raises
    ------
    pynumerics.RuntimeLoadError
        If an assembly with the same name from another file is already loaded in the process.
    """
    import System  # noqa: PLC0415 - the runtime is running

    try:
        assembly = clr.AddReference(str(file))
    except System.IO.FileLoadException as error:
        raise RuntimeLoadError(
            f"Cannot load '{file}': an assembly named '{file.stem}' from another file is already loaded in this "
            "process. Restart the process and call pynumerics.initialize() before any other CLR use."
        ) from error
    location = Path(str(assembly.Location))
    if os.path.normcase(str(location.resolve())) != os.path.normcase(str(file.resolve())):
        raise RuntimeLoadError(
            f"Cannot load '{file}': an assembly named '{file.stem}' is already loaded in this process from "
            f"'{location}'. Restart the process and call pynumerics.initialize() before any other CLR use."
        )
    return assembly


def _guard_against_foreign_runtime() -> None:
    """Fail loudly when ``clr`` was imported before this loader configured CoreCLR.

    A pre-existing ``clr`` module means other code initialized a CLR first; on Windows a
    bare ``import clr`` selects the legacy .NET Framework irreversibly, and the net10.0
    assemblies cannot load there.
    """
    if "clr" not in sys.modules:
        return

    import pythonnet

    info = pythonnet.get_runtime_info()
    if info is None or not _is_coreclr(info):
        raise RuntimeLoadError(
            "The 'clr' module was imported before pynumerics.initialize() and the "
            "process is not hosting CoreCLR. On Windows a bare 'import clr' initializes "
            "the legacy .NET Framework irreversibly; restart the process and call "
            "pynumerics.initialize() before any other CLR use."
        )


def _assert_coreclr(info: Any) -> None:
    """Raise when the loaded runtime is not CoreCLR (for example .NET Framework or Mono)."""
    if info is None or not _is_coreclr(info):
        kind = getattr(info, "kind", None)
        raise RuntimeLoadError(
            f"pythonnet is hosting '{kind}' but the net10.0 Numerics assembly requires "
            "CoreCLR. Restart the process so pynumerics can load CoreCLR before any "
            "other CLR use."
        )


def _is_coreclr(info: Any) -> bool:
    """Return True when a pythonnet RuntimeInfo describes CoreCLR."""
    return str(getattr(info, "kind", "")).lower().replace(" ", "") == "coreclr"


def _resolve_dll_dir(dll_dir: str | os.PathLike[str] | None) -> Path:
    """Resolve the assembly directory through the documented discovery order."""
    candidates: list[tuple[str, Path, bool]] = []
    if dll_dir is not None:
        candidates.append(("explicit dll_dir", Path(dll_dir), True))
    env_dir = os.environ.get(ENV_DLL_DIR)
    if env_dir:
        candidates.append((f"{ENV_DLL_DIR} environment variable", Path(env_dir), True))
    candidates.append(("packaged _bin", _PACKAGED_BIN, False))

    for label, candidate, must_exist in candidates:
        if candidate.is_dir() and (candidate / ASSEMBLY_FILE).is_file():
            return candidate.resolve()
        if must_exist:
            if not candidate.is_dir():
                raise RuntimeLoadError(f"The {label} '{candidate}' is not a directory.")
            raise RuntimeLoadError(f"The {label} '{candidate}' does not hold {ASSEMBLY_FILE}.")

    raise RuntimeLoadError(
        f"No directory holding {ASSEMBLY_FILE} found: pass dll_dir, set {ENV_DLL_DIR}, "
        "install the pynumerics wheel, or stage the assembly with "
        "python/pynumerics/scripts/build_wheel.py --stage-only."
    )


def _resolve_dotnet_root() -> Path:
    """Return the .NET installation root: ``DOTNET_ROOT`` when set, else clr_loader's discovery.

    clr_loader probes the default installation directory on Windows and macOS and the
    real parent of the ``dotnet`` executable on ``PATH`` elsewhere.
    """
    env_root = os.environ.get("DOTNET_ROOT")
    if env_root:
        root = Path(env_root)
        if not root.is_dir():
            raise RuntimeLoadError(f"DOTNET_ROOT '{root}' is not a directory.")
        return root
    from clr_loader.util.find import find_dotnet_root  # noqa: PLC0415 - pythonnet dependency

    try:
        root = Path(find_dotnet_root())
    except RuntimeError as error:
        raise RuntimeLoadError(
            "No .NET installation found. Install the .NET 10 runtime and put 'dotnet' on PATH, "
            "or set DOTNET_ROOT to the installation directory."
        ) from error
    if not root.is_dir():
        raise RuntimeLoadError(
            f"The .NET installation directory '{root}' does not exist. Install the .NET 10 runtime "
            "or set DOTNET_ROOT to the installation directory."
        )
    return root


def _resolve_runtime_config(directory: Path) -> Path:
    """Return the runtime configuration: the directory's own, else the packaged one."""
    for candidate in (directory / RUNTIME_CONFIG_FILE, _PACKAGED_BIN / RUNTIME_CONFIG_FILE):
        if candidate.is_file():
            return candidate
    raise RuntimeLoadError(f"No {RUNTIME_CONFIG_FILE} beside {ASSEMBLY_FILE} or in the packaged _bin directory.")
