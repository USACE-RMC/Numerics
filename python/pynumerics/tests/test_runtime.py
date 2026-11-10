"""The loader: fresh-process start, the clr-first refusal, discovery and the load identity."""

from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

import pynumerics
from pynumerics import runtime
from pynumerics.identity import file_sha256

SRC = Path(pynumerics.__file__).resolve().parents[1]
"""The import root of the package under test (the source tree or an installed wheel)."""


def _run(code: str) -> subprocess.CompletedProcess[str]:
    """Run code in a fresh interpreter with the package source importable."""
    env_code = f"import sys; sys.path.insert(0, {str(SRC)!r})\n" + textwrap.dedent(code)
    return subprocess.run([sys.executable, "-c", env_code], capture_output=True, text=True, timeout=300)


def test_fresh_process_loads_coreclr_and_numerics() -> None:
    # Arrange / Act
    result = _run("""
        import sys
        import pynumerics
        assert "clr" not in sys.modules, "importing the package must not start the runtime"
        assert not pynumerics.is_initialized()
        pynumerics.initialize()
        print(pynumerics.get_runtime_info().kind)
        print(pynumerics.numerics_identity().version)
    """)
    # Assert
    assert result.returncode == 0, result.stderr
    kind, version = result.stdout.split()
    assert kind == "CoreCLR"
    assert version == pynumerics.NUMERICS_ASSEMBLY_VERSION


def test_bare_import_clr_first_is_refused() -> None:
    # Arrange / Act: a bare 'import clr' picks the default runtime (.NET Framework on Windows)
    result = _run("""
        try:
            import clr  # noqa: F401
        except Exception as exc:
            print("NOCLR", type(exc).__name__)
            raise SystemExit(0)
        import pynumerics
        try:
            pynumerics.initialize()
        except pynumerics.RuntimeLoadError as exc:
            print("REFUSED", exc)
        else:
            print("LOADED", pynumerics.get_runtime_info().kind)
    """)
    # Assert
    assert result.returncode == 0, result.stderr
    if result.stdout.startswith("NOCLR"):
        pytest.skip(f"no default CLR to import on this platform ({result.stdout.strip()})")
    assert result.stdout.startswith("REFUSED"), result.stdout
    assert "import clr" in result.stdout


def _copies(tmp_path: Path) -> tuple[Path, Path]:
    """Copy the loaded Numerics.dll into two directories, ``shadow`` and ``chosen``."""
    source = pynumerics.runtime_dir() / runtime.ASSEMBLY_FILE
    directories = (tmp_path / "shadow", tmp_path / "chosen")
    for directory in directories:
        directory.mkdir()
        (directory / runtime.ASSEMBLY_FILE).write_bytes(source.read_bytes())
    return directories


def test_resolved_directory_is_loaded_even_with_another_copy_earlier_on_sys_path(tmp_path: Path) -> None:
    # Arrange
    shadow, chosen = _copies(tmp_path)
    # Act
    result = _run(f"""
        import sys
        sys.path.insert(0, {str(shadow)!r})
        import pynumerics
        pynumerics.initialize({str(chosen)!r})
        print(pynumerics.numerics_identity().location)
    """)
    # Assert
    assert result.returncode == 0, result.stderr
    assert Path(result.stdout.strip()).resolve() == (chosen / runtime.ASSEMBLY_FILE).resolve()


def test_numerics_already_loaded_from_another_file_is_refused(tmp_path: Path) -> None:
    # Arrange
    shadow, chosen = _copies(tmp_path)
    # Act
    result = _run(f"""
        import pythonnet
        import pynumerics
        from pynumerics import runtime
        pythonnet.load("coreclr", runtime_config=str(runtime._PACKAGED_BIN / runtime.RUNTIME_CONFIG_FILE))
        import clr
        clr.AddReference({str(shadow / runtime.ASSEMBLY_FILE)!r})
        try:
            pynumerics.initialize({str(chosen)!r})
        except pynumerics.RuntimeLoadError as exc:
            print("REFUSED", exc)
    """)
    # Assert
    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("REFUSED"), result.stdout
    assert "already loaded" in result.stdout


def test_numerics_of_another_version_already_loaded_is_refused(tmp_path: Path) -> None:
    # Arrange: any packaged Numerics.dll of a different identity in the NuGet global packages folder
    loaded = file_sha256(pynumerics.runtime_dir() / runtime.ASSEMBLY_FILE)
    packages = Path.home() / ".nuget" / "packages" / "rmc.numerics"
    other = next((path for path in sorted(packages.glob("*/lib/net*/Numerics.dll"))
                  if file_sha256(path) != loaded), None)
    if other is None:
        pytest.skip("no other Numerics.dll in the NuGet global packages folder")
    _, chosen = _copies(tmp_path)
    # Act
    result = _run(f"""
        import pythonnet
        import pynumerics
        from pynumerics import runtime
        pythonnet.load("coreclr", runtime_config=str(runtime._PACKAGED_BIN / runtime.RUNTIME_CONFIG_FILE))
        import clr
        clr.AddReference({str(other)!r})
        try:
            pynumerics.initialize({str(chosen)!r})
        except pynumerics.RuntimeLoadError as exc:
            print("REFUSED", exc)
    """)
    # Assert
    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("REFUSED"), result.stdout
    assert "already loaded" in result.stdout


def test_explicit_directory_that_is_missing_is_an_error(tmp_path: Path) -> None:
    # Act / Assert
    with pytest.raises(pynumerics.RuntimeLoadError, match="not a directory"):
        runtime._resolve_dll_dir(tmp_path / "absent")
    with pytest.raises(pynumerics.RuntimeLoadError, match="does not hold"):
        runtime._resolve_dll_dir(tmp_path)


def test_environment_variable_directory_wins_over_packaged(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # Arrange
    (tmp_path / "Numerics.dll").write_bytes(b"")
    monkeypatch.setenv(runtime.ENV_DLL_DIR, str(tmp_path))
    # Act / Assert
    assert runtime._resolve_dll_dir(None) == tmp_path.resolve()
    assert runtime._resolve_runtime_config(tmp_path).name == runtime.RUNTIME_CONFIG_FILE


def test_dotnet_root_that_is_missing_is_an_error(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # Arrange
    monkeypatch.setenv("DOTNET_ROOT", str(tmp_path / "absent"))
    # Act / Assert
    with pytest.raises(pynumerics.RuntimeLoadError, match="DOTNET_ROOT"):
        runtime._resolve_dotnet_root()


def test_add_assembly_refuses_a_second_numerics_and_missing_files(tmp_path: Path) -> None:
    # Arrange
    other = tmp_path / "Numerics.dll"
    other.write_bytes(b"")
    # Act / Assert
    with pytest.raises(pynumerics.RuntimeLoadError, match="one Numerics assembly"):
        pynumerics.add_assembly(other)
    with pytest.raises(pynumerics.RuntimeLoadError, match="does not exist"):
        pynumerics.add_assembly(tmp_path / "Hydrologics.dll")


def test_initialize_is_memoized() -> None:
    # Act
    first = pynumerics.initialize()
    second = pynumerics.initialize("ignored after the first call")
    # Assert
    assert first == second == pynumerics.runtime_dir()


def test_load_identity_reports_version_and_file_digest() -> None:
    # Act
    identity = pynumerics.numerics_identity()
    # Assert
    assert identity.version == pynumerics.NUMERICS_ASSEMBLY_VERSION
    assert identity.informational_version.startswith(pynumerics.NUMERICS_VERSION)
    assert identity.location.resolve().parent == pynumerics.runtime_dir().resolve()
    assert identity.sha256 == file_sha256(identity.location)


def test_packaged_assembly_is_the_pinned_nuget_bytes() -> None:
    # Arrange
    identity = pynumerics.numerics_identity()
    if identity.location.resolve().parent != runtime._PACKAGED_BIN.resolve():
        pytest.skip("Numerics.dll was loaded from an override directory, not the packaged _bin")
    # Act / Assert
    assert identity.sha256 == pynumerics.NUMERICS_SHA256
    assert identity.is_packaged
    assert pynumerics.assert_numerics_identity() == identity
    with pytest.raises(pynumerics.RuntimeLoadError, match="SHA-256"):
        pynumerics.assert_numerics_identity(sha256="0" * 64)
    with pytest.raises(pynumerics.RuntimeLoadError, match="version"):
        pynumerics.assert_numerics_identity(version="0.0.0.0")


def test_distribution_version_equals_numerics_version() -> None:
    # Assert
    numerics = pynumerics.NUMERICS_VERSION
    base, _, dev = numerics.partition("-dev.")
    assert pynumerics.__version__.split(".post")[0] == (f"{base}.dev{dev}" if dev else numerics)
    if dev:
        # The interim local build keeps the 2.2.0.0 assembly version and names its source commit.
        assert pynumerics.NUMERICS_ASSEMBLY_VERSION == "2.2.0.0"
        assert len(pynumerics.NUMERICS_SOURCE_COMMIT) == 40
        assert numerics in pynumerics.NUMERICS_BUILD_COMMAND
    else:
        assert pynumerics.NUMERICS_ASSEMBLY_VERSION.startswith(numerics)
        assert pynumerics.NUMERICS_SOURCE_COMMIT == pynumerics.NUMERICS_BUILD_COMMAND == ""
