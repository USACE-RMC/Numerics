"""scripts/build_wheel.py: the identity read from the assembly and package, and the identity written into the wheel."""

from __future__ import annotations

import importlib.util
import io
import zipfile
from pathlib import Path
from types import ModuleType

import pytest

import pynumerics

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_wheel.py"


@pytest.fixture(scope="module")
def build_wheel() -> ModuleType:
    """Import scripts/build_wheel.py as a module."""
    if not SCRIPT.is_file():
        pytest.skip(f"{SCRIPT} is not next to the tests")
    spec = importlib.util.spec_from_file_location("build_wheel", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_version_resource_matches_the_loaded_assembly(build_wheel: ModuleType) -> None:
    # Arrange: the oracle is the CoreCLR reflection identity of the same file
    identity = pynumerics.numerics_identity()
    assembly = identity.location.read_bytes()
    # Act
    assembly_version = build_wheel.version_resource_string(assembly, "Assembly Version")
    informational = build_wheel.version_resource_string(assembly, "ProductVersion")
    # Assert
    assert assembly_version == identity.version
    assert informational == identity.informational_version


def test_version_resource_missing_key_is_refused(build_wheel: ModuleType) -> None:
    # Act / Assert
    with pytest.raises(SystemExit, match="no 'Assembly Version'"):
        build_wheel.version_resource_string(b"MZ" + bytes(64), "Assembly Version")


def test_nuspec_version_is_read_from_the_package_root(build_wheel: ModuleType) -> None:
    # Arrange: a minimal package with a namespaced nuspec and a nested decoy
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, "w") as archive:
        archive.writestr("RMC.Numerics.nuspec", '﻿<?xml version="1.0" encoding="utf-8"?>'
                         '<package xmlns="http://schemas.microsoft.com/packaging/2013/05/nuspec.xsd">'
                         "<metadata><id>RMC.Numerics</id><version> 2.3.0 </version></metadata></package>")
        archive.writestr("lib/other.nuspec", "<package><metadata><version>9.9.9</version></metadata></package>")
    # Act
    with zipfile.ZipFile(payload) as archive:
        version = build_wheel._nuspec_version(archive)
    # Assert
    assert version == "2.3.0"


def test_write_identity_replaces_each_constant_and_keeps_line_endings(build_wheel: ModuleType) -> None:
    # Arrange
    source = ('"""Doc."""\r\n\r\nNUMERICS_VERSION = "2.3.0"\r\n"""Its doc."""\r\n'
              'NUMERICS_SHA256 = ""\r\nNUMERICS_SHA256_NOTE = "kept"\r\n')
    # Act
    written = build_wheel.write_identity(source, {"NUMERICS_VERSION": "2.3.1", "NUMERICS_SHA256": "ab" * 32})
    # Assert
    assert written == ('"""Doc."""\r\n\r\nNUMERICS_VERSION = "2.3.1"\r\n"""Its doc."""\r\n'
                       f'NUMERICS_SHA256 = "{"ab" * 32}"\r\nNUMERICS_SHA256_NOTE = "kept"\r\n')
    assert build_wheel._pin("NUMERICS_SHA256", written) == "ab" * 32


def test_write_identity_refuses_a_missing_constant_or_an_unsafe_value(build_wheel: ModuleType) -> None:
    # Act / Assert
    with pytest.raises(SystemExit, match="NUMERICS_BUILD_COMMAND once"):
        build_wheel.write_identity('NUMERICS_VERSION = "2.3.0"\n', {"NUMERICS_BUILD_COMMAND": "dotnet build"})
    with pytest.raises(SystemExit, match="plain string literal"):
        build_wheel.write_identity('NUMERICS_VERSION = "2.3.0"\n', {"NUMERICS_VERSION": 'x"y'})


def test_tracked_version_file_defines_every_identity_constant_once(build_wheel: ModuleType) -> None:
    # Arrange
    if not build_wheel.VERSION_FILE.is_file():
        pytest.skip("no source checkout next to the tests")
    text = build_wheel.VERSION_FILE.read_text(encoding="utf-8")
    values = {name: build_wheel._pin(name, text) for name in build_wheel.IDENTITY_NAMES}
    # Act
    written = build_wheel.write_identity(text, values)
    # Assert: writing the tracked values back is the identity transformation
    assert written == text
