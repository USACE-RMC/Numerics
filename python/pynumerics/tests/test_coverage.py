"""Every public Numerics type and member maps to a pynumerics wrapper, an exclusion or a tracked gap.

The map is ``pynumerics._coverage.COVERAGE``. A type with no entry, an entry for a type that no
longer exists, a wrapper path that does not resolve, or a declared member of a wrapped type that is
not accounted for fails. Remaining ``Todo`` entries make the completeness test an expected failure
(not strict) whose message counts them per pynumerics module.
"""

from __future__ import annotations

import dataclasses
import importlib
from collections import Counter
from typing import Any

import pytest

from pynumerics import _coverage
from pynumerics._coverage import COVERAGE, MODULES, Excluded, Todo, Wrapped, net_members
from pynumerics._util import namespace


def _exported() -> dict[str, Any]:
    """Return the public types of the loaded Numerics assembly, by full name."""
    namespace("Numerics")
    import System  # noqa: PLC0415

    assembly = next(a for a in System.AppDomain.CurrentDomain.GetAssemblies() if a.GetName().Name == "Numerics")
    return {str(t.FullName): t for t in assembly.GetExportedTypes()}


def _attribute(obj: Any, name: str) -> Any:
    """Return an attribute, or the field of a dataclass whose fields have no class default."""
    try:
        return getattr(obj, name)
    except AttributeError:
        if dataclasses.is_dataclass(obj):
            field = next((f for f in dataclasses.fields(obj) if f.name == name), None)
            if field is not None:
                return field
        raise


def _resolve(path: str) -> Any:
    """Import the longest module prefix of a dotted path and walk the remaining attributes."""
    parts = path.split(".")
    for index in range(len(parts), 0, -1):
        try:
            obj = importlib.import_module(".".join(parts[:index]))
        except ImportError:
            continue
        for part in parts[index:]:
            obj = _attribute(obj, part)
        return obj
    raise LookupError(path)


def _is_member_type(clr_type: Any) -> bool:
    """Return whether members of the type are checked (classes and structs, not enums, interfaces or delegates)."""
    import clr  # noqa: PLC0415
    import System  # noqa: PLC0415

    return not (clr_type.IsEnum or clr_type.IsInterface
                or clr_type.IsSubclassOf(clr.GetClrType(System.MulticastDelegate)))


def _member_gaps() -> Counter:
    """Return the number of declared members of wrapped types still marked Todo, per pynumerics module."""
    exported, gaps = _exported(), Counter()
    for name, entry in COVERAGE.items():
        if not isinstance(entry, Wrapped) or name not in exported or not _is_member_type(exported[name]):
            continue
        for member in net_members(exported[name]):
            state = entry.members.get(member, entry.rest)
            if isinstance(state, Todo):
                gaps[state.module] += 1
    return gaps


def test_every_public_type_has_an_entry_and_every_entry_a_type() -> None:
    # Arrange
    exported = _exported()
    # Act
    missing = sorted(set(exported) - set(COVERAGE))
    stale = sorted(set(COVERAGE) - set(exported))
    # Assert
    assert not missing, f"{len(missing)} public Numerics types have no coverage entry: {missing}"
    assert not stale, f"{len(stale)} coverage entries name no public Numerics type: {stale}"


def test_entries_are_well_formed() -> None:
    # Act / Assert
    for name, entry in COVERAGE.items():
        assert isinstance(entry, (Wrapped, Excluded, Todo)), name
        if isinstance(entry, Todo):
            assert entry.module in MODULES, f"{name}: unknown module {entry.module}"
        if isinstance(entry, Excluded):
            assert entry.reason.strip(), f"{name}: an exclusion needs a reason"
        if isinstance(entry, Wrapped):
            assert entry.rest is None or isinstance(entry.rest, (Todo, Excluded)), name
            if isinstance(entry.rest, Todo):
                assert entry.rest.module in MODULES, name


def test_wrapped_targets_and_members_resolve() -> None:
    # Arrange
    unresolved = []
    # Act
    for name, entry in COVERAGE.items():
        if not isinstance(entry, Wrapped):
            continue
        try:
            target = _resolve(entry.target)
        except (LookupError, AttributeError):
            unresolved.append(f"{name} -> {entry.target}")
            continue
        for member, state in entry.members.items():
            if not isinstance(state, str):
                continue
            try:
                _resolve(state) if state.startswith("pynumerics.") else _attribute(target, state)
            except (LookupError, AttributeError):
                unresolved.append(f"{name}.{member} -> {entry.target}.{state}")
    # Assert
    assert not unresolved, "unresolved wrapper paths:\n" + "\n".join(unresolved)


def test_wrapped_types_account_for_every_declared_member() -> None:
    # Arrange
    exported = _exported()
    unaccounted, unknown = [], []
    # Act
    for name, entry in COVERAGE.items():
        if not isinstance(entry, Wrapped) or not _is_member_type(exported[name]):
            continue
        declared = set(net_members(exported[name]))
        unknown += [f"{name}.{m}" for m in entry.members if m not in declared]
        if entry.rest is None:
            unaccounted += [f"{name}.{m}" for m in sorted(declared - set(entry.members))]
    # Assert
    assert not unknown, "members listed that the type does not declare:\n" + "\n".join(unknown)
    assert not unaccounted, "declared members with no entry:\n" + "\n".join(unaccounted)


TYPE_TODOS = _coverage.todo_counts()


@pytest.mark.xfail(condition=sum(TYPE_TODOS.values()) > 0, strict=False,
                   reason=f"{sum(TYPE_TODOS.values())} public types still to wrap: "
                          + ", ".join(f"{module} {count}" for module, count in sorted(TYPE_TODOS.items())))
def test_coverage_is_complete() -> None:
    # Act
    members = _member_gaps()
    # Assert
    types = ", ".join(f"{module} {count}" for module, count in sorted(TYPE_TODOS.items()))
    gaps = ", ".join(f"{module} {count}" for module, count in sorted(members.items()))
    assert not TYPE_TODOS and not members, f"Todo types: {types or 'none'}; Todo members of wrapped types: {gaps}"
