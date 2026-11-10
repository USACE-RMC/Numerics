"""Every public object carries a NumPy-style docstring numpydoc accepts."""

from __future__ import annotations

import enum
import importlib
import inspect

import pytest
from numpydoc import validate as numpydoc_validate

PUBLIC_MODULES = (
    "pynumerics",
    "pynumerics.runtime",
    "pynumerics.marshal",
    "pynumerics.identity",
    "pynumerics.distributions",
    "pynumerics.multivariate",
    "pynumerics.copulas",
    "pynumerics.regional",
    "pynumerics.data",
    "pynumerics.statistics",
    "pynumerics.sampling",
    "pynumerics.mcmc",
    "pynumerics.mathematics",
    "pynumerics.bootstrap",
    "pynumerics.integration",
    "pynumerics.differentiation",
    "pynumerics.roots",
    "pynumerics.ode",
    "pynumerics.machine_learning",
    "pynumerics.hypothesis_tests",
    "pynumerics.goodness_of_fit",
    "pynumerics.probability",
    "pynumerics.optimization",
    "pynumerics.graphs",
    "pynumerics.special",
    "pynumerics.linalg",
    "pynumerics.functions",
    "pynumerics.timeseries",
    "pynumerics.download",
    "pynumerics.errors",
    "pynumerics.callbacks",
    "pynumerics.enums",
)

ENFORCED_CHECKS = frozenset({
    "GL02", "GL03", "GL05", "GL06", "GL07", "GL08", "GL09", "GL10",
    "SS01", "SS02", "SS03", "SS04", "SS05", "SS06",
    "PR01", "PR02", "PR03", "PR04", "PR05", "PR06", "PR07", "PR08", "PR09", "PR10",
    "RT01", "RT02", "RT03", "RT04", "RT05",
    "YD01",
})
"""numpydoc checks that fail the build (GL01, ES01, SA01 and EX01 are not required)."""


def _public_names() -> list[str]:
    """Return the dotted names of every public function, class and class member."""
    names: set[str] = set()
    for module_name in PUBLIC_MODULES:
        module = importlib.import_module(module_name)
        exported = getattr(module, "__all__", None)
        if exported is None:
            exported = [n for n, obj in vars(module).items()
                        if not n.startswith("_") and (inspect.isfunction(obj) or inspect.isclass(obj))
                        and getattr(obj, "__module__", "") == module_name]
        for name in exported:
            obj = getattr(module, name)
            if not (inspect.isfunction(obj) or inspect.isclass(obj)):
                continue
            if inspect.isclass(obj) and issubclass(obj, enum.Enum):
                continue
            home = f"{obj.__module__}.{obj.__qualname__}"
            names.add(home)
            if inspect.isclass(obj) and not issubclass(obj, BaseException):
                for member, value in vars(obj).items():
                    if member.startswith("_") or member in {"name", "value"}:
                        continue
                    if inspect.isfunction(value) or isinstance(value, (property, classmethod, staticmethod)):
                        names.add(f"{home}.{member}")
    return sorted(names)


@pytest.mark.parametrize("name", _public_names())
def test_public_docstring(name: str) -> None:
    # Act
    report = numpydoc_validate.validate(name)
    # Assert
    errors = [f"{code}: {message}" for code, message in report["errors"] if code in ENFORCED_CHECKS]
    assert not errors, f"{name}\n" + "\n".join(errors)
