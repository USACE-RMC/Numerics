"""Every Python code block of every documentation page runs on its own."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from test_docstrings import PUBLIC_MODULES

from pynumerics import callbacks

DOCS = Path(__file__).resolve().parents[1] / "docs"
BLOCK = re.compile(r"^```python\n(.*?)^```", re.S | re.M)
SNIPPETS = [
    (page, index, code)
    for page in sorted(DOCS.glob("*.md"))
    for index, code in enumerate(BLOCK.findall(page.read_text(encoding="utf-8")))
]


ROSLYN_MARKERS = ("=>", "compile_csharp", "distribution_log_likelihood")
"""Text that marks a block compiling C#: Python has no ``=>`` token, so ``=>`` is C# lambda source."""


def needs_roslyn(code: str) -> bool:
    """Whether a block compiles C# source without checking that Roslyn is staged."""
    return any(marker in code for marker in ROSLYN_MARKERS) and "roslyn_available()" not in code


def test_api_reference_lists_every_public_module() -> None:
    # Arrange
    listed = re.findall(r"^   (pynumerics[\w.]*)$", (DOCS / "api.rst").read_text(encoding="utf-8"), re.M)
    # Assert
    assert listed == list(PUBLIC_MODULES)


def test_every_page_has_snippets() -> None:
    # Assert
    pages_with_code = {page.name for page, _, _ in SNIPPETS}
    expected = {page.name for page in DOCS.glob("*.md")} - {"index.md"}
    assert pages_with_code == expected


@pytest.mark.parametrize(("page", "index", "code"), SNIPPETS, ids=[f"{p.stem}[{i}]" for p, i, _ in SNIPPETS])
def test_snippet_runs(page: Path, index: int, code: str) -> None:
    # Arrange
    if needs_roslyn(code) and not callbacks.roslyn_available():
        pytest.skip("Roslyn assemblies not staged")
    # Act / Assert
    exec(compile(code, f"{page.name}[{index}]", "exec"), {"__name__": "__docs__"})
