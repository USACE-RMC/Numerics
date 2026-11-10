#!/usr/bin/env bash
# Linux checks of RMC.Numerics and pynumerics, run inside the linux-test image.
#
# Inputs (environment):
#   SRC             the repository tree, mounted read-only (default /src)
#   CHECKS          the steps to run, any of: dotnet wheel pytest (default: all three)
#   NUMERICS_NUPKG  build the wheel from this RMC.Numerics package instead of the local net10.0 build
#   OUT             if this directory exists, the wheel is copied into it (default /out)
#   WHEEL_DIR       test the wheel in this directory instead of the one the wheel step builds
#
# The tree is copied to /work and its tracked files are checked out again with LF line endings, so the
# run sees what a Linux clone sees even when SRC is a Windows working tree with CRLF files.

set -euo pipefail

SRC=${SRC:-/src}
CHECKS=${CHECKS:-dotnet wheel pytest}
OUT=${OUT:-/out}
WORK=/work/numerics
PYTHON=python3.12
LOGS=/work/logs

has() { [[ " $CHECKS " == *" $1 "* ]]; }
step() { printf '\n==== %s ====\n' "$*"; }

step "Copy the tree from $SRC"
rm -rf "$WORK" "$LOGS"
mkdir -p "$WORK" "$LOGS"
tar -C "$SRC" \
    --exclude=./.git --exclude=bin --exclude=obj --exclude=TestResults --exclude=dist --exclude=build \
    --exclude=_build --exclude=generated --exclude=.venv --exclude=__pycache__ --exclude='*.egg-info' \
    --exclude=.pytest_cache --exclude=.ruff_cache --exclude=graphify-out --exclude='_bin/*.dll' \
    -cf - . | tar -C "$WORK" --no-same-owner -xf -
cd "$WORK"
# A CI checkout is owned by the runner's user, not root; keep git's ownership check from rejecting the copy.
git config --global --add safe.directory "$WORK"
# A repository with the GitHub remote, so SourceLink stamps the commit into the informational version.
git init -q
git config core.autocrlf input
git config core.safecrlf false
git remote add origin https://github.com/USACE-RMC/Numerics.git
git add -A
git -c user.name=linux-test -c user.email=linux-test@localhost commit -q -m "Linux check snapshot"
git ls-files -z | xargs -0 rm -f
git checkout -q -- .
echo "snapshot commit $(git rev-parse HEAD)"

if has dotnet; then
    step "dotnet build Test_Numerics (net10.0, Release)"
    dotnet build Test_Numerics/Test_Numerics.csproj -c Release -f net10.0 2>&1 | tee "$LOGS/dotnet-build.log"
    grep -Eq '^ +0 Warning\(s\)' "$LOGS/dotnet-build.log" || { echo "the build reported warnings"; exit 1; }

    step "dotnet test Test_Numerics (net10.0, Release)"
    dotnet test Test_Numerics/Test_Numerics.csproj -c Release -f net10.0 --no-build \
        --filter "TestCategory!=LiveServiceIntegration" 2>&1 | tee "$LOGS/dotnet-test.log"
fi

if has wheel; then
    cd "$WORK/python/pynumerics"
    $PYTHON -m venv /work/venv-build
    /work/venv-build/bin/python -m pip install -q build
    if [[ -n "${NUMERICS_NUPKG:-}" ]]; then
        step "Build the wheel from $NUMERICS_NUPKG"
        /work/venv-build/bin/python scripts/build_wheel.py --nupkg "$NUMERICS_NUPKG" --download --out /work/dist
    else
        step "Build the wheel from the local net10.0 Numerics.dll"
        dll="$WORK/Numerics/bin/Release/net10.0/Numerics.dll"
        if [[ ! -f "$dll" ]]; then
            dotnet build "$WORK/Numerics/Numerics.csproj" -c Release -f net10.0 2>&1 | tee "$LOGS/numerics-build.log"
        fi
        /work/venv-build/bin/python scripts/build_wheel.py --dll "$dll" --download --out /work/dist
    fi
    if [[ -d "$OUT" ]]; then
        cp /work/dist/pynumerics-*-py3-none-any.whl "$OUT"/
    fi
fi

if has pytest; then
    wheel_dir=${WHEEL_DIR:-/work/dist}
    step "Install the wheel from $wheel_dir into a clean environment"
    $PYTHON -m venv /work/venv-test
    /work/venv-test/bin/python -m pip install -q "$wheel_dir"/pynumerics-*-py3-none-any.whl pytest numpydoc ruff scipy \
        statsmodels scikit-learn shapely "arviz>=0.17,<1" "dataretrieval>=1.4" "wetterdienst>=0.184" "kiwis-pie>=0.10.5"

    step "ruff"
    (cd "$WORK/python/pynumerics" && /work/venv-test/bin/python -m ruff check src tests scripts)

    step "pytest against the installed wheel"
    rm -rf /work/pytest
    mkdir -p /work/pytest
    cp -r "$WORK/python/pynumerics/tests" "$WORK/python/pynumerics/docs" "$WORK/python/pynumerics/scripts" \
        "$WORK/python/pynumerics/pyproject.toml" /work/pytest/
    (cd /work/pytest && /work/venv-test/bin/python -m pytest -q -p no:cacheprovider -m "not benchmark" tests 2>&1 | tee "$LOGS/pytest.log")
fi

step "Done: $CHECKS"
