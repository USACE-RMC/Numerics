#!/usr/bin/env bash
# Installs a pynumerics wheel into a clean environment and checks the identity of the loaded Numerics.dll.
#
# Usage (a fresh container, the wheel's directory mounted read-only):
#   docker run --rm -v "$PWD/dist:/dist:ro" numerics-linux-test identity-check.sh /dist
#
# The installed package must load its own packaged Numerics.dll on CoreCLR, and that file must match the
# identity the wheel build wrote into _version.py (assembly version, informational version, SHA-256).

set -euo pipefail

DIST=${1:-/dist}
python3.12 -m venv /tmp/identity-venv
/tmp/identity-venv/bin/python -m pip install -q "$DIST"/pynumerics-*-py3-none-any.whl
cd /tmp
/tmp/identity-venv/bin/python -I - <<'EOF'
import pynumerics
from pynumerics import runtime

identity = pynumerics.assert_numerics_identity()
info = pynumerics.get_runtime_info()
assert str(info.kind) == "CoreCLR", info
assert identity.location.parent.resolve() == runtime._PACKAGED_BIN.resolve(), identity.location
assert identity.informational_version.split("+")[0] == pynumerics.NUMERICS_VERSION, identity
assert identity.is_packaged
print("pynumerics", pynumerics.__version__)
print("runtime", info.kind)
print("NUMERICS_VERSION", pynumerics.NUMERICS_VERSION)
print("NUMERICS_ASSEMBLY_VERSION", pynumerics.NUMERICS_ASSEMBLY_VERSION, "loaded", identity.version)
print("NUMERICS_SHA256", pynumerics.NUMERICS_SHA256)
print("loaded sha256   ", identity.sha256)
print("informational", identity.informational_version)
print("NUMERICS_SOURCE_COMMIT", repr(pynumerics.NUMERICS_SOURCE_COMMIT))
print("location", identity.location)
print("identity check passed")
EOF
