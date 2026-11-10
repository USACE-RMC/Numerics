# Linux test image for RMC.Numerics and pynumerics: the .NET 10 SDK on Ubuntu 24.04 with Python 3.12.
#
# Build from the repository root:
#   docker build -t numerics-linux-test .github/docker/linux-test
# Run every check against a read-only mount of the working tree (see run-checks.sh):
#   docker run --rm -v "$PWD:/src:ro" numerics-linux-test
#
# The image holds only the toolchains. The checks copy the tree into the container, so build
# output never lands in the mounted checkout, and each Python environment is created fresh.

FROM mcr.microsoft.com/dotnet/sdk:10.0-noble

ENV DEBIAN_FRONTEND=noninteractive \
    DOTNET_CLI_TELEMETRY_OPTOUT=1 \
    DOTNET_NOLOGO=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends python3.12 python3.12-venv \
    && rm -rf /var/lib/apt/lists/*

COPY run-checks.sh identity-check.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/run-checks.sh /usr/local/bin/identity-check.sh

WORKDIR /work
CMD ["run-checks.sh"]
