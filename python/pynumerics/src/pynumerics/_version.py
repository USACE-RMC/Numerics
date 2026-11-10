"""Version of the pynumerics distribution and the identity of the Numerics assembly it ships.

Interim (until RMC.Numerics 2.3.0 is published): the wheel ships one locally built ``Numerics.dll``
from the release branch instead of a NuGet package assembly. It was built once, on Windows, from a
clean committed tree with::

    dotnet build Numerics/Numerics.csproj -c Release -f net10.0 -p:Version=2.3.0-dev.3
        -p:ContinuousIntegrationBuild=true

(.NET SDK 10.0.400 on its bundled runtime 10.0.11), and is staged with
``scripts/build_wheel.py --stage-only --dll PATH``. Its assembly version stays 2.2.0.0, so no
consumer pin moves in the interim; the SHA-256 and the informational version tell it apart. At the
release every value below moves to the published 2.3.0 package and ``NUMERICS_SOURCE_COMMIT``
becomes empty.
"""

__version__ = "2.3.0.dev3"
"""The distribution version; it equals the RMC.Numerics version of the bundled assembly (PEP 440 form)."""

NUMERICS_VERSION = "2.3.0-dev.3"
"""The RMC.Numerics version of the bundled ``Numerics.dll`` (a NuGet package version, or the local build's)."""

NUMERICS_ASSEMBLY_VERSION = "2.2.0.0"
"""The assembly version recorded in the bundled ``Numerics.dll`` (2.2.0.0 for the interim local build)."""

NUMERICS_SHA256 = "fb1a02dbd108dd8c135735253a66f60573cc6ad0914c31158550e6197ac9b91c"
"""SHA-256 of the bundled ``Numerics.dll`` (the interim local build of ``NUMERICS_SOURCE_COMMIT``)."""

NUMERICS_SOURCE_COMMIT = "ab0983a93e97a920379fdfdac4705713c90c374c"
"""The Numerics commit a locally built ``Numerics.dll`` was built from; empty for a NuGet package assembly."""

NUMERICS_BUILD_COMMAND = ("dotnet build Numerics/Numerics.csproj -c Release -f net10.0 -p:Version=2.3.0-dev.3 "
                          "-p:ContinuousIntegrationBuild=true "
                          "(Windows, .NET SDK 10.0.400 on its bundled runtime 10.0.11)")
"""How a locally built ``Numerics.dll`` was built (receipts record it); empty for a NuGet package assembly."""

ROSLYN_VERSION = "4.14.0"
"""The Microsoft.CodeAnalysis.Common and .CSharp package version bundled for ``callbacks.compile_csharp``."""

ROSLYN_TARGET_FRAMEWORK = "net9.0"
"""The ``lib`` folder of the Roslyn packages the bundled assemblies come from (they run on the .NET 10 runtime)."""

ROSLYN_CODEANALYSIS_SHA256 = "395736b2d0e02be77bedc3847b269da7a93bf8e49cd3431a4a8e68d2c52d5846"
"""SHA-256 of the bundled ``Microsoft.CodeAnalysis.dll`` (Microsoft.CodeAnalysis.Common ``ROSLYN_VERSION``)."""

ROSLYN_CSHARP_SHA256 = "53302acca07700e9bb10dcf21034435d6e0454c8eebf2481164ec13e520004c9"
"""SHA-256 of the bundled ``Microsoft.CodeAnalysis.CSharp.dll`` (Microsoft.CodeAnalysis.CSharp ``ROSLYN_VERSION``)."""
