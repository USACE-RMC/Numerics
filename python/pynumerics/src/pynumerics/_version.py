"""Version of the pynumerics distribution and the identity of the Numerics assembly it ships.

The wheel ships ``lib/net10.0/Numerics.dll`` of the RMC.Numerics NuGet package of the same version.
``scripts/build_wheel.py`` reads the package version, the assembly version and the SHA-256 from the
package and writes them into the wheel's copy of this file; it refuses a package whose version or
assembly version differs from the values below. In a source checkout ``NUMERICS_SHA256`` is empty,
because the digest of a published package is known only when it is built. ``NUMERICS_SOURCE_COMMIT``
and ``NUMERICS_BUILD_COMMAND`` are written only for a wheel built from a local ``Numerics.dll``
(``build_wheel.py --dll``, as in the pull request check); they are empty in a published wheel.
"""

__version__ = "2.3.0"
"""The distribution version; it equals the RMC.Numerics version of the bundled assembly (PEP 440 form)."""

NUMERICS_VERSION = "2.3.0"
"""The RMC.Numerics version of the bundled ``Numerics.dll`` (a NuGet package version, or the local build's)."""

NUMERICS_ASSEMBLY_VERSION = "2.3.0.0"
"""The assembly version recorded in the bundled ``Numerics.dll``."""

NUMERICS_SHA256 = ""
"""SHA-256 of the bundled ``Numerics.dll``, written when the wheel is built (empty in a source checkout)."""

NUMERICS_SOURCE_COMMIT = ""
"""The Numerics commit a locally built ``Numerics.dll`` was built from; empty for a NuGet package assembly."""

NUMERICS_BUILD_COMMAND = ""
"""How a locally built ``Numerics.dll`` was built (receipts record it); empty for a NuGet package assembly."""

ROSLYN_VERSION = "4.14.0"
"""The Microsoft.CodeAnalysis.Common and .CSharp package version bundled for ``callbacks.compile_csharp``."""

ROSLYN_TARGET_FRAMEWORK = "net9.0"
"""The ``lib`` folder of the Roslyn packages the bundled assemblies come from (they run on the .NET 10 runtime)."""

ROSLYN_CODEANALYSIS_SHA256 = "395736b2d0e02be77bedc3847b269da7a93bf8e49cd3431a4a8e68d2c52d5846"
"""SHA-256 of the bundled ``Microsoft.CodeAnalysis.dll`` (Microsoft.CodeAnalysis.Common ``ROSLYN_VERSION``)."""

ROSLYN_CSHARP_SHA256 = "53302acca07700e9bb10dcf21034435d6e0454c8eebf2481164ec13e520004c9"
"""SHA-256 of the bundled ``Microsoft.CodeAnalysis.CSharp.dll`` (Microsoft.CodeAnalysis.CSharp ``ROSLYN_VERSION``)."""
