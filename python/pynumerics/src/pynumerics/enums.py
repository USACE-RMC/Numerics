"""Python enums mirroring the Numerics enums the wrappers accept.

Every wrapper argument typed with one of these enums also accepts the member value or
the .NET member name as a string, matched case-insensitively with underscores ignored. The
enums themselves resolve the same spellings: ``Transform("normal_z")``, ``Transform("NORMAL_Z")``
and ``Transform("NormalZ")`` are all ``Transform.NORMAL_Z``. Enums added elsewhere in the package
derive from :class:`NetEnum` to get that lookup.
"""

from __future__ import annotations

import enum
from typing import Any

from pynumerics._util import fold_name


class NetEnum(enum.Enum):
    """Base of the pynumerics enums: values are .NET member names, looked up by folded spelling.

    A value that is not a member value exactly is matched against the member values and names
    with underscores removed and case ignored, the rule ``net_enum`` applies to .NET enums, so
    snake-case spellings kept in configuration files resolve.
    """

    @classmethod
    def _missing_(cls, value: Any) -> NetEnum | None:
        if not isinstance(value, str):
            return None
        key = fold_name(value)
        for member in cls:
            if fold_name(str(member.value)) == key or fold_name(member.name) == key:
                return member
        return None


class EstimationMethod(NetEnum):
    """Parameter estimation methods (``Numerics.Distributions.ParameterEstimationMethod``)."""

    MAXIMUM_LIKELIHOOD = "MaximumLikelihood"
    METHOD_OF_MOMENTS = "MethodOfMoments"
    METHOD_OF_LINEAR_MOMENTS = "MethodOfLinearMoments"
    METHOD_OF_PERCENTILES = "MethodOfPercentiles"


class SortOrder(NetEnum):
    """Sort orders (``Numerics.Data.SortOrder``)."""

    ASCENDING = "Ascending"
    DESCENDING = "Descending"
    NONE = "None"


class Transform(NetEnum):
    """Interpolation transforms (``Numerics.Data.Transform``)."""

    NONE = "None"
    LOGARITHMIC = "Logarithmic"
    NORMAL_Z = "NormalZ"


class Extrapolation(NetEnum):
    """Extrapolation sides (``Numerics.Data.ExtrapolationSides``)."""

    NONE = "None"
    BELOW = "Below"
    ABOVE = "Above"
    BOTH = "Both"


class PlottingPosition(NetEnum):
    """Plotting-position formulas (``PlottingPositions.PlottingPostionType``)."""

    WEIBULL = "Weibull"
    MEDIAN = "Median"
    BLOM = "Blom"
    CUNNANE = "Cunnane"
    GRINGORTEN = "Gringorten"
    HAZEN = "Hazen"
    KAPLAN_MEIER = "KaplanMeier"


class ImportanceDistribution(NetEnum):
    """Stratification importance distributions (``Stratify.ImportanceDistribution``)."""

    GUMBEL = "Gumbel"
    NORMAL = "Normal"
    UNIFORM = "Uniform"
    LOG10_UNIFORM = "Log10Uniform"


class WeightType(NetEnum):
    """Weight interpretations for weighted statistics (``Numerics.Data.Statistics.WeightType``)."""

    FREQUENCY = "Frequency"
    RELIABILITY = "Reliability"


class KernelType(NetEnum):
    """Kernel density kernels (``Numerics.Distributions.KernelDensity.KernelType``)."""

    EPANECHNIKOV = "Epanechnikov"
    GAUSSIAN = "Gaussian"
    TRIANGULAR = "Triangular"
    UNIFORM = "Uniform"


class AutocorrelationType(NetEnum):
    """Autocorrelation function kinds (``Numerics.Data.Statistics.Autocorrelation.Type``)."""

    CORRELATION = "Correlation"
    COVARIANCE = "Covariance"
    PARTIAL = "Partial"


__all__ = [
    "AutocorrelationType", "EstimationMethod", "Extrapolation", "ImportanceDistribution", "KernelType", "NetEnum",
    "PlottingPosition", "SortOrder", "Transform", "WeightType",
]
