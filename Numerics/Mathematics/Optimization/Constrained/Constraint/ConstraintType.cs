using System;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// Enumeration of constraint types. 
    /// </summary>
    ///  <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public enum ConstraintType
    {
        /// <summary>
        /// Equality constraint, h(x) = 0
        /// </summary>
        EqualTo,

        /// <summary>
        /// Inequality constraint for greater than or equal to, h(x) &gt;= 0
        /// </summary>
        GreaterThanOrEqualTo,

        /// <summary>
        /// Inequality constraint for lesser than or equal to, h(x) &lt;= 0
        /// </summary>
        LesserThanOrEqualTo
    }
}
