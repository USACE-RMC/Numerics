using System;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// A class for base functionality for constraints.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class Constraint : IConstraint
    {
        /// <summary>
        /// Construct a new constraint. 
        /// </summary>
        /// <param name="constraintFunction">The constraint function.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="value">The value on the right hand side of the constraint equation. </param>
        /// <param name="type">The constraint type. </param>
        /// <param name="tolerance">The violation tolerance for the constraint. Default = 1E-8.</param>
        public Constraint(Func<double[], double> constraintFunction, int numberOfParameters, double value, ConstraintType type, double tolerance = 1E-8)
        {
            Function = constraintFunction;
            NumberOfParameters = numberOfParameters;
            Value = value;
            Type = type;
            Tolerance = tolerance;
        }

        /// <summary>
        /// The constraint type. 
        /// </summary>
        public ConstraintType Type { get; }

        /// <summary>
        /// The number of parameters to evaluate in the function.
        /// </summary>
        public int NumberOfParameters { get; }

        /// <summary>
        /// Calculates the left hand side of the constraint.
        /// </summary>
        public Func<double[], double> Function { get; }

        /// <summary>
        /// The value on the right hand side of the constraint equation. 
        /// </summary>
        public double Value { get; }

        /// <summary>
        /// The violation tolerance for the constraint. Default = 1E-8. 
        /// </summary>
        public double Tolerance { get; }

    }
}
