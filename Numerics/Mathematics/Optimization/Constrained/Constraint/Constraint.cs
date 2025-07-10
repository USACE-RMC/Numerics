/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
