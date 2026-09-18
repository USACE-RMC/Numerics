using System;
using System.Xml.Linq;

namespace Numerics.Functions
{
    /// <summary>
    /// Factory for creating <see cref="IUnivariateFunction"/> instances from
    /// <see cref="UnivariateFunctionType"/> enum values or from serialized <see cref="XElement"/>
    /// representations, and for mapping live instances back to their enum type.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// Mirrors the <see cref="LinkFunctionFactory"/> and
    /// <see cref="Numerics.Distributions.UnivariateDistributionFactory"/> patterns: a closed
    /// switch per surface, with serialization dispatch keyed on the element's local name (which
    /// by contract is the concrete class name). The enum deliberately stays off
    /// <see cref="IUnivariateFunction"/> — external implementors of the interface exist, and the
    /// net481 target has no default interface members, so the serialization surface lives on the
    /// concrete classes (<c>ToXElement()</c> instance methods with static <c>FromXElement</c>
    /// counterparts) and on this factory.
    /// </para>
    /// </remarks>
    public static class UnivariateFunctionFactory
    {
        /// <summary>
        /// Creates an <see cref="IUnivariateFunction"/> instance of the specified type with its
        /// default parameters.
        /// </summary>
        /// <param name="type">The univariate function type.</param>
        /// <returns>A new <see cref="IUnivariateFunction"/> instance.</returns>
        /// <exception cref="NotSupportedException">
        /// Thrown for <see cref="UnivariateFunctionType.Tabular"/>: a tabular function requires
        /// uncertain ordered paired data and cannot be created without parameters.
        /// </exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="type"/> is not a recognized function type.</exception>
        public static IUnivariateFunction CreateFunction(UnivariateFunctionType type)
        {
            switch (type)
            {
                case UnivariateFunctionType.Linear:
                    return new LinearFunction();
                case UnivariateFunctionType.Power:
                    return new PowerFunction();
                case UnivariateFunctionType.Tabular:
                    throw new NotSupportedException("Function type " + type + " requires uncertain ordered paired data and cannot be created without parameters.");
                case UnivariateFunctionType.SegmentedPower:
                    return new SegmentedPowerFunction();
                case UnivariateFunctionType.Composite:
                    throw new NotSupportedException("Function type " + type + " requires child functions and cannot be created without parameters.");
                default:
                    throw new ArgumentOutOfRangeException(nameof(type), type, "The function type is not defined.");
            }
        }

        /// <summary>
        /// Creates an <see cref="IUnivariateFunction"/> instance from a serialized
        /// <see cref="XElement"/>. The element name determines the function type.
        /// </summary>
        /// <param name="xElement">The XElement representing the function. The element name must match a known function class name.</param>
        /// <returns>A new <see cref="IUnivariateFunction"/> instance.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="NotSupportedException">Thrown when the element name does not correspond to a known function type.</exception>
        public static IUnivariateFunction CreateFromXElement(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));
            switch (xElement.Name.LocalName)
            {
                case nameof(LinearFunction):
                    return LinearFunction.FromXElement(xElement);
                case nameof(PowerFunction):
                    return PowerFunction.FromXElement(xElement);
                case nameof(TabularFunction):
                    return TabularFunction.FromXElement(xElement);
                case nameof(SegmentedPowerFunction):
                    return SegmentedPowerFunction.FromXElement(xElement);
                case nameof(CompositeFunction):
                    return CompositeFunction.FromXElement(xElement);
                default:
                    throw new NotSupportedException("Unknown function type: '" + xElement.Name.LocalName + "'.");
            }
        }

        /// <summary>
        /// Gets the <see cref="UnivariateFunctionType"/> of a live function instance.
        /// </summary>
        /// <param name="function">The function instance.</param>
        /// <returns>The function type.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="function"/> is null.</exception>
        /// <exception cref="NotSupportedException">Thrown when the instance is not one of the library's concrete function types.</exception>
        public static UnivariateFunctionType GetFunctionType(IUnivariateFunction function)
        {
            if (function == null) throw new ArgumentNullException(nameof(function));
            if (function is LinearFunction) return UnivariateFunctionType.Linear;
            if (function is PowerFunction) return UnivariateFunctionType.Power;
            if (function is TabularFunction) return UnivariateFunctionType.Tabular;
            if (function is SegmentedPowerFunction) return UnivariateFunctionType.SegmentedPower;
            if (function is CompositeFunction) return UnivariateFunctionType.Composite;
            throw new NotSupportedException("The function type '" + function.GetType().Name + "' is not a library function type.");
        }
    }
}
