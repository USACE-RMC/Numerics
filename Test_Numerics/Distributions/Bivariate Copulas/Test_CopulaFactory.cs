using System;
using System.Xml.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions.Copulas;

namespace Distributions.BivariateCopulas
{
    /// <summary>
    /// Unit tests for the bivariate copula factory and the copula XML serialization surface.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_CopulaFactory
    {
        /// <summary>
        /// Pins the CopulaType member names and implicit values. The enumeration is serialized
        /// by name and has no explicit values, so its order is append-only contract: a new
        /// member may only ever be added after the last.
        /// </summary>
        [TestMethod]
        public void Test_CopulaType_ValuesArePinned()
        {
            var values = (CopulaType[])Enum.GetValues(typeof(CopulaType));
            CollectionAssert.AreEqual(new[]
            {
                "AliMikhailHaq", "Clayton", "Frank", "Gumbel", "Joe", "Normal", "StudentT", "Independence"
            }, Array.ConvertAll(values, v => v.ToString()));

            // The implicit values are contiguous from zero, so position and value agree.
            for (int i = 0; i < values.Length; i++)
                Assert.AreEqual(i, (int)values[i]);
        }

        /// <summary>
        /// Test that the factory creates every defined copula type and each instance reports
        /// the type it was created from.
        /// </summary>
        [TestMethod]
        public void Test_CreateCopula_EnumRoundTrip()
        {
            foreach (CopulaType type in Enum.GetValues(typeof(CopulaType)))
            {
                var copula = CopulaFactory.CreateCopula(type);
                Assert.IsNotNull(copula, $"{type} must be constructible.");
                Assert.AreEqual(type, copula.Type, $"{type} must round-trip through the factory.");
            }

            Assert.Throws<NotSupportedException>(() => CopulaFactory.CreateCopula((CopulaType)99));
        }

        /// <summary>
        /// Test that ToXElement followed by CreateCopula round-trips every family bitwise:
        /// parameters are written as pipe-delimited invariant-culture "G17" strings, which
        /// round-trip IEEE doubles exactly, so the reconstructed parameters must be
        /// bit-identical — asserted on deliberately irrational parameter values.
        /// </summary>
        [TestMethod]
        public void Test_ToXElement_RoundTrip_Bitwise()
        {
            var copulas = new BivariateCopula[]
            {
                new AMHCopula(0.57721566490153287),
                new ClaytonCopula(2.7182818284590452),
                new FrankCopula(-3.1415926535897931),
                new GumbelCopula(2.7182818284590452),
                new JoeCopula(3.1415926535897931),
                new NormalCopula(-0.57721566490153287),
                new StudentTCopula(0.31415926535897931, 7.3890560989306495),
                new IndependenceCopula()
            };

            foreach (var copula in copulas)
            {
                XElement element = copula.ToXElement();
                Assert.AreEqual("Copula", element.Name.LocalName);
                Assert.AreEqual(copula.Type.ToString(), element.Attribute("Type")?.Value);
                Assert.IsNotNull(element.Attribute("Parameters"));

                var restored = CopulaFactory.CreateCopula(element);
                Assert.AreEqual(copula.Type, restored.Type);
                Assert.IsTrue(restored.ParametersValid);
                CollectionAssert.AreEqual(copula.GetCopulaParameters, restored.GetCopulaParameters,
                    $"{copula.Type} parameters must round-trip bit-for-bit.");

                // Marginals are never part of the serialized form
                Assert.IsNull(restored.MarginalDistributionX);
                Assert.IsNull(restored.MarginalDistributionY);
            }
        }

        /// <summary>
        /// Test the zero-parameter round trip: the Independence copula writes an empty
        /// Parameters attribute, and the reader also accepts a missing Parameters attribute
        /// for a zero-parameter copula while rejecting it for a parameterized one.
        /// </summary>
        [TestMethod]
        public void Test_ZeroParameter_RoundTrip()
        {
            var copula = new IndependenceCopula();
            XElement element = copula.ToXElement();
            Assert.AreEqual("", element.Attribute("Parameters")?.Value);

            var restored = CopulaFactory.CreateCopula(element);
            Assert.AreEqual(CopulaType.Independence, restored.Type);
            Assert.IsEmpty(restored.GetCopulaParameters);

            var bare = new XElement("Copula");
            bare.SetAttributeValue("Type", nameof(CopulaType.Independence));
            Assert.AreEqual(CopulaType.Independence, CopulaFactory.CreateCopula(bare).Type);

            var bareClayton = new XElement("Copula");
            bareClayton.SetAttributeValue("Type", nameof(CopulaType.Clayton));
            Assert.Throws<ArgumentException>(() => CopulaFactory.CreateCopula(bareClayton));
        }

        /// <summary>
        /// Test that malformed serialized forms throw: a null element, a missing or unknown
        /// or undefined type, a parameter count that does not match the copula type, an
        /// unparseable or non-finite parameter, and parameters outside the copula's valid
        /// domain.
        /// </summary>
        [TestMethod]
        public void Test_CreateCopula_Malformed_Throws()
        {
            Assert.Throws<ArgumentNullException>(() => CopulaFactory.CreateCopula((XElement)null!));

            Assert.Throws<ArgumentException>(() => CopulaFactory.CreateCopula(new XElement("Copula")));

            var bogusType = new XElement("Copula");
            bogusType.SetAttributeValue("Type", "Bogus");
            Assert.Throws<ArgumentException>(() => CopulaFactory.CreateCopula(bogusType));

            var undefinedType = new XElement("Copula");
            undefinedType.SetAttributeValue("Type", "99");
            Assert.Throws<ArgumentException>(() => CopulaFactory.CreateCopula(undefinedType));

            var wrongCount = new XElement("Copula");
            wrongCount.SetAttributeValue("Type", nameof(CopulaType.Clayton));
            wrongCount.SetAttributeValue("Parameters", "1|2");
            Assert.Throws<ArgumentException>(() => CopulaFactory.CreateCopula(wrongCount));

            var unparseable = new XElement("Copula");
            unparseable.SetAttributeValue("Type", nameof(CopulaType.Clayton));
            unparseable.SetAttributeValue("Parameters", "abc");
            Assert.Throws<ArgumentException>(() => CopulaFactory.CreateCopula(unparseable));

            var nonFinite = new XElement("Copula");
            nonFinite.SetAttributeValue("Type", nameof(CopulaType.Clayton));
            nonFinite.SetAttributeValue("Parameters", "NaN");
            Assert.Throws<ArgumentException>(() => CopulaFactory.CreateCopula(nonFinite));

            // Gumbel requires theta >= 1, so a structurally well-formed but out-of-domain
            // parameter must be rejected by the post-assignment validity check.
            var outOfDomain = new XElement("Copula");
            outOfDomain.SetAttributeValue("Type", nameof(CopulaType.Gumbel));
            outOfDomain.SetAttributeValue("Parameters", "0.5");
            Assert.Throws<ArgumentException>(() => CopulaFactory.CreateCopula(outOfDomain));
        }
    }
}
