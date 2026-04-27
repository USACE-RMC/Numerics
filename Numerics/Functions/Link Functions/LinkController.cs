using System;
using System.Linq;
using System.Xml.Linq;
using Numerics.Mathematics.LinearAlgebra;

namespace Numerics.Functions
{
    /// <summary>
    /// Controller for managing independent link functions applied to each element of a parameter vector.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Each parameter index can have its own link function, or null for no transformation (identity).
    ///     This enables variance stabilization, constraint enforcement, and improved MCMC mixing
    ///     by applying per-parameter transformations.
    /// </para>
    /// <para>
    ///     The Jacobian and log-determinant methods support the change-of-variables formula
    ///     for transformed-space probability calculations: p(&#966;) = p(&#952;) |&#8706;&#952;/&#8706;&#966;|.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class LinkController
    {
        /// <summary>
        /// The array of link functions, one per parameter index. Null entries indicate identity (no transformation).
        /// </summary>
        public ILinkFunction?[] Links { get; private set; }

        /// <summary>
        /// Initializes a new instance with no link functions (all parameters untransformed).
        /// </summary>
        public LinkController()
        {
            Links = Array.Empty<ILinkFunction?>();
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="LinkController"/> class with the specified link functions.
        /// </summary>
        /// <param name="links">
        /// An array of link functions, one per parameter index. Null entries indicate identity (no transformation).
        /// If no arguments are provided, all parameters pass through untransformed.
        /// </param>
        public LinkController(params ILinkFunction?[] links)
        {
            Links = links ?? Array.Empty<ILinkFunction?>();
        }

        /// <summary>
        /// Initializes a new instance from a serialized <see cref="XElement"/>.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize from.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <remarks>
        /// <para>
        ///     The XElement is expected to contain child elements named "Link" with an "Index" attribute
        ///     and optional child elements representing the link function type.
        ///     Uses <see cref="LinkFunctionFactory.CreateFromXElement(XElement)"/> for deserialization,
        ///     which supports only standard Numerics link types. For BestFit-specific types (SESLink, etc.),
        ///     use the <see cref="LinkController(ILinkFunction?[])"/> constructor with custom deserialization.
        /// </para>
        /// </remarks>
        public LinkController(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));
            var slots = xElement.Elements("Link").ToList();
            Links = new ILinkFunction?[slots.Count];
            foreach (var slot in slots)
            {
                var indexAttr = slot.Attribute("Index");
                if (indexAttr == null) continue;
                if (!int.TryParse(indexAttr.Value, out int index)) continue;
                if (index < 0 || index >= Links.Length) continue;
                var child = slot.Elements().FirstOrDefault();
                if (child != null)
                    Links[index] = LinkFunctionFactory.CreateFromXElement(child);
            }
        }

        /// <summary>
        /// Gets the number of link functions registered.
        /// </summary>
        public int Count => Links.Length;

        /// <summary>
        /// Gets the link function at the specified parameter index, or null if no link is assigned.
        /// </summary>
        /// <param name="index">The zero-based parameter index.</param>
        /// <returns>The link function at the index, or null if the index is out of range or no link is assigned.</returns>
        public ILinkFunction? this[int index] =>
            index >= 0 && index < Links.Length ? Links[index] : null;

        /// <summary>
        /// Applies the link functions element-wise: &#951;[i] = h_i(x[i]).
        /// </summary>
        /// <param name="x">The parameter vector in real-space coordinates.</param>
        /// <returns>The transformed parameter vector in link-space.</returns>
        public double[] Link(double[] x)
        {
            var eta = (double[])x.Clone();
            int n = Math.Min(x.Length, Links.Length);
            for (int i = 0; i < n; i++)
            {
                if (Links[i] != null)
                    eta[i] = Links[i]!.Link(x[i]);
            }
            return eta;
        }

        /// <summary>
        /// Applies the inverse link functions element-wise: x[i] = h_i&#8315;&#185;(&#951;[i]).
        /// </summary>
        /// <param name="eta">The parameter vector in link-space coordinates.</param>
        /// <returns>The parameter vector in real-space.</returns>
        public double[] InverseLink(double[] eta)
        {
            var x = (double[])eta.Clone();
            int n = Math.Min(eta.Length, Links.Length);
            for (int i = 0; i < n; i++)
            {
                if (Links[i] != null)
                    x[i] = Links[i]!.InverseLink(eta[i]);
            }
            return x;
        }

        /// <summary>
        /// Computes the diagonal Jacobian matrix of the link transformation.
        /// </summary>
        /// <param name="x">The parameter vector in real-space coordinates.</param>
        /// <returns>A diagonal matrix with elements d&#951;_i/d&#952;_i for each parameter.</returns>
        /// <remarks>
        /// <para>
        ///     The Jacobian is diagonal because each parameter has an independent link function.
        ///     Off-diagonal elements are zero.
        /// </para>
        /// </remarks>
        public Matrix LinkJacobian(double[] x)
        {
            int p = x.Length;
            var G = Matrix.Identity(p);
            int n = Math.Min(p, Links.Length);
            for (int i = 0; i < n; i++)
            {
                if (Links[i] != null)
                    G[i, i] = Links[i]!.DLink(x[i]);
            }
            return G;
        }

        /// <summary>
        /// Computes the log-determinant of the inverse Jacobian |&#8706;&#952;/&#8706;&#966;| for the change-of-variables formula.
        /// </summary>
        /// <param name="phi">The parameter vector in link-space (transformed coordinates).</param>
        /// <returns>log|det(&#8706;&#952;/&#8706;&#966;)|, the log absolute determinant of the inverse Jacobian.</returns>
        /// <remarks>
        /// <para>
        ///     Used in transformed-space probability calculations: p(&#966;) = p(&#952;) |&#8706;&#952;/&#8706;&#966;|.
        /// </para>
        /// <para>
        ///     For diagonal link Jacobians: log|det| = &#8722;&#8721; log|d&#951;_j/d&#952;_j|.
        /// </para>
        /// </remarks>
        public double LogDetJacobian(double[] phi)
        {
            const double tiny = 1e-16;
            double sum = 0.0;
            int n = Math.Min(phi.Length, Links.Length);
            for (int i = 0; i < n; i++)
            {
                if (Links[i] != null)
                {
                    double theta_i = Links[i]!.InverseLink(phi[i]);
                    double dEta_dTheta = Links[i]!.DLink(theta_i);
                    // Use Math.Abs to handle both increasing and decreasing link functions.
                    sum -= Math.Log(Math.Max(Math.Abs(dEta_dTheta), tiny));
                }
            }
            return sum;
        }

        /// <summary>
        /// Serializes the link controller to an <see cref="XElement"/>.
        /// </summary>
        /// <returns>An XElement containing the serialized link functions with their indices.</returns>
        public XElement ToXElement()
        {
            var element = new XElement(nameof(LinkController));
            for (int i = 0; i < Links.Length; i++)
            {
                var slot = new XElement("Link");
                slot.SetAttributeValue("Index", i);
                if (Links[i] != null)
                    slot.Add(Links[i]!.ToXElement());
                element.Add(slot);
            }
            return element;
        }

        /// <summary>
        /// Creates a <see cref="LinkController"/> for the standard 3-parameter (location, scale, shape) case.
        /// </summary>
        /// <param name="locationLink">Optional link function for the location parameter.</param>
        /// <param name="scaleLink">Optional link function for the scale parameter.</param>
        /// <param name="shapeLink">Optional link function for the shape parameter.</param>
        /// <returns>A new <see cref="LinkController"/> configured for three parameters.</returns>
        public static LinkController ForLocationScaleShape(
            ILinkFunction? locationLink = null,
            ILinkFunction? scaleLink = null,
            ILinkFunction? shapeLink = null)
        {
            return new LinkController(locationLink, scaleLink, shapeLink);
        }
    }
}
