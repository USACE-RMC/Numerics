namespace Numerics.MachineLearning
{
    /// <summary>
    /// A decision node class.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class DecisionNode
    {
        /// <summary>
        /// The feature index. 
        /// </summary>
        public int FeatureIndex { get; set; } = -1;

        /// <summary>
        /// The threshold used to split the node.
        /// </summary>
        public double Threshold { get; set; } = double.NaN;

        /// <summary>
        /// Nodes to the left of the threshold.
        /// </summary>
        public DecisionNode? Left { get; set; } = null;

        /// <summary>
        /// Nodes to the right of the threshold.
        /// </summary>
        public DecisionNode? Right { get; set; } = null;

        /// <summary>
        /// The leaf node value.
        /// </summary>
        public double Value { get; set; } = double.NaN;

        /// <summary>
        /// Determines if this is a leaf node.
        /// </summary>
        public bool IsLeafNode { get; set; } = false;
    }
}
