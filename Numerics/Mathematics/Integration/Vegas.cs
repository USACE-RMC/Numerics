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

using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Integration
{
    /// <summary>
    /// A class for adaptive Monte Carlo integration for multidimensional integration.
    /// Enhanced with Power Transform for rare event simulation.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This method aims to reduce error in Monte Carlo simulations by using a probability distribution function to concentrate the search
    /// in those areas of the integrand that make the greatest contribution. 
    /// </para>
    /// <para>
    /// <b> Power Transform Enhancement: </b>
    /// The Power Transform (γ parameter) enables efficient sampling of rare tail events without changing the integrand.
    /// For rare events (p &lt; 1e-4), set TailFocusParameter &gt; 1 to concentrate samples in tail regions.
    /// This maintains numerical stability in high dimensions (unlike z-space methods).
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// </description></item>
    /// <item><description>
    /// "Numerical Recipes: The art of Scientific Computing, Third Edition. Press et al. 2017.
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/VEGAS_algorithm"/>
    /// </description></item>
    /// </list>
    /// </remarks>
    [Serializable]
    public class Vegas : Integrator
    {
        /// <summary>
        /// Creates a new Vegas class for adaptive Monte Carlo integration for multidimensional integration.
        /// </summary>
        /// <param name="function">The multidimensional function to integrate.</param>
        /// <param name="dimensions">The number of dimensions in the function to evaluate.</param>
        /// <param name="min">The minimum values under which the integral must be computed.</param>
        /// <param name="max">The maximum values under which the integral must be computed.</param>
        public Vegas(Func<double[], double, double> function, int dimensions, IList<double> min, IList<double> max)
        {
            if (dimensions < 1) throw new ArgumentOutOfRangeException(nameof(dimensions), "There must be at least 1 dimension to evaluate.");

            // Check if the length of the min and max values equal the number of dimensions
            if (min.Count != dimensions || max.Count != dimensions)
            {
                throw new ArgumentOutOfRangeException(nameof(min), "The minimum and maximum values must be the same length as the number of dimensions.");
            }

            if (dimensions > MaxDimensions)
                throw new ArgumentOutOfRangeException(nameof(dimensions), "The maximum number of dimensions is 20.");

            // Check if the minimum values are less than the maximum values
            for (int i = 0; i < min.Count; i++)
            {
                if (max[i] <= min[i])
                {
                    throw new ArgumentOutOfRangeException(nameof(max), "The maximum values cannot be less than or equal to the minimum values.");
                }
            }

            Function = function ?? throw new ArgumentNullException(nameof(function), "The function cannot be null.");
            Dimensions = dimensions;
            Min = min.ToArray();
            Max = max.ToArray();
            Random = new MersenneTwister();
            _sobol = new SobolSequence(Dimensions);

            RelativeTolerance = 1E-3;
            InitializeParameters();

        }

        private double[] _region = null!;
        private int _numberOfBins = 50;
        private double _standardError;
        private double _chiSquared;
        private SobolSequence _sobol;

        // Constants
        private const int MaxDimensions = 20;
        private const double TinyValue = 1.0e-30;

        // Algorithm Variables

        private int i, it, j, k, gridRefinementFlag = 1, currentBinCount, previousBinCount = 1, stratificationLevels, samplesPerBinGroup;
        private double totalFunctionCalls, binVolume, binSpacing, functionValue, functionValueSquared, functionValueSquaredSum, functionValueSum, randomPointInBin, totalIntegralEstimate;
        private double totalErrorEstimate, weight, jacobianDeterminant, randomBinIndex, normalizedBinCount, binSpacingDelta, sumChiSquared, sumWeightedResults, sumIntegrationWeights;


        private int[] binIndexes = null!;
        private int[] regionBinCounts = null!;
        private double[] dimensionContributions = null!;
        private double[] dx = null!;
        private double[] binDensities = null!;
        private double[] randomPoint = null!;
        private double[] binEdges = null!;
        private double[,] binContributions = null!;
        private double[,] refinedBinContributions = null!;
        private double[,] stratificationGrid = null!;


        /// <summary>
        /// The multidimensional function to integrate.
        /// </summary>
        public Func<double[], double, double> Function { get; }

        /// <summary>
        /// The number of dimensions in the function to evaluate./>.
        /// </summary>
        public int Dimensions { get; }

        /// <summary>
        /// The minimum values under which the integral must be computed.
        /// </summary>
        public double[] Min { get; }

        /// <summary>
        /// The maximum values under which the integral must be computed. 
        /// </summary>
        public double[] Max { get; }

        /// <summary>
        /// Gets and sets the random number generator to be used within the Monte Carlo integration.
        /// </summary>
        public Random Random { get; set; }

        /// <summary>
        /// Determines whether to use a Sobol sequence or a pseudo-Random number generator. 
        /// </summary>
        public bool UseSobolSequence { get; set; } = true;

        /// <summary>
        /// Determines whether to check convergence and exit when integrating.
        /// </summary>
        public bool CheckConvergence { get; set; } = true;

        /// <summary>
        /// Determines how to initialize the Vegas routine. 
        /// </summary>
        /// <remarks>      
        /// If 0, then Vegas enters on a cold start. If Initialize 1, then inherit the grid from a previous call, but not its answers. If 2, then inherit the previous grid and its answers.
        /// </remarks>
        public int Initialize { get; set; } = 0;

        /// <summary>
        /// Gets and sets the number of statistically independent evaluations of the integral, per iteration. 
        /// </summary>
        public int IndependentEvaluations { get; set; } = 1000;

        /// <summary>
        /// Gets and sets the number of function evaluations within each independent evaluation. Default = 10,000.
        /// </summary>
        public int FunctionCalls { get; set; } = 10000;

        /// <summary>
        /// The refinement damping parameter for the grid. The default = 1.5.
        /// </summary>
        public double Alpha { get; set; } = 1.5;

        /// <summary>
        /// Gets and sets the number of stratification bins for each dimension. The default = 50.
        /// </summary>
        public int NumberOfBins
        {
            get { return _numberOfBins; }
            set
            {
                _numberOfBins = value;
                InitializeParameters();
            }
        }

        /// <summary>
        /// Power transform parameter for tail-focused rare event sampling.
        /// Default = 1.0 (standard uniform sampling, backward compatible).
        /// Set γ > 1 to focus sampling on upper tail (p → 1) for rare events.
        /// Recommended values: γ=2 (moderate focus), γ=4 (strong focus), γ=10 (very strong focus for p &lt; 1e-6).
        /// </summary>
        /// <remarks>
        /// The power transform uses p' = 1 - (1-p)^γ to concentrate samples in the upper tail.
        /// Unlike z-space transforms, this maintains numerical stability in high dimensions.
        /// Weights are corrected by Jacobian: dp'/dp = γ(1-p)^(γ-1), which stays O(1).
        /// </remarks>
        public double TailFocusParameter { get; set; } = 1.0;

        /// <summary>
        /// Gets the stratification grid boundaries. 
        /// </summary>
        public double[,] Grid => stratificationGrid;

        /// <summary>
        /// Gets integration standard error. 
        /// </summary>
        public double StandardError => _standardError;

        /// <summary>
        /// Gets the Chi-Squared statistic
        /// </summary>
        public double ChiSquared => _chiSquared;


        /// <summary>
        /// Initialize the parameter arrays.
        /// </summary>
        private void InitializeParameters()
        {
            binIndexes = new int[Dimensions];
            regionBinCounts = new int[Dimensions];
            dimensionContributions = new double[Dimensions];
            dx = new double[Dimensions];
            binDensities = new double[NumberOfBins];
            randomPoint = new double[Dimensions];
            binEdges = new double[NumberOfBins];
            binContributions = new double[NumberOfBins, Dimensions];
            refinedBinContributions = new double[NumberOfBins, Dimensions];
            stratificationGrid = new double[Dimensions, NumberOfBins];

            // Create region array
            _region = new double[2 * Dimensions];
            for (int i = 0; i < Dimensions; i++)
                _region[i] = Min[i];
            for (int i = Dimensions; i < 2 * Dimensions; i++)
                _region[i] = Max[i - Dimensions];

            Iterations = 0;
        }

        /// <summary>
        /// Apply power transform to probability for tail-focused sampling.
        /// p' = 1 - (1-p)^γ concentrates samples in upper tail when γ > 1.
        /// When γ = 1, this is identity transform (backward compatible).
        /// </summary>
        private double ApplyPowerTransform(double p)
        {
            double gamma = TailFocusParameter;

            // Identity transform when γ = 1 (standard Vegas behavior)
            if (Math.Abs(gamma - 1.0) < 1e-10)
            {
                return p;
            }

            // Power transform for upper tail focus
            // Maps [0,1] → [0,1] but concentrates samples near 1
            return 1.0 - Math.Pow(1.0 - p, gamma);
        }

        /// <summary>
        /// Compute Jacobian for power transform: dp'/dp = γ(1-p)^(γ-1).
        /// This weight correction ensures unbiased integration.
        /// When γ = 1, returns 1.0 (identity Jacobian).
        /// </summary>
        private double PowerTransformJacobian(double p)
        {
            double gamma = TailFocusParameter;

            // Identity Jacobian when γ = 1
            if (Math.Abs(gamma - 1.0) < 1e-10)
            {
                return 1.0;
            }

            // Jacobian: dp'/dp = γ(1-p)^(γ-1)
            // Clamp (1-p) to avoid numerical issues
            double oneMinusP = Math.Max(1.0 - p, 1e-15);
            return gamma * Math.Pow(oneMinusP, gamma - 1.0);
        }

        /// <summary>
        /// Configure Vegas for rare tail events.
        /// Automatically sets TailFocusParameter based on target probability.
        /// </summary>
        /// <param name="targetProbability">Target rare event probability (e.g., 1e-6)</param>
        public void ConfigureForRareEvents(double targetProbability)
        {
            // Choose γ so that target probability appears in ~5% of transformed samples
            // Solving: (1 - 0.95)^γ ≈ targetProbability
            // γ ≈ ln(targetProbability) / ln(0.05)

            double gamma = Math.Log(targetProbability) / Math.Log(0.05);
            gamma = Math.Max(1.0, Math.Min(gamma, 20.0));  // Clamp to [1, 20]

            this.TailFocusParameter = gamma;
            this.NumberOfBins = Math.Max(100, this.NumberOfBins);
            this.Alpha = 1.8;  // More aggressive grid adaptation
        }

        /// <inheritdoc/>
        public override void Integrate()
        {
            Validate();

            try
            {
                // Compute
                vegas();

                // Update status
                if (FunctionEvaluations >= MaxFunctionEvaluations)
                {
                    Status = IntegrationStatus.MaximumFunctionEvaluationsReached;
                }
                else
                {
                    Status = IntegrationStatus.Success;
                }

            }
            catch (Exception)
            {
                Status = IntegrationStatus.Failure;
                if (ReportFailure) throw;
            }

        }

        /// <summary>
        /// Helper function for Integrate(). This performs the Vegas algorithm/.
        /// </summary>
        private void vegas()
        {
            if (Initialize <= 0)
            {
                // Normal entry. Enter here on a cold start. 
                // Change flag = 0 to disable stratified sampling, i.e. use importance sampling only. 
                gridRefinementFlag = previousBinCount = 1;
                for (j = 0; j < Dimensions; j++)
                {
                    stratificationGrid[j, 0] = 1.0;
                }
            }
            if (Initialize <= 1)
            {
                // Enter here to inherit the grid from a previous call, but not its answers.
                sumWeightedResults = sumIntegrationWeights = sumChiSquared = 0.0;
            }
            if (Initialize <= 2)
            {
                // Inherit the grid and previous results for further refinement.
                currentBinCount = NumberOfBins;
                stratificationLevels = 1;

                // Set up stratification grid based on the number of function calls and dimensions
                if (gridRefinementFlag != 0)
                {
                    stratificationLevels = (int)Math.Pow(FunctionCalls / 2.0 + 0.25, 1.0 / Dimensions);
                    gridRefinementFlag = 1;

                    if ((2 * stratificationLevels - NumberOfBins) >= 0)
                    {
                        gridRefinementFlag = -1;
                        samplesPerBinGroup = stratificationLevels / NumberOfBins + 1;
                        currentBinCount = stratificationLevels / samplesPerBinGroup;
                        stratificationLevels = samplesPerBinGroup * currentBinCount;
                    }
                }

                int totalSampleGroups = 1;
                for (i = 0; i < Dimensions; i++)
                {
                    totalSampleGroups *= stratificationLevels;
                }

                samplesPerBinGroup = Math.Max((int)(FunctionCalls / totalSampleGroups), 2);
                totalFunctionCalls = (double)(samplesPerBinGroup) * (double)(totalSampleGroups);

                binVolume = 1;
                binSpacing = 1.0 / stratificationLevels;
                for (i = 0; i < Dimensions; i++)
                {
                    binVolume *= binSpacing;
                }

                binVolume = Tools.Sqr(totalFunctionCalls * binVolume) / samplesPerBinGroup / samplesPerBinGroup / (samplesPerBinGroup - 1.0);
                normalizedBinCount = currentBinCount;
                binSpacing *= normalizedBinCount;

                jacobianDeterminant = 1.0 / totalFunctionCalls;
                for (j = 0; j < Dimensions; j++)
                {
                    dx[j] = Max[j] - Min[j];
                    jacobianDeterminant *= dx[j];
                }

                // Refine the grid binning if necessary
                if (currentBinCount != previousBinCount)
                {
                    for (i = 0; i < Math.Max(currentBinCount, previousBinCount); i++)
                    {
                        binDensities[i] = 1.0;
                    }
                    for (j = 0; j < Dimensions; j++)
                    {
                        RefineGrid(previousBinCount / normalizedBinCount, currentBinCount, binDensities, binEdges, stratificationGrid, j);
                    }
                    previousBinCount = currentBinCount;
                }

            }
            // Main iteration loop. Can enter here (init >= 3) to do additional independent evaluations with all other parameters unchanged.
            for (it = 0; it < IndependentEvaluations; it++)
            {
                Iterations++;

                totalIntegralEstimate = totalErrorEstimate = 0.0;

                // Initialize bin contributions and stratification indexes 
                for (j = 0; j < Dimensions; j++)
                {
                    regionBinCounts[j] = 1;
                    for (i = 0; i < currentBinCount; i++)
                    {
                        binContributions[i, j] = refinedBinContributions[i, j] = 0.0;
                    }
                }

                // Stratified sampling loop
                while (true)
                {
                    functionValueSum = functionValueSquaredSum = 0.0;
                    for (k = 0; k < samplesPerBinGroup; k++)
                    {
                        weight = jacobianDeterminant;

                        // Generate random or Sobol sequences
                        double[] rnd = UseSobolSequence ? _sobol.NextDouble() : Random.NextDoubles(Dimensions);

                        for (j = 0; j < Dimensions; j++)
                        {
                            randomBinIndex = (regionBinCounts[j] - rnd[j]) * binSpacing + 1.0;
                            binIndexes[j] = Math.Max(Math.Min((int)randomBinIndex, NumberOfBins), 1);

                            if (binIndexes[j] > 1)
                            {
                                binSpacingDelta = stratificationGrid[j, binIndexes[j] - 1] - stratificationGrid[j, binIndexes[j] - 2];
                                randomPointInBin = stratificationGrid[j, binIndexes[j] - 2] + (randomBinIndex - binIndexes[j]) * binSpacingDelta;
                            }
                            else
                            {
                                binSpacingDelta = stratificationGrid[j, binIndexes[j] - 1];
                                randomPointInBin = (randomBinIndex - binIndexes[j]) * binSpacingDelta;
                            }

                            // MODIFIED: Apply power transform if enabled (γ > 1)
                            // Grid position is in [0, 1]
                            double p_uniform = randomPointInBin;

                            // Transform to focus on tails (identity if γ = 1)
                            double p_transformed = ApplyPowerTransform(p_uniform);

                            // Scale to [Min, Max] range
                            randomPoint[j] = Min[j] + p_transformed * dx[j];

                            // Weight includes Jacobian correction (1.0 if γ = 1)
                            double jacobian = PowerTransformJacobian(p_uniform);
                            weight *= binSpacingDelta * normalizedBinCount * jacobian;
                        }

                        // Evaluate integrand function
                        functionValue = weight * Function(randomPoint, weight);
                        FunctionEvaluations++;

                        functionValueSum += functionValue;
                        functionValueSquared = functionValue * functionValue;
                        functionValueSquaredSum += functionValueSquared;

                        for (j = 0; j < Dimensions; j++)
                        {
                            refinedBinContributions[binIndexes[j] - 1, j] += functionValue;
                            if (gridRefinementFlag >= 0)
                            {
                                binContributions[binIndexes[j] - 1, j] += functionValueSquared;
                            }
                        }
                    }

                    functionValueSquaredSum = Math.Sqrt(functionValueSquaredSum * samplesPerBinGroup);
                    functionValueSquaredSum = (functionValueSquaredSum - functionValueSum) * (functionValueSquaredSum + functionValueSum);
                    if (functionValueSquaredSum <= 0.0) functionValueSquaredSum = TinyValue;

                    totalIntegralEstimate += functionValueSum;
                    totalErrorEstimate += functionValueSquaredSum;

                    if (gridRefinementFlag < 0)
                    {
                        for (j = 0; j < Dimensions; j++)
                        {
                            binContributions[binIndexes[j] - 1, j] += functionValueSquaredSum;
                        }
                    }
                    // Increment stratification indexes
                    for (k = Dimensions - 1; k >= 0; k--)
                    {
                        regionBinCounts[k] %= stratificationLevels;
                        if (++regionBinCounts[k] != 1) break;
                    }
                    if (k < 0) break;
                }

                // Compute final results for this iteration
                totalErrorEstimate *= binVolume;
                weight = 1.0 / totalErrorEstimate;
                sumIntegrationWeights += weight;
                sumWeightedResults += weight * totalIntegralEstimate;
                sumChiSquared += weight * totalIntegralEstimate * totalIntegralEstimate;

                _result = sumWeightedResults / sumIntegrationWeights;
                _standardError = Math.Sqrt(1.0 / sumIntegrationWeights);
                totalErrorEstimate = Math.Sqrt(totalErrorEstimate);

                // Ensure chi-squared is non-negative
                _chiSquared = (sumChiSquared - sumWeightedResults * _result) / (it + 1);
                if (_chiSquared < 0.0) _chiSquared = 0.0;

                // check convergence
                if ((CheckConvergence && Iterations > 1 && Math.Abs(_standardError / _result) < RelativeTolerance) || FunctionEvaluations >= MaxFunctionEvaluations)
                {
                    break;
                }

                // Refine the grid. Consult references to understand the subtlety of this procedure. The refinement is damped.
                // to avoid rapid, destabilizing changes, and also compressed in range by the exponent ALPH.
                for (j = 0; j < Dimensions; j++)
                {
                    double lowerContribution = binContributions[0, j];
                    double upperContribution = binContributions[1, j];
                    binContributions[0, j] = (lowerContribution + upperContribution) / 2.0;
                    dimensionContributions[j] = binContributions[0, j];
                    for (i = 2; i < currentBinCount; i++)
                    {
                        double regionContribution = lowerContribution + upperContribution;
                        lowerContribution = upperContribution;
                        upperContribution = binContributions[i, j];
                        binContributions[i - 1, j] = (regionContribution + upperContribution) / 3.0;
                        dimensionContributions[j] += binContributions[i - 1, j];
                    }

                    binContributions[currentBinCount - 1, j] = (lowerContribution + upperContribution) / 2.0;
                    dimensionContributions[j] += binContributions[currentBinCount - 1, j];
                }
                for (j = 0; j < Dimensions; j++)
                {
                    double regionContribution = 0.0;
                    for (i = 0; i < currentBinCount; i++)
                    {
                        if (binContributions[i, j] < TinyValue) binContributions[i, j] = TinyValue;
                        binDensities[i] = Math.Pow((1.0 - binContributions[i, j] / dimensionContributions[j]) / (Math.Log(dimensionContributions[j]) - Math.Log(binContributions[i, j])), Alpha);
                        regionContribution += binDensities[i];
                    }
                    RefineGrid(regionContribution / normalizedBinCount, currentBinCount, binDensities, binEdges, stratificationGrid, j);
                }

            }

        }

        /// <summary>
        /// Refines the grid by redistributing bin boundaries based on the contribution densities.
        /// </summary>
        /// <param name="totalDensity">Total density to distribute across bins.</param>
        /// <param name="numberOfBins">Number of bins in the current dimension.</param>
        /// <param name="binDensities">Array containing the density values for each bin.</param>
        /// <param name="newBinEdges">Array to store the adjusted bin edges.</param>
        /// <param name="stratificationGrid">2D grid containing the stratification for all dimensions.</param>
        /// <param name="dimensionIndex">Index of the dimension being refined.</param>
        private void RefineGrid(double totalDensity, int numberOfBins, double[] binDensities, double[] newBinEdges, double[,] stratificationGrid, int dimensionIndex)
        {
            int currentBin = 0;
            double cumulativeDensity = 0.0, upperEdge = 0.0, lowerEdge = 0.0;
            for (int i = 0; i < numberOfBins - 1; i++)
            {
                // Accumulate densities until it exceeds the target for the current boundary
                while (totalDensity > cumulativeDensity)
                {
                    cumulativeDensity += binDensities[(++currentBin) - 1];
                }

                // Compute the new bin edge
                if (currentBin > 1)
                {
                    lowerEdge = stratificationGrid[dimensionIndex, currentBin - 2];
                }
                upperEdge = stratificationGrid[dimensionIndex, currentBin - 1];
                cumulativeDensity -= totalDensity;

                // Calculate the adjusted bin edge based on cumulative density
                newBinEdges[i] = upperEdge - (upperEdge - lowerEdge) * cumulativeDensity / binDensities[currentBin - 1];
            }
            // Update the stratification grid
            for (int i = 0; i < numberOfBins - 1; i++)
            {
                stratificationGrid[dimensionIndex, i] = newBinEdges[i];
            }
            stratificationGrid[dimensionIndex, numberOfBins - 1] = 1.0;
        }
    }
}