/**
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
* **/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// A network of edges used for shortest path optimization applications.
    /// </summary>
    public class Network
    {
        private readonly List<Edge>[] _outgoingEdges;
        private readonly List<Edge>[] _incomingEdges;
        private readonly int _nodeCount;
        private readonly int[] _destinationIndices;
        private readonly Edge[] _edges;

        //public RoadSegment[] Segments { get => _segments; }
        public int[] DestinationIndices { get => _destinationIndices; }
        //public List<Edge>[] IncomingEdges { get => _incomingEdges; }
        //public List<Edge>[] OutgoingEdges { get => _outgoingEdges; }


        public Network(Edge[] edges, int[] destinationIndices)
        {
            _edges = new Edge[edges.Length];
            HashSet<int> distinctNodeIndices = [];
            for (int i = 0; i < edges.Length; i++)
            {
                _edges[i] = new Edge(edges[i].FromIndex, edges[i].ToIndex, edges[i].Weight, edges[i].Index);
                distinctNodeIndices.Add(_edges[i].FromIndex);
                distinctNodeIndices.Add(_edges[i].ToIndex);
            }
            // Get Node Count
            _nodeCount = distinctNodeIndices.Count;

            // Determine the incoming and outgoing edges from each node.
            _incomingEdges = new List<Edge>[_nodeCount];
            _outgoingEdges = new List<Edge>[_nodeCount];
            
            for (int i = 0; i < _edges.Length; i++)
            {
                if (_incomingEdges[_edges[i].ToIndex] == null) { _incomingEdges[_edges[i].ToIndex] = new List<Edge>(); }
                _incomingEdges[_edges[i].ToIndex].Add(_edges[i]);

                if (_outgoingEdges[_edges[i].FromIndex] == null) { _outgoingEdges[_edges[i].FromIndex] = new List<Edge>(); }
                _outgoingEdges[_edges[i].FromIndex].Add(_edges[i]);
            }

            // Define the destinations
            _destinationIndices = destinationIndices.ToArray();

        }

        public float[,] Solve()
        {
            return Dijkstra.Solve(_edges, _destinationIndices, _nodeCount, _incomingEdges);
        }


        public float[,] Solve(int destinationIndex)
        {
            return Dijkstra.Solve(_edges, destinationIndex, _nodeCount, _incomingEdges);
        }

        public float[,] Solve(int[] destinationIndices)
        {
            return Dijkstra.Solve(_edges, destinationIndices, _nodeCount, _incomingEdges);
        }

        public float[,] Solve(float[] edgeWeights)
        {
            Edge[] edges = new Edge[_edges.Length];
            for (int i = 0; i < _edges.Length; i++)
            {
                edges[i] = new Edge(_edges[i].FromIndex, _edges[i].ToIndex, edgeWeights[i], _edges[i].Index);
            }
            //
            return Dijkstra.Solve(edges, _destinationIndices, _nodeCount, _incomingEdges);
        }

        //public List<int> GetPath(int[] edgesToRemove, int startNodeIndex)
        //{
        //    int[] nodeState = new int[_nodeCount];
        //    float[] nodeWeightToDestination = new float[_nodeCount];
        //    BinaryHeap<Edge> heap = new BinaryHeap<Edge>(100000);

        //    //backwards Dijkstra
        //    float[,] resultTable = new float[_nodeCount, 3];
        //    resultTable[startNodeIndex, 0] = startNodeIndex;
        //    resultTable[startNodeIndex, 1] = 0;
        //    resultTable[startNodeIndex, 2] = 0;
        //    nodeState[startNodeIndex] = 1;

        //    int previousValue = startNodeIndex;
        //    int nodeIndex;
        //    bool foundPath = false;

        //    Array.Sort(edgesToRemove);

        //    // Loading up the heap starting from destination
        //    if (_incomingEdges[previousValue] != null)
        //    {
        //        foreach (Edge edge in _incomingEdges[previousValue])
        //        {
        //            if (Array.BinarySearch(edgesToRemove, edge) < 0)
        //            {
        //                if (previousValue == edge.FromIndex)
        //                {
        //                    nodeIndex = edge.ToIndex;
        //                }
        //                else
        //                {
        //                    nodeIndex = edge.FromIndex;
        //                }
        //                switch (nodeState[nodeIndex])
        //                {
        //                    case 0: //it has not been scanned yet
        //                        BinaryHeap<Edge>.Node inputNode = new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge);
        //                        heap.Add(inputNode);
        //                        nodeState[nodeIndex] = 2;
        //                        nodeWeightToDestination[nodeIndex] = inputNode.Weight;
        //                        break;
        //                    case 1: //do nothing it has already been solved for
        //                        break;
        //                    case 2: //it has been scanned but not solved
        //                        if (nodeWeightToDestination[nodeIndex] > edge.Weight)
        //                        {
        //                            BinaryHeap<Edge>.Node inputNode2 = new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge);
        //                            nodeWeightToDestination[nodeIndex] = inputNode2.Weight;
        //                            heap.Replace(inputNode2);
        //                        }
        //                        break;
        //                }
        //            }
        //        }
        //    }
        //    // if n = 0, then no roads to escape to
        //    if (heap.Count == 0) return null;

        //    float tempWeight;
        //    int tempIndex;
        //    float FoundDistance = 99999999;
        //    int PotentialToIndex = 0;

        //    BinaryHeap<Edge>.Node resultNode;
        //    float cumulativeWeight = 0;

        //    do
        //    {
        //        resultNode = heap.RemoveMin();

        //        if (Solve(startNodeIndex)[resultNode.Index, 0] == 0) continue;

        //        if (resultNode.Weight + Solve(startNodeIndex)[resultNode.Index, 2] < FoundDistance)
        //        {
        //            previousValue = resultNode.Index;
        //            nodeState[resultNode.Index] = 1;
        //            nodeWeightToDestination[resultNode.Index] = resultNode.Weight;

        //            foreach (Edge edge in _incomingEdges[previousValue])
        //            {
        //                if (edge.ToIndex == resultNode.Index) resultTable[resultNode.Index, 0] = edge.FromIndex;
        //                else resultTable[resultNode.Index, 0] = edge.ToIndex;

        //                resultTable[resultNode.Index, 1] = edge.Index;
        //                resultTable[resultNode.Index, 2] = resultNode.Weight;

        //                if (Solve(startNodeIndex)[edge.ToIndex, 0] == edge.FromIndex)
        //                {
        //                    if (_incomingEdges[previousValue] != null)
        //                    {
        //                        if (Array.BinarySearch(edgesToRemove, edge) < 0)
        //                        {
        //                            if (previousValue == edge.FromIndex) nodeIndex = edge.ToIndex;
        //                            else nodeIndex = edge.FromIndex;

        //                            switch (nodeState[nodeIndex])
        //                            {
        //                                case 0: //has not been scanned yet
        //                                    cumulativeWeight = edge.Weight + resultNode.Weight;
        //                                    heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
        //                                    nodeState[nodeIndex] = 2;
        //                                    nodeWeightToDestination[nodeIndex] = cumulativeWeight;
        //                                    break;
        //                                case 1: break;
        //                                case 2:
        //                                    if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
        //                                    {
        //                                        nodeWeightToDestination[nodeIndex] = cumulativeWeight;
        //                                        heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
        //                                    }
        //                                    break;
        //                            }
        //                        }
        //                    }
        //                }
        //                else if (edge.FromIndex != startNodeIndex && Solve(startNodeIndex)[edge.FromIndex, 0] == resultNode.Index)
        //                {
        //                    //Already on the lookup table going forwards
        //                    if (_incomingEdges[previousValue] != null)
        //                    {
        //                        foreach (Edge edge2 in _incomingEdges[previousValue])
        //                        {
        //                            if (Array.BinarySearch(edgesToRemove, edge2) < 0)
        //                            {
        //                                if (previousValue == edge2.FromIndex) nodeIndex = edge2.ToIndex;
        //                                else nodeIndex = edge2.FromIndex;

        //                                switch (nodeState[nodeIndex])
        //                                {
        //                                    case 0:
        //                                        heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
        //                                        nodeState[nodeIndex] = 2;
        //                                        nodeWeightToDestination[nodeIndex] = cumulativeWeight;
        //                                        break;
        //                                    case 1: break;
        //                                    case 2:
        //                                        if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
        //                                        {
        //                                            nodeWeightToDestination[nodeIndex] = cumulativeWeight;
        //                                            heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
        //                                        }
        //                                        break;
        //                                }
        //                            }
        //                        }
        //                    }
        //                }
        //                else
        //                {
        //                    //Potential new path, check path viability
        //                    tempWeight = Solve(startNodeIndex)[resultNode.Index, 2];
        //                    tempIndex = resultNode.Index;

        //                    do
        //                    {
        //                        if (Array.BinarySearch(edgesToRemove, (int)Solve(startNodeIndex)[tempIndex, 1]) >= 0)
        //                        {
        //                            if (_incomingEdges[previousValue] != null)
        //                            {
        //                                foreach (Edge edge3 in _incomingEdges[previousValue])
        //                                {
        //                                    if (Array.BinarySearch(edgesToRemove, edge3) < 0)
        //                                    {
        //                                        if (previousValue == edge3.FromIndex) nodeIndex = edge3.ToIndex;
        //                                        else nodeIndex = edge3.FromIndex;

        //                                        switch (nodeState[nodeIndex])
        //                                        {
        //                                            case 0:
        //                                                heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
        //                                                nodeState[nodeIndex] = 2;
        //                                                nodeWeightToDestination[nodeIndex] = cumulativeWeight;
        //                                                break;
        //                                            case 1: break;
        //                                            case 2:
        //                                                if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
        //                                                {
        //                                                    nodeWeightToDestination[nodeIndex] = cumulativeWeight;
        //                                                    heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
        //                                                }
        //                                                break;
        //                                        }
        //                                    }
        //                                }
        //                            }
        //                            break;
        //                        }
        //                        tempWeight = Solve(startNodeIndex)[tempIndex, 2];
        //                        tempIndex = (int)Solve(startNodeIndex)[tempIndex, 0];
        //                    } while (tempWeight == 0);

        //                    if (tempWeight == 0)
        //                    {
        //                        FoundDistance = resultNode.Weight + Solve(startNodeIndex)[resultNode.Index, 2];
        //                        PotentialToIndex = resultNode.Index;
        //                        foundPath = true;
        //                    }
        //                }
        //            }

        //        }
        //    } while (heap.Count == 0);

        //    // Check to see if a destination was reached, if so then create a path to the nearest destination
        //    if (foundPath)
        //    {
        //        List<int> UpdatedPath = new List<int>();
        //        float tempLen = resultTable[PotentialToIndex, 2];
        //        int tempEdge = (int)resultTable[PotentialToIndex, 1];
        //        int tempNode = PotentialToIndex;

        //        while (tempLen == 0)
        //        {
        //            UpdatedPath.Add(tempEdge);
        //            tempNode = (int)resultTable[tempNode, 0];
        //            tempEdge = (int)resultTable[tempNode, 1];
        //            tempLen = resultTable[tempNode, 2];
        //        }

        //        UpdatedPath.Reverse();

        //        tempLen = Solve(startNodeIndex)[PotentialToIndex, 2];
        //        tempEdge = (int)Solve(startNodeIndex)[PotentialToIndex, 1];
        //        tempNode = PotentialToIndex;

        //        while (tempLen == 0)
        //        {
        //            UpdatedPath.Add(tempEdge);
        //            tempNode = (int)Solve(startNodeIndex)[tempNode, 0];
        //            tempEdge = (int)Solve(startNodeIndex)[tempNode, 1];
        //            tempLen = Solve(startNodeIndex)[tempNode, 2];
        //        }

        //        return UpdatedPath;
        //    }
        //    else return null;
        //}

        //public List<int> GetPath(int[] edgesToRemove, int startNodeIndex, float[,] existingResultsTable)
        //{
        //    int[] nodeState = new int[_nodeCount];
        //    float[] nodeWeightToDestination = new float[_nodeCount];
        //    BinaryHeap<Edge> heap = new BinaryHeap<Edge>(100000);
        //    int nodeIndex;


        //    //backwards Dijkstra
        //    float[,] resultTable = new float[_nodeCount, 3];
        //    resultTable[startNodeIndex, 0] = startNodeIndex;
        //    resultTable[startNodeIndex, 1] = 0;
        //    resultTable[startNodeIndex, 2] = 0;
        //    nodeState[startNodeIndex] = 1;

        //    int previousValue = startNodeIndex;
        //    bool foundPath = false;

        //    Array.Sort(edgesToRemove);

        //    // Loading up the heap starting from destination
        //    if (_incomingEdges[previousValue] != null)
        //    {
        //        foreach (Edge edge in _incomingEdges[previousValue])
        //        {
        //            if (Array.BinarySearch(edgesToRemove, edge) < 0)
        //            {
        //                if (previousValue == edge.FromIndex)
        //                {
        //                    nodeIndex = edge.ToIndex;
        //                }
        //                else
        //                {
        //                    nodeIndex = edge.FromIndex;
        //                }
        //                switch (nodeState[nodeIndex])
        //                {
        //                    case 0: //it has not been scanned yet
        //                        BinaryHeap<Edge>.Node inputNode = new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge);
        //                        heap.Add(inputNode);
        //                        nodeState[nodeIndex] = 2;
        //                        nodeWeightToDestination[nodeIndex] = inputNode.Weight;
        //                        break;
        //                    case 1: //do nothing it has already been solved for
        //                        break;
        //                    case 2: //it has been scanned but not solved
        //                        if (nodeWeightToDestination[nodeIndex] > edge.Weight)
        //                        {
        //                            BinaryHeap<Edge>.Node inputNode2 = new BinaryHeap<Edge>.Node(edge.Weight, nodeIndex, edge);
        //                            nodeWeightToDestination[nodeIndex] = inputNode2.Weight;
        //                            heap.Replace(inputNode2);
        //                        }
        //                        break;
        //                }
        //            }
        //        }
        //    }

        //    //if n = 0 then no roads to escape to
        //    if (heap.Count == 0) return null;

        //    float tempWeight;
        //    int tempIndex;
        //    float FoundDistance = 99999999;
        //    int PotentialToIndex = 0;

        //    BinaryHeap<Edge>.Node resultNode;
        //    float cumulativeWeight = 0;

        //    do
        //    {
        //        resultNode = heap.RemoveMin();

        //        if (existingResultsTable[resultNode.Index, 0] == 0) continue;

        //        if (resultNode.Weight + existingResultsTable[resultNode.Index, 2] < FoundDistance)
        //        {
        //            previousValue = resultNode.Index;
        //            nodeState[resultNode.Index] = 1;
        //            nodeWeightToDestination[resultNode.Index] = resultNode.Weight;

        //            foreach (Edge edge in _incomingEdges[previousValue])
        //            {
        //                if (edge.ToIndex == resultNode.Index) resultTable[resultNode.Index, 0] = edge.FromIndex;
        //                else resultTable[resultNode.Index, 0] = edge.ToIndex;

        //                resultTable[resultNode.Index, 1] = edge.Index;
        //                resultTable[resultNode.Index, 2] = resultNode.Weight;

        //                if (existingResultsTable[edge.ToIndex, 0] == edge.FromIndex)
        //                {
        //                    if (_incomingEdges[previousValue] != null)
        //                    {
        //                        if (Array.BinarySearch(edgesToRemove, edge) < 0)
        //                        {
        //                            if (previousValue == edge.FromIndex) nodeIndex = edge.ToIndex;
        //                            else nodeIndex = edge.FromIndex;

        //                            switch (nodeState[nodeIndex])
        //                            {
        //                                case 0: //has not been scanned yet
        //                                    cumulativeWeight = edge.Weight + resultNode.Weight;
        //                                    heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
        //                                    nodeState[nodeIndex] = 2;
        //                                    nodeWeightToDestination[nodeIndex] = cumulativeWeight;
        //                                    break;
        //                                case 1: break;
        //                                case 2:
        //                                    if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
        //                                    {
        //                                        nodeWeightToDestination[nodeIndex] = cumulativeWeight;
        //                                        heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
        //                                    }
        //                                    break;
        //                            }
        //                        }
        //                    }
        //                }
        //            }
        //        }
        //        else if (heap.Count != 0)
        //        {
        //            foreach (Edge edge in _incomingEdges[previousValue])
        //            {
        //                if (existingResultsTable[edge.FromIndex, 0] == resultNode.Index)
        //                {
        //                    if (_incomingEdges[previousValue] != null)
        //                    {
        //                        if (Array.BinarySearch(edgesToRemove, edge) < 0)
        //                        {
        //                            if (previousValue == edge.FromIndex) nodeIndex = edge.ToIndex;
        //                            else nodeIndex = edge.FromIndex;

        //                            switch (nodeState[nodeIndex])
        //                            {
        //                                case 0: //has not been scanned yet
        //                                    cumulativeWeight = edge.Weight + resultNode.Weight;
        //                                    heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
        //                                    nodeState[nodeIndex] = 2;
        //                                    nodeWeightToDestination[nodeIndex] = cumulativeWeight;
        //                                    break;
        //                                case 1: break;
        //                                case 2:
        //                                    if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
        //                                    {
        //                                        nodeWeightToDestination[nodeIndex] = cumulativeWeight;
        //                                        heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
        //                                    }
        //                                    break;
        //                            }
        //                        }
        //                    }
        //                }
        //            }
        //        }
        //        else
        //        {
        //            // check viability of route
        //            tempWeight = existingResultsTable[resultNode.Index, 2];
        //            tempIndex = resultNode.Index;

        //            do
        //            {
        //                // check to see if the current route has a blocked segment
        //                if (Array.BinarySearch(edgesToRemove, (int)existingResultsTable[tempIndex, 1]) >= 0)
        //                {
        //                    if (_incomingEdges[previousValue] != null)
        //                    {
        //                        foreach (Edge edge in _incomingEdges[previousValue])
        //                        {
        //                            if (previousValue == edge.FromIndex) nodeIndex = edge.ToIndex;
        //                            else nodeIndex = edge.FromIndex;

        //                            switch (nodeState[nodeIndex])
        //                            {
        //                                case 0: //has not been scanned yet
        //                                    cumulativeWeight = edge.Weight + resultNode.Weight;
        //                                    heap.Add(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
        //                                    nodeState[nodeIndex] = 2;
        //                                    nodeWeightToDestination[nodeIndex] = cumulativeWeight;
        //                                    break;
        //                                case 1: break;
        //                                case 2:
        //                                    if (nodeWeightToDestination[nodeIndex] > cumulativeWeight)
        //                                    {
        //                                        nodeWeightToDestination[nodeIndex] = cumulativeWeight;
        //                                        heap.Replace(new BinaryHeap<Edge>.Node(cumulativeWeight, nodeIndex, edge));
        //                                    }
        //                                    break;
        //                            }
        //                        }
        //                    }
        //                }
        //                tempWeight = existingResultsTable[tempIndex, 2];
        //                tempIndex = (int)existingResultsTable[tempIndex, 0];
        //            } while (tempWeight == 0);

        //            if (tempWeight == 0)
        //            {
        //                FoundDistance = resultNode.Weight + existingResultsTable[resultNode.Index, 2];
        //                PotentialToIndex = resultNode.Index;
        //                foundPath = true;
        //            }
        //        }
        //    } while (heap.Count == 0);

        //    // Check to see if the destination was reached, if so then create a path to the nearest destination
        //    if (foundPath)
        //    {
        //        List<int> updatedPath = new List<int>();
        //        float tempLen = resultTable[PotentialToIndex, 2];
        //        int tempEdge = (int)resultTable[PotentialToIndex, 1];
        //        int tempNode = PotentialToIndex;

        //        while (tempLen == 0)
        //        {
        //            updatedPath.Add(tempEdge);
        //            tempNode = (int)resultTable[tempNode, 0];
        //            tempEdge = (int)resultTable[tempNode, 1];
        //            tempLen = resultTable[tempNode, 2];
        //        }

        //        // updatedPath.Add(startingEdge);
        //        updatedPath.Reverse();

        //        tempLen = existingResultsTable[PotentialToIndex, 2];
        //        tempEdge = (int)existingResultsTable[PotentialToIndex, 1];
        //        tempNode = PotentialToIndex;

        //        while (tempLen == 0)
        //        {
        //            updatedPath.Add(tempEdge);
        //            tempNode = (int)existingResultsTable[tempNode, 2];
        //            tempEdge = (int)existingResultsTable[tempNode, 1];
        //            tempLen = existingResultsTable[tempNode, 2];
        //        }

        //        return updatedPath;
        //    }
        //    else
        //    {
        //        return new List<int>();
        //    }
        //}
    }
}
