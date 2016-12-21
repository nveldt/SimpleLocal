# SimpleLocal
Implementation of SimpleLocal for Local Cut Improvement in Graphs

This repository includes code and data from the paper:

"A Simple and Strongly-Local Flow-Based Method for Cut Improvement"
 
Nate Veldt, David Gleich, Michael Mahoney 
Proceedings of The 33rd International Conference on Machine Learning, pp. 1938–1947, 2016

Included is an implementation of the SimpleLocal algorithm, and the adjacency matrix subgraph of the brain network used in numerical experiments. The full graph is not uploaded due to size, but the same experiment can be reproduced on the subgraph, which contains the target ventricle of low conductance.

### Testing

1. Download and install Gurobi.
2. Run TestSimpleLocal.m

### Datasets

3. BrainSubgraph.mat: adjacency matrix for a subgraph which contains the target ventricle

### Dependencies

**General Purpose Min-Cut Solver**

The SimpleLocal function relies on a general min-cut/max-flow solver implemented in Gurobi (GurobiMaxFlow.m). An academic lisence can be obtained free of charge at http://www.gurobi.com/.

Alternatively, this can be replaced with any general purpose min-cut/max-flow routine that takes in an adjacency matrix of a graph with a source (1st node) and sink (last node), and outputs a min-cut indicator vector.