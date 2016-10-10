# PowerTrees
DirGraph is a data structure to store the tree graph and corresponding equations. 
Power_eq's_ar populates the equations in a dirgraph object given a graph object and matrix of impedances. (can be modified to only accept the impedance matrix). Assumes lossless zero power-injection, hence ar.
elimVars generates a univariate polynomials for each child of the root that can be used to find all the solutions.
