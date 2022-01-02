This is a django based application for computing disjoint pair of least-cost paths separated by a minimum distance.
We realized that in hilly areas where the ground is higly undulating, the available least-cost path routing software
packages may give suboptimal soulution as they usually find two consecutive least-cost paths separated by the minimum
distance which may be more expensive than some combination of two paths. So we minimize the combined cost of paths 
using a heuristics 'Optimal transmission network topology for resilient power supply' developed by \cite{Zinchenko Y. et al. 2016}.


The algorithm uses the following set up and processes:
1. The project area is represented as a diamond shape biconnected digraph and two such orthogonal graphs are combine as a 3D graph using affine transformations.
   This 3D embedding only allows eligible combination of nodes from 2D graph that provide embedded separation constraint and biconnectivity constraint in the 3D graph edges.
   ![2-Paths plan view](Images/snippet app.jpg?raw=true "Optional Title")
   
2. A single shortest path is computed in the 3D which overcomes the computational intractibility of finding disjoint paths with constraints on a graph in polynomial time.
3. This single least-cost path in 3D is projectes back to 2D grpahs. Following images may help to clear out the process.
   https://github.com/jha-amit/DS2P/issues/1#issue-1091936540
4. Near the terminals we introduced two complete graphs so that we can avoid the intractability due to path separation constraint. We relaxed the path separation constraint
   in the neighbourhood of terminals and computed disjoint paths. Then we computed two least cost paths on the two complete graphs that connect with the disjoint paths at points
   where path separation constraint is relaxed. We replace the portion of separation constraint relaxed disjoint path with the single least-cost paths on the complete graph. 
   
The following instructions are to be followed to run the application
1. Input lat-long-cost data in file provided as Cost_matrix.csv
2. Provide Endgraph 
