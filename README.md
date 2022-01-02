This is a django based application for computing disjoint pair of least-cost paths separated by a minimum distance.
We realized that in hilly areas where the ground is higly undulating, the available least-cost path routing software
packages may give suboptimal soulution as they usually find two consecutive least-cost paths separated by the minimum
distance which may be more expensive than some combination of two paths. So we minimize the combined cost of paths 
using a heuristics 'Optimal transmission network topology for resilient power supply' developed by \cite{Zinchenko Y. et al. 2016}.
https://github.com/jha-amit/DS2P/issues/1#issue-1091936540

The algorithm uses the following set up and processes:
1. The project area is represented as a diamond shape biconnected digraph and two such orthogonal graphs are combine as a 3D graph using affine transformations.
   This 3D embedding only allows eligible combination of nodes from 2D graph that provide embedded separation constraint and biconnectivity constraint in the 3D graph edges.
   
2. A single shortest path is computed in the 3D which overcomes the computational intractibility of finding disjoint paths with constraints on a graph in polynomial time.
3. This single least-cost path in 3D is projectes back to 2D grpahs. Following images may help to clear out the process. 
The following instructions are to be followed to run the application
