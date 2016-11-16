load('power_eqs_ar.sage')
load('DirGraph.sage')
load('elimVars.sage')
G = Graph()
G.add_edges([(0,1),(1,2),(1,3),(1,4),(2,5),(3,6)])
A = G.adjacency_matrix()
adj_mx = G.adjacency_matrix()
B = matrix(ZZ, [[  0,  -1,  -3, -11, -11, -12,  -2], [ -1,   0,  -3,  -2,-10,  -4, -11], [ -3,  -3,   0,  -5,  -4,  -9,  -4], [-11,  -2,  -5,   0, -6,  -4,  -9], [-11, -10,  -4,  -6,   0,  -6,  -9],[-12,  -4,  -9,  -4,  -6,   0,  -9],[ -2, -11,  -4,  -9,  -9 , -9,   0]])
Q = [0,1,2,3,4,5,6]
P = [0,0,0,0,0,0,0]
P_values = P
dirgraph = power_eqs_ar(G,Q,B,P=P_values)
final_eqs = elimVars(dirgraph,dirgraph.nodes[0])
solutions = alt_solve(dirgraph,dirgraph.nodes[0])