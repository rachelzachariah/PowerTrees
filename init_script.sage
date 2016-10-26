load('power_eqs_ar.sage')
load('DirGraph.sage')
load('elimVars.sage')
G = Graph()
G.add_edges([(0,1),(1,2),(1,3),(1,4),(2,5),(3,6)])
A = G.adjacency_matrix()
adj_mx = G.adjacency_matrix()
B = Matrix(7,7)
for i in range(7):
    for j in range(i+1,7):
    	if A[i,j] != 0:
        	B[i,j] = B[j,i] = randint(-10,-1)  
Q = [0,1,2,3,4,5,6]
P = [0,0,0,0,0,0,0]
P_values = P
dirgraph = power_eqs_ar(G,Q,B,P=P_values)