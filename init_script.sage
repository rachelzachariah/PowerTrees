load('power_eqs_abr.sage')
load('DirGraph.sage')
load('elimVarsLossy.sage')
graph = Graph()
graph.add_edges([(0,1),(1,2),(1,3)])
A = graph.adjacency_matrix()
B = Matrix(4,4)
G = Matrix(4,4)
for i in range(4):
    for j in range(i+1,4):
    	if A[i,j] != 0:
           	B[i,j] = B[j,i] = randint(-10,-1)  
        	G[i,j] = G[j,i] = randint(-10,-1) 
Q = [0,0,0,0,1,2,3,4,5,6]
P = [0,0,0,0,0,0,0]
dirgraph = power_eqs_abr(graph,G,B,Q,P)