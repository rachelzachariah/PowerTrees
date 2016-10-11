load('power_eqs_ar.sage')
load('DirGraph.sage')
load('elimVars.sage')
B = Matrix(5,5)
for i in range(5):
    for j in range(i+1,5):
        B[i,j] = B[j,i] = randint(-10,0)  
G = graphs.RandomTree(5)
Q = [0,1,2,3,4]
P = [1,3,5,2,6]
P_values = P
dirgraph = power_eqs_ar(G,Q,B,P=P_values)