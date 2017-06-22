import argparse

def solve_eqs(dirgraph):
	final_eqs = altElim(dirgraph,dirgraph.nodes[0])
	solutions = alt_solve(dirgraph,dirgraph.nodes[0],tol=1.0e-3)
	return solutions

def write_equations(eqs,filename):
	disp = str(len(eqs))+'\n'
	for coeff in eqs:
		str_coeff = str(coeff)
		disp = disp+str_coeff.replace(' ','') + ';\n'
	f = open(filename+'_eqs.txt','w')
	f.write(disp)
	f.close()

parser = argparse.ArgumentParser()
parser.add_argument('-reps', type=int, dest="reps", default = 10)
parser.add_argument('-n', type=int, dest="n", default = 6)
parser.add_argument('-md', type=int, dest="md", default = 5)
parser.add_argument('-w', type=str, dest="w", default = "")

args = vars(parser.parse_args())
reps = args["reps"]
n = args["n"]
max_deg = args["md"]
filename = args["w"]

load('power_eqs_ar.sage')
load('DirGraph.sage')
load('elimVars.sage')
#load('write_eqs.sage')

#for n in range(4,11):
times = []

curr_time = 0
for t in range(reps):
	found = False
	while found == False:
		G = graphs.RandomTree(n)
		if max(G.degree()) <= max_deg:
			found = True	

	gd = G.degree()
	max_gd = max(gd)
	best_v = gd.index(max_gd)
	if best_v != 0:
		G.relabel({0:best_v, best_v:0})

	A = G.adjacency_matrix()
	B = Matrix(ZZ,n)
	P_values = []
	Q = []
	for i in range(n):
		for j in range(i+1,n):
			if A[i,j] != 0:
				B[i,j] = B[j,i] = randint(-50,-1)
		P_values.append(randint(1,10))
		Q.append(randint(1,10))
	dirgraph = power_eqs_ar(G,Q,B,P=P_values)
	if filename != "":
		write_equations(dirgraph.eqs,filename+"t"+str(t))
	curr_time += timeit('solve_eqs(dirgraph)',seconds=True,repeat=1)
curr_time = curr_time/reps
print "Average time for "+str(n)+" nodes : "+str(curr_time)
