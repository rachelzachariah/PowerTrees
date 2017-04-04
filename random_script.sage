import argparse

def solve_eqs(dirgraph):
	final_eqs = altElim(dirgraph,dirgraph.nodes[0])
	solutions = alt_solve(dirgraph,dirgraph.nodes[0],tol=1.0e-3)
	return solutions

parser = argparse.ArgumentParser()
parser.add_argument('-reps', type=int, dest="reps", default = 10)
parser.add_argument('-n', type=int, dest="n", default = 3)
parser.add_argument('-N', type=int, dest="N", default = 6)

args = vars(parser.parse_args())
reps = args["reps"]
nmin = args["n"]
nmax = args["N"]


load('power_eqs_ar.sage')
load('DirGraph.sage')
load('elimVars.sage')
#load('write_eqs.sage')

#for n in range(4,11):
times = []
for n in range(nmin,nmax+1):
	curr_time = 0
	for t in range(reps):
		G = graphs.RandomTree(n)
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
		curr_time += timeit('solve_eqs(dirgraph)',seconds=True,repeat=1)
	curr_time = curr_time/reps
	print "Average time for "+str(n)+" nodes : "+str(curr_time)
