import argparse
import signal
from multiprocessing import Process
from time import sleep

def handler(signum, frame):
	raise Exception("Process is taking too long")

def find_ideal(dirgraph):
	final_eqs = altElim(dirgraph,dirgraph.nodes[0])
	return final_eqs

def grob_ideal(dirgraph):
	S = dirgraph.S
	v = S.gens()
	eqs = dirgraph.eqs
	m = len(dirgraph.nodes[0].children)
	I = S.ideal(eqs)
	J = I.elimination_ideal(v[m:])
	return J

def solve_eqs(dirgraph):
	final_eqs = altElim(dirgraph,dirgraph.nodes[0])
	solutions = alt_solve(dirgraph,dirgraph.nodes[0],tol=1.0e-3)
	return solutions

def write_equations(eqs,filename):
	disp = str(len(eqs))+'\n'
	for coeff in eqs:
		str_coeff = str(coeff)
		disp = disp+str_coeff.replace(' ','') + ';\n'
	f = open('phc_eqs/'+filename+'_eqs.txt','w')
	f.write(disp)
	f.close()

def bounded_tree(n,max_deg):
	G = Graph(n)
	for i in range(1,n):
		adm_verts = [j for j in range(i) if G.degree(i) < max_deg]
		r = randint(0,len(adm_verts)-1)
		k = adm_verts[r]
		G.add_edge(i,k)
		G.add_edge(k,i)
	return G

parser = argparse.ArgumentParser()
parser.add_argument('-reps', type=int, dest="reps", default = 10)
parser.add_argument('-n', type=int, dest="n", default = 6)
parser.add_argument('-md', type=int, dest="md", default = 5)
parser.add_argument('-a', type=int, dest="alarm_val", default=60)
parser.add_argument('-w', type=str, dest="w", default = "")
parser.add_argument('-g',action='store_true', default=False, dest='grobner')
parser.add_argument('-u',action='store_true', default=False, dest='grobner_undo')

args = vars(parser.parse_args())
reps = args["reps"]
n = args["n"]
alarm_val = args["alarm_val"]
max_deg = args["md"]
filename = args["w"]
grob = args["grobner"]
gu = args["grobner_undo"]

load('power_eqs_ar.sage')
load('DirGraph.sage')
load('elimVars.sage')
#load('write_eqs.sage')

#for n in range(4,11):
times = []
grob_times = []

curr_time = 0
bad = []
bad_grob = []

for t in range(reps):
	G = bounded_tree(n,max_deg)
	# found = False
	# while found == False:
	# 	G = graphs.RandomTree(n)
	# 	if max(G.degree()) <= max_deg:
	# 		found = True	

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
		write_equations(dirgraph.eqs,filename+"_"+str(t))

	if grob:
		signal.signal(signal.SIGALRM, handler)	
		signal.alarm(alarm_val)
		curr_time = -1
		grob_time = -1
		try:
			curr_time = timeit('find_ideal(dirgraph)',seconds=True,number=1)
		except Exception, exc:
			bad.append(t)
		signal.alarm(0)
		times.append(curr_time)
		
		if gu==False:
			signal.signal(signal.SIGALRM, handler)	
			signal.alarm(alarm_val)
			try:
				grob_time = timeit('grob_ideal(dirgraph)',seconds=True,number=1)
			except Exception, exc:
				bad_grob.append(t)	
			signal.alarm(0)			
			grob_times.append(grob_time)

	else:
		signal.signal(signal.SIGALRM, handler)	
		signal.alarm(alarm_val)
		curr_time = -1
		try:		
			curr_time = timeit('solve_eqs(dirgraph)',seconds=True,number=1)
		except Exception, exc:
			bad_grob.append(t)	
		signal.alarm(0)	
		times.append(curr_time)


if grob:
	print str(n)+" nodes elim : "+str(times)
	print "Failed executions : " + str(bad)
	if gu==False:
		print str(n)+" nodes Grobner : "+str(grob_times)
		print "Failed Grobner executions : "+str(bad_grob)
else:
	print str(n)+" nodes : "+str(times)
	print "Failed executions : " + str(bad)
