import argparse
import math

def pv_eqs(G):
	n = len(G)
	A = G.adjacency_matrix()
	x_names = ['x'+str(i) for i in range(1,n)]
	y_names = ['y'+str(i) for i in range(1,n)]
	b_names = []
	V_names = []
	for i in range(n):
		for j in range(i,n):
			if A[i,j] != 0:
				b_names.append('b'+str(i)+str(j))
		if i>0:
			V_names.append('V'+str(i))

	S = PolynomialRing(QQ,2*(n-1)+len(b_names)+len(V_names),names=x_names+y_names+b_names+V_names)
	v = S.variable_names()
	B = Matrix(S,n,n)
	V = []
	b = v[2*(n-1):]
	t = 0
	for i in range(n):
		for j in range(i,n):
			if j!=n:
				if A[i,j] != 0:
					bij = b[t]
					B[i,j] = S(bij)
					B[j,i] = S(bij)
					t+=1
	for i in range(len(V_names)):
		Vi = b[t]
		V.append(S(Vi))
		t+=1

	x = [1]+list(S.gens()[0:n-1])
	y = [0]+list(S.gens()[n-1:2*(n-1)])

	f = []
	h = []
	for i in range(1,n):
		p_eq = 0
		pv_eq = x[i]**2+y[i]**2-V[i]
		for j in range(n):
			if A[i,j] != 0:
				p_eq += B[i,j]*(x[j]*y[i]-x[i]*y[j])
		f.append(p_eq)
		h.append(pv_eq)

	eqs = f+h
	return eqs, S, (x,y,B,V)


parser = argparse.ArgumentParser()
parser.add_argument('-iters', type=int, dest="iters", default = [10000])
parser.add_argument('-s', nargs=9, type=float, dest="start_pt", default = [-50,-100,-100,-100,-100,-100,.1,.1,.1])
parser.add_argument('-r', type=float, dest="search_range", default = 1)
parser.add_argument('-v',action='store_true', default=False, dest='verbose')

args = vars(parser.parse_args())
trials = args["iters"]
search_range = args["search_range"]
verbose = args["verbose"]

G = graphs.CompleteGraph(4)
n = len(G)
eqs, S, vars = pv_eqs(G)
x_vars = vars[0]
y_vars = vars[1]
xy_vars_indet = x_vars[1:]+y_vars[1:]
B = vars[2]
b_vars = []
for i in range(n):
	for j in range(i+1,n):
		if B[i,j] != 0: b_vars.append(B[i,j])
C.<x> = ComplexField(100)[]

f_eqs = eqs[0:3]
h_eqs = eqs[3:]
elim_vars = xy_vars_indet[0:3]+xy_vars_indet[4:6]

start_pt = args["start_pt"]
start_d = {}
for i in range(6):
	start_d[b_vars[i]] = start_pt[i]
for i in range(3):
	start_d[P[i]] = start_pt[6+i]

print "Searching for solutions...\n"

g_eqs = [f.subs(start_d) for f in f_eqs]
p = S((S.ideal(g_eqs+h_eqs)).elimination_ideal(elim_vars).gens()[0])
c = p.coefficients()
s = 0
degr = len(c)-1
for i in range(degr+1):
	s += c[i]*x**(degr-i)
sol = s.roots()
max_wrong_l1 = 0
for t in sol:
	if t[0].real_part() < 0:
		l1 = abs(t[0].real()) + abs(t[0].imag())
		if l1 > max_wrong_l1:
			max_wrong_l1 = l1
	if t[0].real_part() >= 0:
		l1 = abs(t[0].imag())
		if l1 > max_wrong_l1:
			max_wrong_l1 = l1
if degr < 6:
	max_wrong_l1 = 1000000000

min_pt = start_d
min_dist = float(max_wrong_l1)

curr_pt = start_d
curr_dist = float(max_wrong_l1)

for T in range(trials):
	d = {}
	for i in range(len(b_vars)):
		d[b_vars[i]] = curr_pt[b_vars[i]] + uniform(-search_range, search_range)
	d[P[0]] = curr_pt[P[0]] + uniform(-search_range, search_range)
	for i in range(1,len(P)):
		d[P[i]] = curr_pt[P[i]] 

	g_eqs = [f.subs(d) for f in f_eqs]
	p = S((S.ideal(g_eqs+h_eqs)).elimination_ideal(elim_vars).gens()[0])

	c = p.coefficients()
	s = 0
	degr = len(c)-1
	for i in range(degr+1):
		s += c[i]*x**(degr-i)
	sol = s.roots()

	max_wrong_l1 = 0
	num_real_sol = 0
	for t in sol:
		if abs(t[0].imag_part()) < 1e-10 and t[0].real_part() >= 0:
			num_real_sol += 1


		if t[0].real_part() < 0:
			l1 = abs(t[0].real()) + abs(t[0].imag())
			if l1 > max_wrong_l1:
				max_wrong_l1 = l1

		if t[0].real_part() >= 0:
			l1 = abs(t[0].imag())
			if l1 > max_wrong_l1:
				max_wrong_l1 = l1

	if degr < 6:
		max_wrong_l1 = 1000000000

	if max_wrong_l1 <= min_dist:
		min_dist = float(max_wrong_l1)
		min_pt = d

	if max_wrong_l1 <= curr_dist :
		accept_prob = 1
	else:
		accept_prob = math.e**(-(float(max_wrong_l1)-curr_dist+0.0)/(T+1.0))

	threshold = uniform(1,1)
	if accept_prob >= threshold:
		curr_pt = d
		curr_dist = max_wrong_l1
		if verbose:
			print "Updated point:"
			pt_str = "Current point :"
			for i in range(6): pt_str +=" "+str(curr_pt[b_vars[i]])
			print pt_str
			print "Current distance : " + str(curr_dist)
			print "Number of nonnegative real solutions : " + str(num_real_sol)
			print ""

print "roots : "+ str(sol) 
print "Number of nonnegative real solutions : " + str(num_real_sol)
print "Min distance point : " + str(min_pt)
print "Min distance : " + str(min_dist)








