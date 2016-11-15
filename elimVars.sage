from itertools import product

#elimVars eliminates variables to give univariate equations for the tree-system
#one univariate equation for each child of the root.

def elimVars(G,Q,B,P_values=None):
	dirgraph = power_eqs_ar(G,Q,B,P=P_values) #initalize req structures
	root=dirgraph.nodes[0]
	elim=elimVarsinternal(root,dirgraph)
	return elim

#Given a directed graph and a node i, this recursively calls to get the variable reduced
#equations from all the node's children. It then uses resultants to eliminate the variables a_{i,j}
#for each child j. Finally, it eliminates R_i and returns the equations
#If the node passed is the root, is returns a list of the resulting equation for each child of the root
#This method also adds an attribute res_eq to each node in the dirgraph that displays what polynomial it passed
def altElim(dirgraph,node):
	i = node.label
	if i==0:
		elim_eqs = [altElim(dirgraph,child) for child in node.children]
		return elim_eqs
	else:
		R = dirgraph.R
		A = dirgraph.A
		g_eq = node.q_eq
		r_eq = node.r_eq

		for child in node.children:
			j = child.label
			elim_eq = altElim(dirgraph,child)
			g_eq = g_eq.resultant(elim_eq,dirgraph.A[i,j])
			g_eq = g_eq/(gcd(g_eq.coefficients()))

		g_eq = g_eq.resultant(r_eq,dirgraph.R[i])
		g_eq = g_eq/(gcd(g_eq.coefficients()))
		node.res_eq = g_eq
		return g_eq

#Converts a polynomial in a real multivariate ring to a real univariate polynomial
#It then finds the roots of said polynomial and returns them
def convert_poly(f):
	var = f.variables()[0]
	degr = f.degree()
	R.<x> = PolynomialRing(RR)
	g = f.constant_coefficient()
	for i in range(1,degr+1):
		g = g + RR(f.coefficient(var^i))*x^i
	return g

def alt_solve(dirgraph,root,tol=1.0e-6):

	beta = dirgraph.beta
	sol_amalg = []

	for node in root.children:
		p = node.label
		branch_solutions = []
		f = node.res_eq
		g = convert_poly(f)
		g_roots = g.roots()
		num_roots = len(g_roots)
		for i in range(num_roots):
			print "Finding solution " + str(i)
			sol_dict = {}
			sol_dict[f.variables()[0]] = g_roots[i][0]
			sol_dict[dirgraph.R[p]] = g_roots[i][0]^2+beta[root.label,node.label]^2

			to_visit = [node] #We perform a breadth-first search on the tree
			for child in node.children:
				to_visit.append(child)

			real_solution = True

			while len(to_visit) > 0:
				current_node = to_visit[0]
				print "Visiting node " + str(current_node.label)
				to_visit.remove(current_node)

				sol_dict = update_sol(dirgraph,current_node,sol_dict,tol=tol)
				if sol_dict==None:
					real_solution = False
					to_visit = []

			if real_solution:
				branch_solutions.append(sol_dict)
			print ""
		sol_amalg.append(branch_solutions)
	return sol_amalg

def update_sol(dirgraph,node,sol_dict,tol):

	if len(node.children)==0: return sol_dict

	child_solutions = []
	child_vars = []

	p = node.label

	num_potential_solutions = 1

	for child in node.children:
		f_child = child.res_eq
		f_child_sub = f_child.subs(sol_dict)
		g_child = convert_poly(f_child_sub)
		child_roots = g_child.roots()
		child_roots = [s[0] for s in child_roots]
		child_solutions.append(child_roots)
		num_potential_solutions = num_potential_solutions*len(child_roots)

	curr_q_eq = (node.q_eq).subs(sol_dict)
	num_children = len(node.children)
	prod = product(*child_solutions)
	found = False

	num_solutions_checked = 0
	print "num_potential_solutions =" + str(num_potential_solutions)

	while found == False and num_solutions_checked < num_potential_solutions:
		num_solutions_checked += 1
		print "num sol checked = " + str(num_solutions_checked)
		child_roots = prod.next()
		child_solutions_dict = {}
		for cc in range(num_children):
			c = node.children[cc].label
			child_solutions_dict[dirgraph.A[p,c]] = child_roots[cc]
		sub_q_eq = curr_q_eq.subs(child_solutions_dict)

		print abs(RR(sub_q_eq))

		if abs(RR(sub_q_eq)) <= tol:
			found = True
			for cc in range(num_children):
				child = node.children[cc]
				c = child.label
				sol_dict[dirgraph.A[p,c]] = child_roots[cc]
				sol_dict[dirgraph.R[c]] = (child_roots[cc]^2+dirgraph.beta[p,c]^2)/sol_dict[dirgraph.R[p]]
			return sol_dict
			
	return None


#elimVarsinternal is recursively called to eliminate the nodes R variable and 
#alpha_ij with its children. 

def elimVarsinternal(node,dirgraph):
	S = dirgraph.S  
	A = dirgraph.A
	B = dirgraph.B
	R = dirgraph.R
	U = dirgraph.U
	children= node.children
	i = node.label
	if node.has_parent(): #non-root nodes
		j1 = node.parent.label 
		if not children:  #leaf node
			#eliminate the nodes local R variable from the Q-eqn
			#, R[j1]*node.q_eq.subs({R[i]:(A[i,j1]^2)/R[j1]})
		 	degree = node.q_eq.degree(R[i])
		 	temp = S(node.q_eq)
		 	elim = S((R[j1]^degree)*temp.subs({R[i]:(A[i,j1]^2)/R[j1]}))
		 	
		else:
			#eliminate the alpha_ij's first
			g=S(node.q_eq)
			counter=0
			for child in node.children:
				j2 = child.label
				p = S(elimVarsinternal(child,dirgraph))
				if counter==0: #the first alpha_ij can just be eliminated by subs
					degree = p.degree(A[i,j2])
					#print 'eliminating',A[i,j2],S((B[i,j2]^degree)*p.subs({A[i,j2]: S((g/S(B[i,j2]))+A[i,j2])}))
					g = S((B[i,j2]^degree)*p.subs({A[i,j2]: (g/S(B[i,j2]))+A[i,j2]}))
					counter = 1
				else:#remaining alpha_ij needs resultants for elimination
					g = S(g.resultant(p,A[i,j2]))

			#eliminate the node's loval R variable next
	        #,R[j1]*g.subs({R[i]:(A[i,j1]^2)/R[j1]})
			degree = g.degree(R[i])
			elim = S((R[j1]^degree)*g.subs({R[i]:(A[i,j1]^2)/R[j1]}))
		if node.parent.label == 0 :
			elimq = elim
			elimr = U(elim)
			elim = []
			elim.append(elimq)
			elim.append(elimr)	
	else: #root node
		elim = []
		elimQ = [] #list to hold univariate eqns from all root's children
		elimR = []
		for child in node.children:
			tempelim = elimVarsinternal(child,dirgraph)
			elimQ.append(tempelim[0])
			elimR.append(tempelim[1])
			#print 'child:',elimVarsinternal(child,dirgraph,A,B,R,S)
		elim.append(elimQ)
		elim.append(elimR)
	return elim 
