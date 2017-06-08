#elimVars eliminates variables to give univariate equations for the tree-system
#one univariate equation for each child of the root.

def elimVars(dirgraph,root):
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

		g_eq = g_eq.resultant(r_eq,R[i])
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

#Converts a univariate polynomial in a real multivariate ring to a real univariate polynomial
#It then finds the roots of said polynomial and returns them.
#Here there is the added ability to specify the precision of the real field you want.
def convert_poly_precision(f,precision):
	var = f.variables()[0]
	degr = f.degree()
	R.<x> = PolynomialRing(RealField(precision))
	g = f.constant_coefficient()
	for i in range(1,degr+1):
		g = g + R(f.coefficient(var^i))*x^i
	return g

#Converts a real polynomial to a complex polynomial.
def convert_complex_poly(f):
	var = f.variables()[0]
	degr = f.degree()
	R.<x> = PolynomialRing(CC)
	g = f.constant_coefficient()
	for i in range(1,degr+1):
		g = g + CC(f.coefficient(var^i))*x^i
	return g

#Given a directed graph with the root node specified, this solves the power flow equations
def alt_solve(dirgraph,root,tol=1.0e-6):

	beta = dirgraph.beta

	#Sol amalg will have n arrays in it
	#It will have one array for each child of the root node
	#The ith element in sol_amalg will contain dictionaries of all the real solutions
	#to that part of the graph
	sol_amalg = []

	for node in root.children:
		p = node.label
		branch_solutions = []
		f = node.res_eq #This is a univariate polynomial in a multivariate ring
		g = convert_poly(f) #We convert it to a polynomial in a real univariate ring
		g_roots = g.roots() #These are all its real roots
		num_roots = len(g_roots)
		print "Node "+str(p)+" branch : "+str(num_roots)+" potential solutions"

		#We now go through all the real roots and see which will result in real solutions
		#to the entire system
		for i in range(num_roots):
			print "Checking solution " + str(i)

			#We create a dictionary and add the alpha and R values for the node
			#just adjacent to the root node
			sol_dict = {}
			sol_dict[f.variables()[0]] = g_roots[i][0]
			sol_dict[dirgraph.R[p]] = g_roots[i][0]^2+beta[root.label,node.label]^2

			to_visit = [node] #We perform a breadth-first search on the tree
			real_solution = True

			while len(to_visit) > 0:
				current_node = to_visit[0]
				to_visit.remove(current_node)
				for child in current_node.children:
					to_visit.append(child)
				#print "Visiting node " + str(current_node.label)

				#We update the solution dictionary to contain the corresponding solutions
				#to the variables below the node
				#If any solution is not real, this returns None
				sol_dict = update_sol(dirgraph,current_node,sol_dict,tol=tol)

				#Here we check if we returned none
				if sol_dict==None:
					real_solution = False
					to_visit = []

			#If we don't return none, then we have a set of real solutions to the system
			if real_solution:
				branch_solutions.append(sol_dict)
		sol_amalg.append(branch_solutions)

	n_branches = len(root.children)
	for l in range(n_branches):
		c = root.children[l].label
		print "Node "+str(c)+" branch : "+str(len(sol_amalg[l]))+" solutions found"

	return sol_amalg

#Given a directed graph and a node, with some values in sol_dict, finds 
#the corresponding solutions to the children of the node
def update_sol(dirgraph,node,sol_dict,tol):

	#If the node has no children then we're done
	if len(node.children)==0: return sol_dict

	child_solutions = []
	child_vars = []

	p = node.label

	num_potential_solutions = 1

	for child in node.children:
		f_child = child.res_eq #This is the polynomial governing the child
		f_child_sub = f_child.subs(sol_dict) #We sub in the values of the variables above it
		g_child = convert_poly(f_child_sub) #We convert it to a real univariate polynomial
		child_roots = g_child.roots() #We find its roots
		child_roots = [s[0] for s in child_roots]
		child_solutions.append(child_roots)
		num_potential_solutions = num_potential_solutions*len(child_roots)

	#We check if the solutions to each of the children satisfy the necessary equation
	curr_q_eq = (node.q_eq).subs(sol_dict)
	num_children = len(node.children)
	prod = product(*child_solutions) #This contains the solutions to the children
	found = False

	num_solutions_checked = 0

	while found == False and num_solutions_checked < num_potential_solutions:
		num_solutions_checked += 1
		#print "num sol checked = " + str(num_solutions_checked)
		child_roots = prod.next()
		child_solutions_dict = {}
		for cc in range(num_children):
			c = node.children[cc].label
			child_solutions_dict[dirgraph.A[p,c]] = child_roots[cc]
		sub_q_eq = curr_q_eq.subs(child_solutions_dict) #We sub in the current children solutions

		#We have approximate values, so we heck if they satisfy the equation up to some tolerance
		if abs(RR(sub_q_eq)) <= tol:
			found = True
			for cc in range(num_children):
				child = node.children[cc]
				c = child.label
				sol_dict[dirgraph.A[p,c]] = child_roots[cc]
				sol_dict[dirgraph.R[c]] = (child_roots[cc]^2+dirgraph.beta[p,c]^2)/sol_dict[dirgraph.R[p]]
			return sol_dict
		
	#If none of the solutions work, return None		
	return None

#Converts solutions to abr in to solutions for xy
def convert_abr_xy(dirgraph,ar_solutions):
	beta = dirgraph.beta
	A = dirgraph.A
	R = dirgraph.R

	xy_solutions = []
	root = dirgraph.nodes[0]
	n_branches = len(root.children)
	for b in range(n_branches):
		branch_top = root.children[b]
		xy_branch_solutions = []
		for l in range(len(ar_solutions[b])):
			ar_sol_dict = ar_solutions[b][l]
			xy_sol_dict = {}
			xy_sol_dict['x0'] = 1
			xy_sol_dict['y0'] = 0
			to_visit = [branch_top]

			while len(to_visit) > 0:
				current_node = to_visit[0]
				to_visit.remove(current_node)
				for child in current_node.children:
					to_visit.append(child)
				c = current_node.label
				p = current_node.parent.label

				eq_matrix = Matrix(RR,2,2)
				eq_matrix[0,0] = xy_sol_dict['x'+str(p)]
				eq_matrix[0,1] = xy_sol_dict['y'+str(p)]
				eq_matrix[1,0] = -xy_sol_dict['y'+str(p)]
				eq_matrix[1,1] = xy_sol_dict['x'+str(p)]

				if det(eq_matrix) == 0:
					xy_sol_dict['x'+str(c)] = 0
					xy_sol_dict['y'+str(c)] = 0
				else:
					a_b_vals = Matrix(RR,2,1)
					a_b_vals[0,0] = ar_sol_dict[A[c,p]]
					a_b_vals[1,0] = beta[c,p]
					current_xy = eq_matrix.inverse()*a_b_vals
					xy_sol_dict['x'+str(c)] = current_xy[0,0]
					xy_sol_dict['y'+str(c)] = current_xy[1,0]
			xy_branch_solutions.append(xy_sol_dict)
		xy_solutions.append(xy_branch_solutions)
	return xy_solutions	


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
			node.res_eq = elimq
		else:
			node.res_eq = elim
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
