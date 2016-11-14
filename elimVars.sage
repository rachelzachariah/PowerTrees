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
	print "i = "+str(i)
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
			print "j = " + str(j)
			elim_eq = altElim(dirgraph,child)
			g_eq = g_eq.resultant(elim_eq,dirgraph.A[i,j])
			g_eq = g_eq/(gcd(g_eq.coefficients()))
			if i==1 and j==3:
				print g_eq

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

def alt_solve(dirgraph,node):
	final_eqs = altElim(dirgraph,node)

	for f in final_eqs:
		g = convert_poly(f)
		g_roots = g.roots()
		num_roots = len(g_roots)
		solutions = []
		for i in range(num_roots):
			d = {}
			d[f.variables()[0]] = g_roots[i]


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
