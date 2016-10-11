#elimVars eliminates variables to give univariate equations for the tree-system
#one univariate equation for each child of the root.

def elimVars(G,Q,B,P_values=None):
	dirgraph = power_eqs_ar(G,Q,B,P=P_values) #initalize req structures
	root=dirgraph.nodes[0]
	elim=elimVarsinternal(root,dirgraph)
	return elim

def altElim(dirgraph,node):
	i = node.label
	if i==0:
		elim_eqs = [altElim(dirgraph,child) for child in node.children]
		return elim_eqs
	else:
		g_eq = node.q_eq
		h_eq = node.h_eq

		for child in node.children:
			j = child.label
			elim_eq = altElim(dirgraph,child)
			g_eq = g_eq.resultant(elim_eq,dirgraph.A[i,j])

		g_eq = g_eq.resultant(h_eq,R[i])
		return g_eq

#Converts a polynomial in a real multivariate ring to a real univariate polynomial
#It then finds the roots of said polynomial
def convert_poly(f):
	var = f.variables()[0]
	degr = f.degree()
	C.<x> = PolynomialRing(RR)
	g = f.constant_coefficient()
	for i in range(1,degr+1):
		g = g + CC(f.coefficient(var^i))*x^i
	return g.roots()

def beta_vals(dirgraph):
	B = dirgraph.B
	n = B.ncols()
	beta = Matrix(RR,n,n)
	for i in range(1,n):
		curr_v = dirgraph.nodes[i]
		parent_v = curr_v.parent
		j = parent_v.label
		P_sum = curr_v.p_value
		for child in curr_v.children:
			P_sum += child.p_value
		beta[i,j] = P_sum/(-B[i,j])
		beta[j,i] = -beta[i,j]
	return beta


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
