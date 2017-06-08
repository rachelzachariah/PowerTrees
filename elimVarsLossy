#elimVars eliminates variables to give univariate equations for the tree-system
#one univariate equation for each child of the root.
#Graph - underlying adjacency matrix for buses in the grid.
#G+ jB - line impedances
#P_values - active power
#Q_values - reactive power

def elimVars(Graph,G,B,Q_values=None,P_values=None):
	dirgraph = power_eqs_abr(Graph,G,Q,B,Q_values,P_values) #initalize req structures
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
		a = dirgraph.a
		g_eq = node.q_eq
		r_eq = node.r_eq

		for child in node.children:
			j = child.label
			elim_eq = altElim(dirgraph,child)
			g_eq = g_eq.resultant(elim_eq,dirgraph.A[i,j])
			g_eq = g_eq/(gcd(g_eq.coefficients()))
			if i==1 and j==3:
				print g_eq

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

#Given a directed graph with some susceptances and P values, 
#this returns a matrix beta whose i,j entry is the value of beta_{i,j}
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
	a = dirgraph.a
	b = dirgraph.b
	G = dirgraph.G
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
		 	elim = S((R[j1]^degree)*temp.subs({R[i]:(a[i,j1]^2)/R[j1]}))
		 	
		else:
			#eliminate the alpha/beta_ij's first
			g1=S(node.q_eq)
			g2=S(node.p_eq)
			counter=0
			for child in node.children:
				j2 = child.label
				p = S(elimVarsinternal(child,dirgraph)[0])
				q = S(elimVarsinternal(child,dirgraph)[1])
				if counter==0: #the first alpha_ij can just be eliminated by subs
					Gc = G[i,j2]
					Bc = B[i,j2]
					g1 = S((Gc*g2-Bc*g1)/(Gc^2+Bc^2))
					g2 = S((Bc*g2+Gc*g1)/(Gc^2+Bc^2))
					degree1 = p.degree(a[i,j2])
					degree2 = q.degree(b[i,j2])
					#print 'eliminating',A[i,j2],S((B[i,j2]^degree)*p.subs({A[i,j2]: S((g/S(B[i,j2]))+A[i,j2])}))
					pnew = S(((Gc^2+Bc^2)^degree1)*p.subs({a[i,j2]: g1-a[i,j2], b[i,j2]: g2-b[i,j2]}))
					qnew = S(((Gc^2+Bc^2)^degree2)*q.subs({a[i,j2]: g1-a[i,j2], b[i,j2]: g2-b[i,j2]}))
					counter = 1
				else:#remaining alpha_ij needs resultants for elimination
					p1 = S(pnew.resultant(qnew,b[i,j2])) #whats the right choices here???
					p2 = S(qnew.resultant(q,b[i,j2]))
					p3 = S(p.resultant(q,b[i,j2]))
					pnew = S(p1.resultant(p2,a[i,j2]))	
					qnew = S(p2.resultant(p3,a[i,j2]))

			#eliminate the node's local R variable next
	        #,R[j1]*g.subs({R[i]:(A[i,j1]^2)/R[j1]})
			degreep = pnew.degree(R[i])
			degreeq = qnew.degree(R[i])
			pnew = S((R[j1]^degreep)*pnew.subs({R[i]:(a[i,j1]^2+b[i,j1]^2)/R[j1]}))
			qnew = S((R[j1]^degreeq)*qnew.subs({R[i]:(a[i,j1]^2+b[i,j1]^2)/R[j1]}))
		if node.parent.label == 0 :
			elimq = pnew.resultant(qnew,b[i,j2])
			elimr = U(elimq)
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
