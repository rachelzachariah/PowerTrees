#elimVars eliminates variables to give univariate equations for the tree-system
#one univariate equation for each child of the root.

def elimVars(G,Q,B):
	dirgraph = power_eqs_ar(G,Q,B) #initalize req structures
	root=dirgraph.nodes[0]
	elim=elimVarsinternal(root,dirgraph)
	return elim

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
