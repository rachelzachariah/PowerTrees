#This function takes as input a graph G and a list Q of Q_i values
#It returns the equations in terms of alpha_{i,j} and R_i
import numpy
def power_eqs_ar(G,Q,B,P=None):
	adj_mx = G.adjacency_matrix() #We compute the adjacency matrix of G
	n = len(G) #This is the number of buses
	dirgraph = DirGraph() 

	#This creates the names of variables we will use
	#For example, a might equal ['a01','a12','a23'] if G is a path graph
	#aij corresponds to alpha_{i,j} in our equations
	#bij are the susceptances
	#Ri are the R_i = x_i^2 + y_i^2

	a_names = []
	#b_names = []
	R_names = ['R'+str(i) for i in range(1,n)]

	for i in range(n):
		for j in range(i+1,n):
			if adj_mx[i,j] != 0:
				a_names.append('a'+str(i)+str(j))
				#b_names.append('b'+str(i)+str(j))


	#We now create a polynomial ring with the given variables above
	S = PolynomialRing(QQ,len(a_names)+n-1,names=a_names+R_names)
	U = PolynomialRing(RR,len(a_names)+n-1,names=a_names+R_names)
	v = S.gens_dict() #This v takes as input the variable names and outputs the variable objects


	#We now construct matrices of the variable objects for the a, b
	#We also construct a list of the Rs
	#Note that R[0] = 1 since R_0 = 1.
	#Now A[i,j] = a_{i,j} and B[i,j] = b_{i,j}

	A = Matrix(S,n,n)
	#B = Matrix(S,n,n)
	R = [1]+[v['R'+str(i)] for i in range(1,n)]

	for i in range(n):
		for j in range(i+1,n):
			if adj_mx[i,j] != 0:
				A[i,j] = A[j,i] = v['a'+str(i)+str(j)]
				#B[i,j] = B[j,i] = v['b'+str(i)+str(j)]

    #Create a rooted tree object where every node stores the r_eq(between node and parent),
    #the node's own q_eq, and A,B,R,S

	dirgraph.A = A
	dirgraph.B = B
	dirgraph.R = R
	dirgraph.S = S
	dirgraph.U = U
	root = Node(0,0,0)
	dirgraph.add_node(root)
	populate_node(root,dirgraph,adj_mx,P=P)


	beta = beta_vals(dirgraph) #This solves for the beta_{i,j} values

	q_eqs = [] #These will be the Q equations
	r_eqs = [] #These equations give the relations a_{i,j}^2 = R_iR_j

	for i in range(1,n): #Note that i ranges over [1,2,...,n-1]
		#We will create the Q equation at bus i
		#We exclude 0 as a value for i since 0 is relaxed and has no Q equation
		q_eq = Q[i]

		#Note that j ranges over [0,1,...,n-1]
		#We will now loop over all other nodes (including 0) to see if there is an edge (i,j)
		for j in range(n):

			#We now test to see if i has an edge to j
			if adj_mx[i,j] != 0:

				#This adds the relevant term to the Q equation
				q_eq += B[i,j]*(R[i]-A[i,j])

				#This creates the a_{i,j}^2 = R_iR_J equation
				#We restrict to j < i so we don't count the same equation twice
				#Note that since R[0] = 1, when j = 0 and (i,j) is an edge
				#we get the equation A[i,j]**2 - R[i]*1
				if j < i:
					#We add the R equation to the list of R equations
					r_eq = A[i,j]**2+beta[i,j]**2-R[i]*R[j]
					dirgraph.nodes[i].r_eq = r_eq
					r_eqs.append(r_eq)

		#We've now added all relevant parts to the Q equation, so we add it to the list of Q equations		
		q_eqs.append(q_eq)

	#This concatenates the lists to create all the equations
	eqs = q_eqs + r_eqs

	#Add the equations to dirgraph
	dirgraph.eqs = eqs 

	return dirgraph



#If you want to use this method, you might do the following in sage:

#load('sage_example.sage') 
#G = graphs.PathGraph(3)
#Q = [1,1,1]
#eqs, S, vars = power_eqs_ar(G,Q)

#If I want to compute an ideal, I could do:
#I = S.ideal(eqs)

#If I want to plug in some susceptances, I could do:
#susceptances = {B[0,1]:1, B[1,2]:2}
#J = I.subs(susceptances)

#J is now the ideal with B[0,1] = 1, B[1,2] = 2
#If you'd like to extract the equations of J, you can do

#J_eqs = J.gens()


#populate_node adds nodes beginning with "node" and recursively fills in their parents, children 
# Q-eqns and R-eqns based on adjacency_matrix.

def populate_node(node, dirgraph, adjacency_matrix, P=None):
	A = dirgraph.A
	B = dirgraph.B
	R = dirgraph.R
	S = dirgraph.S
	U = dirgraph.U
	i = node.label
	q=0
	#the node's parent (if any) is dealt with first
	if node.has_parent():
		parent = node.parent  
		j = parent.label 
		q = B[i,j]*(R[i]- A[i,j]) + Q[i] #parent's contribution to Q-eqn 
		node.r_eq = A[i,j]^2-R[i]*R[j] #input R-eqn between node and parent

	n=len(numpy.array(adjacency_matrix)) #number of nodes

	#now we deal with the children
	for j in range(n):
		if ((node.has_parent() and j != parent.label) or not node.has_parent()) and adjacency_matrix[i,j]!=0:
			if node.has_parent(): 
				q += B[i,j]*(R[i]- A[i,j]) #update children's contibution to node's Q-eqn
				                           #root has no Q-eqn, it is set to 0 when initialized 

			#initialize a new child node and populate the necessary fields	                           
			chld = Node(0,0,j) 
			chld.add_parent(node) 
			dirgraph.add_node(chld) #update the node-holder
			node.add_child(chld) #update the node
			populate_node(chld,dirgraph,adjacency_matrix,P=P)
	#set node's Q-eqn that has been built now
	node.q_eq=q #set node's Q-eqn that has been built now


	#Add a P value to the node
	#If no P values are input, it defaults to lossless
	if P:
		node.p_value = P[i]
	else:
		node.p_value = 0

	return


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



