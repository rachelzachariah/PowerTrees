# class DirGraph is a node holder.

class DirGraph:

	def __init__(self,nodes=dict()):
		self.nodes = nodes

	def add_node(self,node):
		self.nodes[node.label] = node

# Node class carries the bus's Q-eqn and R-eqn (between node and parent).
# It keeps track of its parent and children if any.  

class Node:

	def __init__(self,r_eq,q_eq,label):
		self.r_eq = r_eq  #R-eqn between node and parent
		self.q_eq = q_eq   #Q-eqn at node 
		self.label = label  #number of the bus
		self.children = []	 #children is initialized to an empty list	
		self.parent = [] #parent is intialized to an empty list

	def add_parent(self,node): #Adds parent node to the node's "parent" field
		self.parent = node  

	def add_child(self,child): #Adds a child node to the list of nodes in "children" field
		self.children.append(child)

	def has_parent(self):  #node can be queried about if it is root or not
		if not self.parent:
			return false
		else:
			return true





