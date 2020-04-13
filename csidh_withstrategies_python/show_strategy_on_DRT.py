# -----------------------------------------------------------------------------
# Multiplicative strategy written as a graph (a triangle à la SIDH)

# Python-code implementation perfomed by 
# Jesús-Javier Chi-Domínguez <jjchi@computacion.c.cinvestav.mx, chidoys@gmail>
#                                Computer science department, Cinvestav - IPN.
# -----------------------------------------------------------------------------

from math import sqrt
import sys

import matplotlib.pyplot as plt
import networkx as nx

def DRT(n):

	vertexes = dict()      # list of the position of each node
	vertex_colors = []    # color of each node: red for the leaves, otherwise color is set to white
	acc = 0

	# Different shape of the isogeny graph
	for i in range(n):
		for j in range(n - 1 - i):
			vertex_colors.append('black')
			vertexes[acc] = (i, -j)
			acc += 1

	
	return vertexes, vertex_colors, [], []

def strategy_evaluation(strategy, n):

	vertexes = dict()      # list of the position of each node
	edges    = []          # edges of the isogeny triangle

	edge_colors   = []    # color of each edge: blue for scalar multiplications, and orange for isogeny evalutions
	vertex_colors = []    # color of each node: red for the leaves, otherwise color is set to white

	vertexes[0] = (0.0,0.0)	# Root of the strategy
	ramifications = [] 	# nodes having two different types of edges
	moves = [0] 		# 
	k = 0       		# current element of the strategy
	t = 0				# current vertex

	# Strategy evaluation starts next
	ramifications.append([0, vertexes[0]])							# The root is always a ramification
	for i in range(len(strategy)):
	    
	    # Getting the new possible ramifications
	    while sum(moves) < (n - 1 - i):

	    	vertex = ramifications[-1][1]

	    	moves.append(strategy[k])					# Increasing moves
	    	k += 1                                  	# Next element of the strategy

	    	t += 1										# Next vertex
	    	edges.append( (t - 1, t) )					# New edge to be added
	    	edge_colors.append('tab:blue')					# Color of this type of edge is always blue

	    	vertexes[t] = (i, vertex[1] - strategy[k-1])	# New vertex to be added
	    	ramifications.append([t, vertexes[t]])			# New ramification
	    	vertex_colors.append('black')
	    
	    # Removing the ramifications (it is not required more!)
	    ramifications.pop()
	    moves.pop()
	    
	    # Updating the ramifications
	    for j in range(len(ramifications)):
	    	t += 1
	    	vertexes[t] = ( ramifications[j][1][0] + 1.0, ramifications[j][1][1] )
	    	edges.append( (ramifications[j][0], t))

	    	ramifications[j] = [t, vertexes[t]]
	    	edge_colors.append('tab:red')
	
	return vertexes, vertex_colors, edges, edge_colors
    


BATCHES = lambda n, sigma: [ len([ j for j in range(i, n, sigma) ]) for i in range(sigma) ]
def simba(n, sigma):

	vertexes = dict()      # list of the position of each node
	edges    = []          # edges of the isogeny triangle

	edge_colors   = []    # color of each edge: blue for scalar multiplications, and orange for isogeny evalutions
	vertex_colors = []    # color of each node: red for the leaves, otherwise color is set to white

	batches = BATCHES(n, sigma)
	t = -1
	acc = 0
	for i in range(sigma):

		t += 1
		vertexes[t] = (acc, 0)					# Initial vertex of the current batch
		vertex_colors.append('black')

		t += 1
		vertexes[t] = (acc, - n + batches[i])	# Root of the current batch
		vertex_colors.append('black')

		edges.append( (t - 1, t) )				# New edge to be added
		edge_colors.append('tab:blue')				# Color of this type of edge is always blue

		vertex = vertexes[t]
		for j in range(batches[i] - 1):

			t += 1
			vertexes[t] = (vertexes[t-1][0], vertexes[t-1][1] - batches[i] + j + 1)	 # Next leaf
			vertex_colors.append('black')

			edges.append( (t - 1, t) )				# New edge to be added
			edge_colors.append('tab:blue')				# Color of this type of edge is always blue

			t += 1
			vertexes[t] = (vertexes[t-2][0] + 1, vertexes[t-2][1])
			vertex_colors.append('black')

			edges.append( (t - 2, t) )				# New edge to be added
			edge_colors.append('tab:red')				# Color of this type of edge is always blue

		acc += batches[i]
	return vertexes, vertex_colors, edges, edge_colors



# ---------------------------------------------------------------------------------------
n = 74  # Number of primes: {l_1, l_2, ..., l_n}


# OAYT

S = [
	 43,  16,   7,   4,   3,   2,   1,   3,   2,   1,   6,   5,   4,   3,   2,   1,  10,   5,   4,   3,   2,   1,   9,   8,   7,
	  6,   5,   4,   3,   2,   1,  31,  11,  10,   9,   8,   7,   6,   5,   4,   3,   2,   1,  22,   8,   7,   6,   5,   4,   3,
	  2,   1,  21,  20,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,   9,   8,   7,   6,   5,   4,   3,   2,   1
    ]

'''
# MCR
S = [
43,  16,   7,   4,   2,   1,   1,   2,   1,   1,   4,   2,   1,   3,   2,   1,  10,   4,   1,   3,   2,   1,   7,   2,   1,
	  5,   1,   4,   3,   2,   1,  31,   9,   2,   1,   7,   1,   6,   5,   4,   3,   2,   1,  24,   6,   5,   4,   3,   2,   1,
	 19,   4,   3,   2,   1,  15,   3,   2,   1,  12,   2,   1,  10,   1,   9,   8,   7,   6,   5,   4,   3,   2,   1
    ]
'''
'''
# Dummy-free
S = [
	 42,  16,   7,   4,   4,   3,   2,   1,   3,   2,   1,   6,   5,   4,   3,   2,   1,  10,   5,   4,   3,   2,   1,   9,   8,
	  7,   6,   5,   4,   3,   2,   1,  30,  11,  10,   9,   8,   7,   6,   5,   4,   3,   2,   1,  21,   8,   7,   6,   5,   4,
	  3,   2,   1,  20,  19,  18,  17,  16,  15,  14,  13,  12,  11,  10,   9,   8,   7,   6,   5,   4,   3,   2,   1
    ]
'''

#S = list(range(n - 1, 0, -1))
# Multiplicative strategy is written as a graph
vertexes, vertex_colors, edges, edge_colors = strategy_evaluation(S, n)
#vertexes, vertex_colors, edges, edge_colors = simba(n, 3)

#vertexes, vertex_colors, edges, edge_colors = DRT(n)
G = nx.Graph()

# Adding nodes in specific positions
G.add_nodes_from(range(len(vertexes)))
nx.set_node_attributes(G, 'pos', vertexes)

# Adding edges with specific colors
for i in range(len(edges)):
    G.add_edge(edges[i][0], edges[i][1], color=edge_colors[i])

# Setting variables for a pretty plot of the graph
edges = G.edges()
edge_colors = [G[u][v]['color'] for u,v in edges]
weights = [6 for u,v in edges]
vertex_sizes = [24]*len(vertexes)

# Finally, the graph will be plotted
plt.figure(1,figsize=(17,17))
nx.draw(G, vertexes, node_color=vertex_colors, node_size=vertex_sizes, edge_color=edge_colors, width=weights, edge_labels=True)

# Saving the graph as a .PNG figure
plt.savefig('strategy.png')

plt.show()
