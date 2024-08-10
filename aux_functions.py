import gurobipy as gp
import networkx as nx
import math
import time
import os
import csv
from gurobipy import GRB

def get_state_codes():
	return  {
	'WA': '53', 'DE': '10', 'WI': '55', 'WV': '54', 'HI': '15', 'SD': '46',
	'FL': '12', 'WY': '56', 'NJ': '34', 'NM': '35', 'TX': '38', 'OR': '41',
	'LA': '22', 'NC': '37', 'ND': '38', 'NE': '31', 'TN': '47', 'NY': '36',
	'PA': '42', 'AK': '02', 'NV': '32', 'NH': '33', 'VA': '51', 'CO': '08',
	'CA': '06', 'AL': '01', 'AR': '05', 'VT': '50', 'IL': '17', 'GA': '13',
	'IN': '18', 'IA': '19', 'MA': '25', 'AZ': '04', 'ID': '16', 'CT': '09',
	'ME': '23', 'MD': '24', 'OK': '40', 'OH': '39', 'UT': '49', 'MO': '29',
	'MN': '27', 'MI': '26', 'RI': '44', 'KS': '20', 'MT': '30', 'MS': '28',
	'SC': '45', 'KY': '21'
}

def get_congressional_codes():
	return {
	'WA': 10, 'DE': 1, 'WI': 8, 'WV': 2, 'HI': 2,
	'FL': 28, 'WY': 1, 'NJ': 12, 'NM': 3, 'TX': 38,
	'LA': 6, 'NC': 14, 'ND': 1, 'NE': 3, 'TN': 9, 'NY': 26,
	'PA': 17, 'AK': 1, 'NV': 4, 'NH': 2, 'VA': 11, 'CO': 8,
	'CA': 52, 'AL': 7, 'AR': 4, 'VT': 1, 'IL': 17, 'GA': 14,
	'IN': 9, 'IA': 4, 'MA': 9, 'AZ': 9, 'ID': 2, 'CT': 5,
	'ME': 2, 'MD': 8, 'OK': 5, 'OH': 15, 'UT': 4, 'MO': 8,
	'MN': 8, 'MI': 13, 'RI': 2, 'KS': 4, 'MT': 2, 'MS': 4,
	'SC': 7, 'KY': 6, 'OR': 6, 'SD': 1
}

def read_gerrychain(state, parcel_level, s, num_minority_districts):
	with open('./results/gerrychain/' + state + '_' + parcel_level + '_(' + s + ',' + num_minority_districts + ').pckl', 'rb') as doc:
		gerrychain_result = pickle.load(doc)
		minority_districts = gerrychain_result[0]
		majority_districts = gerrychain_result[1]

		return minority_districts, majority_districts


def most_possible_nodes_in_one_district(population, U):
	cumulative_population = 0
	num_nodes = 0
	for ipopulation in sorted(population):
		cumulative_population += ipopulation
		num_nodes += 1
		if cumulative_population > U:
			return num_nodes - 1


def stabliity_number_callback(model, where):
	if where == gp.GRB.Callback.MIPNODE:
		obj = model.cbGet(gp.GRB.Callback.MIPNODE_OBJBST)

		if obj > model._k + 1:
			model.terminate()


def s_call_back(m, where):

	if where == gp.GRB.Callback.MIPSOL:

		x_sol = m.cbGetSolution(m._X)
		z_sol = m.cbGetSolution(m._Z)

		for (u, v) in m._complement_power_graph.edges():
			for j in range(m._upper_bound_minority_districts):
				if x_sol[u, j] + x_sol[v, j] > z_sol[j]:
					m.cbLazy(m._X[u, j] + m._X[v, j] <= m._Z[j])


def compuete_stability_number(method, graph, state, k, s):
	
	if method == 'read':
		if os.path.exists('./results/max_independent_set/' + state + '_' + str(s) + '.csv'):
			with open('./results/max_independent_set/' + state + '_' + str(s) + '.csv') as csvfile:
				csvreader = csv.reader(csvfile, delimiter=',')
			
				for row in csvreader:
					
					res = row[0].strip('][').split(', ')
					max_independent_set = [int(i) for i in res]
			
			return max_independent_set
			
	graph = nx.power(graph, s)
	start_time = time.time()

	m = gp.Model()
	m._k = k
	#m.Params.Cutoff = 10
	m._X = m.addVars(graph.nodes(), vtype=gp.GRB.BINARY)

	


	m.setObjective(m._X.sum(), gp.GRB.MAXIMIZE)
	######REVISIT WITH HAMID
	#m.addConstr(gp.quicksum(m._X) <= k + 1)
	m.addConstrs(m._X[i] + m._X[j] <= 1 for (i,j) in graph.edges())
	#for (i, j) in graph.edges():
	#	conflict_constr = m.addConstr(m._X[i] + m._X[j] <= 1)
	#	conflict_constr.Lazy = 3

	
	m.optimize(stabliity_number_callback)
	
	max_independent_set = []


	if m.status == gp.GRB.TIME_LIMIT:
		output_string = 'timed_out'
		
	else:
		output_string = m.getAttr(gp.GRB.Attr.ObjVal)
		for v in graph.nodes():
			if m._X[v].X == 1:
				max_independent_set.append(v)
		print(max_independent_set)

	end_time = time.time()

	with open('./results/max_independent_set/' + state + '_' + str(s) + '.csv', 'w') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		writer.writerow([max_independent_set, end_time - start_time])
		
	#csvwriter.write([str(max_independent_set), str(end_time - start_time)])

	return output_string, max_independent_set


def convert_interesting_s_vals_to_str(s_vals):
	ouput_str = ''
	for s_val in s_vals:
		ouput_str += ('(' + str(s_val[0]) + '|' + str(s_val[1]) + ') ')
	return ouput_str


def find_fischetti_separator(DG, component, b, n):
	if n == 'default':
		n = len(DG.nodes())

	neighbors_component = [False for i in range(n)]

	for i in nx.node_boundary(DG, component, None):
		#print(i)
		neighbors_component[i] = True

	visited = [False for i in range(n)]
	child = [b]
	visited[b] = True

	while child:
		parent = child
		child = []
		for i in parent:
			if not neighbors_component[i]:
				for j in DG.neighbors(i):
					if not visited[j]:
						child.append(j)
						visited[j] = True

	C = [ i for i in DG.nodes if neighbors_component[i] and visited[i] ]
	return C


def compactness_callback(m, where):

	if where == gp.GRB.Callback.MIPSOL:

		selector = m._option

		district_vars_exist = m._district_vars_exist
		G = m._G
		DG = nx.to_directed(G)
		s = m._s
		k = m._k
		population = m._population
		L = m._L
		n = m._n

		if selector != 'majority':
			xval = m.cbGetSolution(m._X)
			if district_vars_exist:
				zval = m.cbGetSolution(m._Z)
			minority_district_range = m._minority_district_range
		if selector != 'minority':
			yval = m.cbGetSolution(m._Y)
			if district_vars_exist:
				wval = m.cbGetSolution(m._W)
			majority_district_range = m._majority_district_range

		district_labels = []
		if selector == 'minority':
			for j in minority_district_range:
				district_labels.append(j)
		else:
			for j in range(k):
				district_labels.append(j)

		#if selector == 'both':
		#	district_labels = []
		#	for j in range(k):
		#		district_labels.append(j)
		#print(district_labels)
		##print(majority_district_range[-1])
		#print(k)
		###
		"""
		if selector == 'minority':
			list_of_nodes = []
			for vertex in G.nodes():
				for j in minority_district_range:
					if xval[vertex, j] > 0.5:
						list_of_nodes.append(vertex)
			list_of_majority_verteicies = list(set(G.nodes()) - set(list_of_nodes))
			majority_subgraph = G.subgraph(list_of_majority_verteicies)
			smallest = min(nx.connected_components(majority_subgraph), key = len)
			pop_of_smallest = 0
			for vertex in smallest:
				pop_of_smallest += m._population[vertex]
			if pop_of_smallest < L:
				m.cbLazy(gp.quicksum(1 - gp.quicksum(m._X[v, j] for j in minority_district_range) for v in smallest) <= len(smallest) - 1 + gp.quicksum(1 - gp.quicksum(m._X[bv, j] for j in minority_district_range) for bv in nx.node_boundary(G, smallest)))
		"""
		###

		for j in district_labels:
			# vertices assigned to this district (label j)
			if selector == 'minority':
				V_j = [v for v in G.nodes() if xval[v, j] > 0.5]

			elif selector == 'majority':
				V_j = [v for v in G.nodes() if yval[v, j] > 0.5]
			else:
				if j in minority_district_range:
					V_j1 = [v for v in G.nodes() if xval[v, j] > 0.5]
				else:
					V_j1 = []

				if j in majority_district_range:
					V_j2 = [v for v in G.nodes() if yval[v, j] > 0.5]
				else:
					V_j2 = []

				V_j = V_j1 + V_j2

			subgraph_j = DG.subgraph(V_j)
			if not subgraph_j:
				continue
			if not nx.is_strongly_connected(subgraph_j):
				smallest = min(nx.strongly_connected_components(subgraph_j), key = len)
				b = list(smallest)[0]
				for component in nx.strongly_connected_components(subgraph_j):

					if b in component: continue
					a = list(component)[0]
					C = find_fischetti_separator(DG, component, b, n)
					for (u,v) in DG.edges():
						DG[u][v]['sep-weight'] = 1
					# "remove" vertices from C
					for c in C:
						for neighbor in DG.neighbors(c):
							DG[c][neighbor]['sep-weight'] = s+1
					# check if C\{c} is a length-s a,b-separator
					remove_from_C=[]
					for c in C:
						for neighbor in DG.neighbors(c):
							DG[c][neighbor]['sep-weight'] = 1
						# calculate distance between a and b
						distance_from_a = nx.single_source_dijkstra_path_length(DG, a, weight='sep-weight')
						#distance_from_c = nx.single_source_dijkstra_path_length(DG, c, weight='sep-weight')
						if distance_from_a[b] > s:
							remove_from_C.append(c)
						else:
							for neighbor in DG.neighbors(c):
								DG[c][neighbor]['sep-weight'] = s+1
					# create minimal separator
					minC = [c for c in C if c not in remove_from_C]

					if selector != 'majority':
						for index in minority_district_range:
							if district_vars_exist:
								m.cbLazy(m._X[a,index] + m._X[b,index] <= m._Z[index] + gp.quicksum(m._X[c,index] for c in minC))
							else:
								m.cbLazy(m._X[a,index] + m._X[b,index] <= 1 + gp.quicksum(m._X[c,index] for c in minC))
					if selector != 'minority':
						#below may be wrong
						for index in majority_district_range:
						#above may be wrong
							if district_vars_exist:
								m.cbLazy(m._Y[a,index] + m._Y[b,index] <= m._W[index] + gp.quicksum(m._Y[c,index] for c in minC))
							else:
								m.cbLazy(m._Y[a,index] + m._Y[b,index] <= 1 + gp.quicksum(m._Y[c,index] for c in minC))

			else:
				G_j = G.subgraph(V_j)
				sub_diam = nx.diameter(G_j)
				if sub_diam <= s: continue
				# We first minimalize V \ V_j. Start with S' = V \ V_j.
				S_prime = list(set(G.nodes) - set(V_j))
				for a in V_j:
					distance_from_a = nx.single_source_dijkstra_path_length(G_j, a)
					for vertex in distance_from_a:
						b = vertex
						if a < b and distance_from_a[b] > s:

							distance_from_a_in_G = nx.single_source_dijkstra_path_length(G, a)
							distance_from_b_in_G = nx.single_source_dijkstra_path_length(G, b)
							remove_from_S_prime = []
							for v in S_prime:
								if distance_from_a_in_G[v] + distance_from_b_in_G[v] > s:
									remove_from_S_prime.append(v)
							minimal_S_prime = [node for node in S_prime if node not in remove_from_S_prime]

							# Now we define C as a length-s a,b separator and we minimalize it
							C = minimal_S_prime
							remove_from_C = []

							for u, v in G.edges():
								G[u][v]['sep-weight'] = 1

							for c in C:
								for neighbor in G.neighbors(c):
									G[c][neighbor]['sep-weight'] = s + 1

							for c in C:
								for neighbor in G.neighbors(c):
									G[c][neighbor]['sep-weight'] = 1

								#G_prime = G.subgraph(list(set(V_j).union({c})))
								#print ("Is c in G_prime?", c in G_prime.nodes)
								###distance_from_a_in_G_prime = nx.single_source_dijkstra_path_length(G_prime, a)
								###distance_from_b_in_G_prime = nx.single_source_dijkstra_path_length(G_prime, b)
								distance_from_a = nx.single_source_dijkstra_path_length(G, a, weight='sep-weight')

								if distance_from_a[b] > s:
									remove_from_C.append(c)
								else:
									for neighbor in G.neighbors(c):
										G[c][neighbor]['sep-weight'] = s + 1

							minC = [c for c in C if c not in remove_from_C]
							#print("Lazy cut 2 is added")
							if selector != 'majority':
								for index in minority_district_range:
									if district_vars_exist:
										m.cbLazy(m._X[a,index] + m._X[b,index] <= m._Z[index] + gp.quicksum(m._X[c,index] for c in minC))
									else:
										m.cbLazy(m._X[a,index] + m._X[b,index] <= 1 + gp.quicksum(m._X[c,index] for c in minC))
							if selector != 'minority':
								for index in majority_district_range:
									if district_vars_exist:
										m.cbLazy(m._Y[a,index] + m._Y[b,index] <= m._W[index] + gp.quicksum(m._Y[c,index] for c in minC))
									else:
										m.cbLazy(m._Y[a,index] + m._Y[b,index] <= 1 + gp.quicksum(m._Y[c,index] for c in minC))


def continuous_ONLY_callback(m, where):

	if where == gp.GRB.Callback.MIPSOL:

		selector = m._option

		district_vars_exist = m._district_vars_exist
		G = m._G
		DG = nx.to_directed(G)
		k = m._k
		n = len(G.nodes())

		if selector != 'majority':
			xval = m.cbGetSolution(m._X)
			if district_vars_exist:
				zval = m.cbGetSolution(m._Z)
			minority_district_range = m._minority_district_range
		if selector != 'minority':
			yval = m.cbGetSolution(m._Y)
			if district_vars_exist:
				wval = m.cbGetSolution(m._W)
			majority_district_range = m._majority_district_range

		district_labels = []
		if selector != 'majority':
			for j in minority_district_range:
				district_labels.append(j)
		if selector != 'minority':
			for j in majority_district_range:
				district_labels.append(j)

		if selector == 'both':
			district_labels = []
			for j in majority_district_range:
				district_labels.append(j)

		for j in district_labels:
			# vertices assigned to this district (label j)
			if selector == 'minority':
				V_j = [v for v in G.nodes() if xval[v, j] > 0.5]

			elif selector == 'majority':
				V_j = [v for v in G.nodes() if yval[v, j] > 0.5]
			else:
				if j in minority_district_range:
					V_j = [v for v in G.nodes() if xval[v, j] > 0.5]
				else:
					V_j = [v for v in G.nodes() if yval[v, j] > 0.5]

			subgraph_j = DG.subgraph(V_j)
			if not subgraph_j:
				continue
			if not nx.is_strongly_connected(subgraph_j):
				smallest = min(nx.strongly_connected_components(subgraph_j), key = len)
				b = list(smallest)[0]
				for component in nx.strongly_connected_components(subgraph_j):

					if b in component: continue
					a = list(component)[0]
					minC = find_fischetti_separator(DG, component, b, n)


					if selector != 'majority':
						for index in minority_district_range:
							if district_vars_exist:
								m.cbLazy(m._X[a,index] + m._X[b,index] <= m._Z[index] + gp.quicksum(m._X[c,index] for c in minC))
							else:
								m.cbLazy(m._X[a,index] + m._X[b,index] <= 1 + gp.quicksum(m._X[c,index] for c in minC))
					if selector != 'minority':
						#below may be wrong
						for index in majority_district_range:
						#above may be wrong
							if district_vars_exist:
								m.cbLazy(m._Y[a,index] + m._Y[b,index] <= m._W[index] + gp.quicksum(m._Y[c,index] for c in minC))
							else:
								m.cbLazy(m._Y[a,index] + m._Y[b,index] <= 1 + gp.quicksum(m._Y[c,index] for c in minC))


def break_point_plot(title, break_points, filename):
	x = []
	y = []
	for point in break_points:
		x.append(point[0])
		y.append(point[1])

	plt.figure()
	plt.title(title)
	plt.ylabel('Majority-minority district number')
	plt.xlabel('s value')
	plt.scatter(x, y)
	plt.savefig(filename + title + 's_val_plot.png')


def export_to_png(G, threshold, df, districts, filename1, filename2):

	assignment = [ -1 for u in G.nodes ]

	minority_district = [0 for u in G.nodes]

	for j in range(len(districts)):
		total_population = 0
		minority_population = 0
		for i in districts[j]:
			geoID = G.nodes[i]['GEOID20']
			total_population += G.node[i]['VAP']
			minority_population += G.node[i]['BVAP']
			#minority_population += G.node[i]['HVAP']
			for u in G.nodes:
				if geoID == df['GEOID20'][u]:
					assignment[u] = j



		if minority_population > threshold*total_population:
			print('Threshold is ', minority_population/total_population, ' for district ', j)
			for v in districts[j]:
				geoID = G.nodes[v]['GEOID20']
				for u in G.nodes:
					if geoID == df['GEOID20'][u]:
						minority_district[u] = 1
		else:
			print('Threshold is ', minority_population/total_population, ' for district ', j)


	if min(assignment[v] for v in G.nodes) < 0:
		print('Error: did not assign all nodes in district map png.')
	else:
		df['assignment'] = assignment
		my_fig = df.plot(column='assignment').get_figure()
		#RESIZE_FACTOR = 3
		#my_fig.set_size_inches(my_fig.get_size_inches()*RESIZE_FACTOR)
		plt.axis('off')
		my_fig.savefig(filename1)

		df['minority'] = minority_district
		cmap = LinearSegmentedColormap.from_list('minority', [(0, 'white'), (1, 'lightgray')])
		splot = df.plot(cmap=cmap, column='minority',figsize=(10, 10), linewidth=1, edgecolor='0.25').get_figure()  # display the S map
		plt.axis('off')
		splot.savefig(filename2)


def upper_bound_s_callback(m, where):

	if where == gp.GRB.Callback.MIPSOL:
		DG = m._DG
		v = m._v
		tval = m.cbGetSolution(m._t)
		gval = m.cbGetSolution(m._g)
		selected_arcs = [key for key in gval if gval[key] > 0.5]
		selected_verticies = [key for key in tval if tval[key]  > 0.5]
		solution_subgraph = DG.edge_subgraph(selected_arcs)

		for connected_component in nx.weakly_connected_components(solution_subgraph):
			if v not in connected_component:
				

				component_complement = [i for i in DG.nodes if i not in connected_component]
				boundary_edges = list(nx.edge_boundary( DG, component_complement, connected_component))
				#m.cbLazy( gp.quicksum(m._g[i,j] for i,j in boundary_edges) >= gp.quicksum(m._t[u] for u in selected_verticies))
				for u in connected_component:
					m.cbLazy(gp.quicksum(m._g[i,j] for i,j in boundary_edges) >= m._t[u])
   

def DFJ_callback(m, where):
    
    # check if LP relaxation at this BB node is integer
    if where == gp.GRB.Callback.MIPSOL: 
        
        # retrieve the LP relaxation solution at this BB node
        qval = m.cbGetSolution(m._q)
        
        # which edges are selected in the LP solution?
        chosen_edges = [(i,j) for i,j in m._DG.edges if qval[i,j] > 0.5]
        
        # if the solution is not a tree, it will have multiple pieces
        for component in nx.weakly_connected_components(m._DG.edge_subgraph(chosen_edges)):
            
            # each piece that does not contain the tree root r is a cycle (w/ too many edges)
            if m._r_q not in component:
                
                # must pick fewer than |component| interior edges
                interior_edges = nx.edge_boundary(m._DG, component, component)
                
                m.cbLazy(gp.quicksum(m._q[i,j] for i,j in interior_edges) <= len(component) - 1)


def recursive_bipartition_heuristic(G, gingles_districts=list()):
    
    assigned_nodes = [ i for j in range(len(gingles_districts)) for i in gingles_districts[j] ]
    unassigned_nodes = [ i for i in G.nodes if i not in assigned_nodes ]
    
    clusters = [ unassigned_nodes ]
    sizes = [ G._k - len(gingles_districts) ]
    districts = gingles_districts.copy()
    num_multidistricts = 0
    DG = nx.DiGraph(G)

    while len(sizes) > 0:

        # pick a cluster
        cluster = clusters.pop()
        size = sizes.pop()

        size1 = math.floor( size / 2 )
        size2 = math.ceil( size / 2 )

        DH = DG.subgraph(cluster)
        DH._L = size1 * G._L
        DH._U = size1 * G._U
        DH._CL = size2 * G._L
        DH._CU = size2 * G._U

        print("Using one split county, attempting to bipartition cluster into sizes:",size1,size2)

        # first, try to bipartition using 1 split county
        m = build_single_district_mip(DH, objective='cut_edges', split_counties_limit=1, deviation_penalty=0.001, verbose=False,
                                      contiguity='cut', complement_contiguity='cut', complement_balance=True)

        m.Params.TimeLimit = 3600
        m.optimize(m._callback)

        if m.solCount == 0:

            # otherwise, allow any number of split counties
            print("Without limiting splits, attempting to bipartition cluster")
            m = build_single_district_mip(DH, objective='cut_edges', deviation_penalty=0.001, verbose=False,
                                      contiguity='cut', complement_contiguity='cut', complement_balance=True)

            m.Params.TimeLimit = 3600
            m.optimize(m._callback)

            if m.solCount == 0:
                print("Unable to bipartition. Keeping as multidistrict of size =",size)
                districts.append(cluster)
                num_multidistricts += 1
                continue

        cluster1 = [ i for i in DH.nodes if m._x[i].x > 0.5 ]
        cluster2 = [ i for i in DH.nodes if m._x[i].x < 0.5 ]

        if size1 == 1:
            districts.append(cluster1)
        else:
            clusters.append(cluster1)
            sizes.append(size1)

        if size2 == 1:
            districts.append(cluster2)
        else:
            clusters.append(cluster2)
            sizes.append(size2)
            
    print(f"After recursive partitioning, we have {len(districts)} districts.")
    if num_multidistricts > 0:
        print(f"This includes {num_multidistricts} multidistricts.")
    return districts

def build_single_district_mip(DG, objective='polsby_popper', contiguity=None, root=None, verbose=True,
                              split_counties_limit=None, deviation_penalty=0.0,
                              minority=None, mvap_lower=0.5, mvap_upper=1.0, mvap_excess_penalty=0.0, 
                              complement_contiguity=None, complement_balance=False, complement_compactness_penalty=0.0):
    
    # sanity check the inputs
    assert objective in {'polsby_popper', 'cut_edges'}
    assert minority in { None, 'Asian', 'Black', 'Hispanic', 'Native' }
    assert contiguity in { None, 'shir', 'cut' }
    assert complement_contiguity in { None, 'cut' } # shir not supported
    assert objective=='polsby_popper' or complement_compactness_penalty==0.0, "Can use complement_compactness_penalty only for PP objective"
    if contiguity == 'shir':
        assert root is not None
        assert root in DG.nodes
        
    ##################################
    # CREATE MODEL AND MAIN VARIABLES
    ##################################
    
    m = gp.Model()
    if not verbose:
        m.Params.OutputFlag = 0
    
    # x[i] equals one when node i is selected in the district
    m._x = m.addVars(DG.nodes, name='x', vtype=GRB.BINARY)

    # y[u,v] equals one when arc (u,v) is cut because u (but not v) is selected in the district
    m._y = m.addVars(DG.edges, name='y', vtype=GRB.BINARY)

    ###########################
    # ADD MAIN CONSTRAINTS
    ###########################
    
    # add constraints saying that the district has population at least L and at most U
    m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * m._x[i] for i in DG.nodes) >= DG._L )
    m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * m._x[i] for i in DG.nodes) <= DG._U )
    
    # add constraints saying that the complement has population at least CL and at most CU
    if complement_balance:
        m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * (1-m._x[i]) for i in DG.nodes) >= DG._CL )
        m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * (1-m._x[i]) for i in DG.nodes) <= DG._CU )

    # add constraints saying that edge {u,v} is cut if u (but not v) is selected in the district
    m.addConstrs( m._x[u] - m._x[v] <= m._y[u,v] for u,v in DG.edges )
    
    ###########################
    # ADD OBJECTIVE 
    ###########################
    
    if objective == 'polsby_popper':
        # z is inverse Polsby-Popper score for the district
        m._z = m.addVar(name='z')

        # objective is to minimize the inverse Polsby-Popper score
        m.setObjective( m._z, GRB.MINIMIZE )

        if complement_compactness_penalty > 0:
            m._zc = m.addVar(name='zc')
            m._zc.obj = complement_compactness_penalty
        
    elif objective == 'cut_edges':
        undirected_edges = [ (i,j) for i,j in DG.edges if i < j ]
        m.addConstrs( m._y[i,j] == m._y[j,i] for i,j in undirected_edges )
        m.setObjective( gp.quicksum( m._y[i,j] for i,j in undirected_edges ), GRB.MINIMIZE )
    
    ###################################
    # ADD POLSBY-POPPER CONSTRAINTS 
    ###################################
    
    if objective == 'polsby_popper':
        
        # A = area of the district
        m._A = m.addVar(name='A')

        # P = perimeter of the district
        m._P = m.addVar(name='P')

        # add SOCP constraint relating inverse Polsby-Popper score z to area and perimeter
        m.addConstr( m._P * m._P <= 4 * math.pi * m._A * m._z )

        # add constraint on area A
        m.addConstr( m._A == gp.quicksum( DG.nodes[i]['area'] * m._x[i] for i in DG.nodes ) )

        # add constraint on perimeter P
        m.addConstr( m._P == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v] for u,v in DG.edges )
                     + gp.quicksum( DG.nodes[i]['boundary_perim'] * m._x[i] for i in DG.nodes if DG.nodes[i]['boundary_node'] ) )

        if complement_compactness_penalty > 0:
            m._Ac = m.addVar(name='Ac')
            m._Pc = m.addVar(name='Pc')
            m.addConstr( m._Pc * m._Pc <= 4 * math.pi * m._Ac * m._zc )
            m.addConstr( m._Ac == gp.quicksum( DG.nodes[i]['area'] * (1-m._x[i]) for i in DG.nodes ) )
            m.addConstr( m._Pc == gp.quicksum( DG.edges[u,v]['shared_perim'] * m._y[u,v] for u,v in DG.edges )
                         + gp.quicksum( DG.nodes[i]['boundary_perim'] * (1-m._x[i]) for i in DG.nodes if DG.nodes[i]['boundary_node'] ) )
        
        m.update()
    
    ###################################
    # ADD MINORITY CONSTRAINT(S)
    ###################################
    
    if minority is not None:
        
        codes = get_census_codes(minority)
        
        for i in DG.nodes:
        
            # voting age population (VAP)
            DG.nodes[i]['VAP'] = DG.nodes[i]['P0030001']
        
            # minority voting age population (MVAP)
            DG.nodes[i]['MVAP'] = sum( DG.nodes[i][code] for code in codes )
            
        # Idea: impose mvap >= 0.5 * vap
        m._mvap = m.addVar(name='mvap')
        m._vap = m.addVar(name='vap')

        m.addConstr( m._mvap == gp.quicksum( DG.nodes[i]['MVAP'] * m._x[i] for i in DG.nodes ) )
        m.addConstr( m._vap == gp.quicksum( DG.nodes[i]['VAP'] * m._x[i] for i in DG.nodes ) )

        excess = m.addVar(name='excess')
        excess.obj = mvap_excess_penalty
        m.addConstr( m._mvap - excess == mvap_lower * m._vap )
        
        # don't exceed, say, 80% BVAP
        if mvap_upper < 1:
            m.addConstr( m._mvap <= mvap_upper * m._vap )
        
    ###################################
    # ADD CONTIGUITY CONSTRAINTS
    ###################################
    
    m._callback = None
    m._numCallbacks = 0
    m._numLazyCuts = 0
    
    if contiguity == 'shir':
        
        m._x[root].LB = 1
        M = DG.number_of_nodes() - 1
        
        # Add flow variables: f[u,v] = amount of flow sent across arc uv 
        m._f = m.addVars( DG.edges, name='f' )
        
        # if selected but not a root, consume one unit of flow
        m.addConstrs( gp.quicksum( m._f[j,i] - m._f[i,j] for j in DG.neighbors(i) ) == m._x[i] for i in DG.nodes if i != root )

        # flow can only enter selected nodes
        m.addConstrs( gp.quicksum( m._f[j,i] for j in DG.neighbors(i) ) <= M * m._x[i] for i in DG.nodes if i != root )
        
    if contiguity == 'cut' or complement_contiguity == 'cut':
        
        m.Params.LazyConstraints = 1
        m._DG = DG
        m._root = root
        m._contiguity = contiguity
        m._complement_contiguity = complement_contiguity
        m._callback = short_burst_cut_callback
        
    ##########################
    # LIMIT SPLIT COUNTIES?
    ##########################
    
    if split_counties_limit is not None:
        
        fips = list( { DG.nodes[i]['GEOID20'][0:5] for i in DG.nodes } )
        is_all = m.addVars(fips, vtype=GRB.BINARY) # is all of county selected?
        is_some = m.addVars(fips, vtype=GRB.BINARY) # is (at least) some of county selected?
        is_split = m.addVars(fips, vtype=GRB.BINARY) # is county split?

        m.addConstr( gp.quicksum( is_split ) <= split_counties_limit )
        m.addConstrs( is_split[c] == is_some[c] - is_all[c] for c in fips )

        for i in DG.nodes:
            c = DG.nodes[i]['GEOID20'][0:5]
            m.addConstr( is_all[c] <= m._x[i] )
            m.addConstr( m._x[i] <= is_some[c] )
            
    ###################################
    # IDEAL PENALTY
    ###################################
    
    if deviation_penalty > 0.0:
        
        below = m.addVar()
        below.obj = deviation_penalty
        
        above = m.addVar()
        above.obj = deviation_penalty
        
        total = sum( DG.nodes[i]['TOTPOP'] for i in DG.nodes )
        ideal = total * ( DG._L + DG._U ) / ( DG._L + DG._U + DG._CL + DG._CU )
        m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * m._x[i] for i in DG.nodes) == ideal + above - below )
        
    ###################################
    # SOLVE PARAMETERS
    ###################################
    
    m.Params.MIPGap = 0.00
    m.Params.FeasibilityTol = 1e-7
    m.Params.IntFeasTol = 1e-7
    m.update()
    
    return m

def short_burst_cut_callback(m, where):
    
    if where == GRB.Callback.MIPSOL:
        m._numCallbacks += 1 
        DG = m._DG
        xval = m.cbGetSolution(m._x)

        ##########################################
        # ADD CUT FOR COMPLEMENT?
        ##########################################
        
        if m._contiguity == 'cut':
            # vertices assigned to this district 
            S = [ v for v in DG.nodes if xval[v] > 0.5 ]

            # what shall we deem as the "root" of this district? call it b
            b = m._root # possibly None
            
            # for each component that doesn't contain b, add a cut
            for component in sorted( nx.strongly_connected_components( DG.subgraph(S) ), key=len, reverse=True ):

                # what is the maximum population node in this component?
                maxp = max( DG.nodes[v]['TOTPOP'] for v in component)
                mpv = [ v for v in component if DG.nodes[v]['TOTPOP'] == maxp ][0]

                # if no root 'b' has been selected yet, pick one
                if b is None:
                    # find some vertex "b" that has largest population in this component
                    b = mpv

                if b in component: 
                    continue

                # find some vertex "a" that has largest population in this component
                a = mpv

                # get minimal a,b-separator
                C = short_burst_find_fischetti_separator(DG, component, b)

                # add lazy cut
                m.cbLazy( m._x[a] + m._x[b] <= 1 + gp.quicksum( m._x[c] for c in C ) )
                m._numLazyCuts += 1
            
        ##########################################
        # ADD CUT FOR COMPLEMENT?
        ##########################################
        
        if m._complement_contiguity == 'cut':
            
            # vertices assigned to the complement
            S = [ v for v in DG.nodes if xval[v] < 0.5 ]

            # what shall we deem as the "root" of the complement? call it b
            b = None
            
            # for each component that doesn't contain b, add a cut
            for component in sorted( nx.strongly_connected_components( DG.subgraph(S) ), key=len, reverse=True ):

                # what is the maximum population node in this component?
                maxp = max( DG.nodes[v]['TOTPOP'] for v in component )
                mpv = [ v for v in component if DG.nodes[v]['TOTPOP'] == maxp ][0]

                # if no root 'b' has been selected yet, pick one
                if b is None:
                    # find some vertex "b" that has largest population in this component
                    b = mpv
                    continue

                # find some vertex "a" that has largest population in this component
                a = mpv

                # get minimal a,b-separator
                C = short_burst_find_fischetti_separator(DG, component, b)

                # add lazy cut
                # replace x by 1-x in:
                #      m.cbLazy( m._x[a] + m._x[b] <= 1 + gp.quicksum( m._x[c] for c in C ) )
                # i.e., if neither a nor b is picked in district, then not all of C can be picked in district
                m.cbLazy( gp.quicksum( m._x[c] for c in C ) + 1 <= m._x[a] + m._x[b] + len(C) )
                m._numLazyCuts += 1
                
    return

def short_burst_find_fischetti_separator(DG, component, b):
    neighbors_component = { i : False for i in DG.nodes }
    for i in nx.node_boundary(DG, component, None):
        neighbors_component[i] = True
    
    visited = { i : False for i in DG.nodes }
    child = [ b ]
    visited[b] = True
    
    while child:
        parent = child
        child = []
        for i in parent:
            if not neighbors_component[i]:
                for j in DG.neighbors(i):
                    if not visited[j]:
                        child.append(j)
                        visited[j] = True
    
    C = [ i for i in DG.nodes if neighbors_component[i] and visited[i] ]
    return C

def fips_support(G, districts):
    fips = { G.nodes[i]['GEOID20'][0:5] for i in G.nodes }
    fs = { f : list() for f in fips }
    for j in range(len(districts)):
        for i in districts[j]:
            f = G.nodes[i]['GEOID20'][0:5]
            if j not in fs[f]:
                fs[f].append(j)
    return fs

def number_of_counties_split(G, districts, verbose=False):
    fs = fips_support(G, districts)
    if verbose:
        print("\nCounties split (by fips code):", [ f for f in fs.keys() if len(fs[f]) > 1 ])
    return sum( 1 for f in fs.keys() if len(fs[f]) > 1 )
    
def number_of_county_splits(G, districts, verbose=False):
    fs = fips_support(G, districts)
    if verbose:
        print("County splits (by fips code):", { f : len(fs[f])-1  for f in fs.keys() if len(fs[f]) > 1 })
    return sum( ( len(fs[f]) - 1 ) for f in fs.keys() )
    
def polsby_popper(G, district, label):
    area = sum( G.nodes[i]['area'] for i in district )
    perim = sum( G.edges[u,v]['shared_perim'] for u in district for v in G.neighbors(u) if label[u]!=label[v] )
    perim += sum( G.nodes[i]['boundary_perim'] for i in district if G.nodes[i]['boundary_node'] ) 
    return 4 * math.pi * area / ( perim * perim )

def average_polsby_popper(G, districts, verbose=False):
    label = { i : j for j in range(len(districts)) for i in districts[j] }
    if verbose:
        print("\nDistrict Polsby-Popper scores:")
        for p in range(len(districts)):
            print(p, round(polsby_popper(G, districts[p], label),4) )
    return sum( polsby_popper(G, district, label) for district in districts ) / len(districts) 

def number_of_gingles_districts(G, districts, minority, verbose=False):
    
    gingles_count = 0
    if verbose:
        print(f"\nDistrict {minority} percentages:")
    for p in range(len(districts)):
        district = districts[p]
        vap = sum( G.nodes[i]['VAP'] for i in district )
        mvap = sum( G.nodes[i]['MVAP'] for i in district )
        if mvap >= 0.5 * vap:
            gingles_count += 1
        if verbose:
            flag = '***' if mvap >= 0.5 * vap else ''
            print(p, round(100*mvap/vap,2), flag)
    return gingles_count

def report_metrics(G, districts, minority=None, verbose=False):
    if minority is not None:
        gingles = number_of_gingles_districts(G, districts, minority, verbose=verbose)
        #gingles = 100
        print(f"-> {gingles} majority-{minority} districts")
    print(f"-> {number_of_counties_split(G, districts, verbose=verbose)} counties split a total of {number_of_county_splits(G, districts, verbose=verbose)} times")
    avepp = round( average_polsby_popper(G, districts, verbose=verbose), 4)
    print("-> average Polsby-Popper score of",avepp)
    return