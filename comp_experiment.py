import sys
import gurobipy as gp
import networkx as nx
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib
import math
import pandas
import json
import os
import time
import subprocess
import random
import csv
import pickle
import numpy as np

from gurobipy import GRB
from gerrychain import MarkovChain
from gerrychain import (GeographicPartition, Graph, updaters, constraints, accept, proposals)

from gerrychain import Partition
from gerrychain import constraints as constraints_class
from gerrychain.optimization import Gingleator
from gerrychain.tree import recursive_tree_part
from gerrychain.proposals import recom
from functools import partial
from matplotlib.colors import LinearSegmentedColormap
from aux_functions import *
pandas.options.mode.chained_assignment = None


class problem_instance(object):

	def __init__(self, state, parcel_level, s):

		self.state = state
		self.parcel_level = parcel_level
		self.model_time_limit = 3600

		self.error_message = ''
		self.f = 0.5
		
		self.population_deviation = .01

		if parcel_level == 'county':
			parcel_level_plural = 'counties'

		if parcel_level == 'tract':
			parcel_level_plural = 'tracts'

		self.k = congressional_districts[state.upper()]

		print('Attempting to read in json data')
		if os.path.exists('./raw_data/' + parcel_level + '/json/' + state + '_' + parcel_level_plural + '.json'):
			self.G = Graph.from_json('./raw_data/' + parcel_level + '/json/' + state + '_' + parcel_level_plural + '.json')

			f = open('./raw_data/' + parcel_level + '/json/' + state + '_' + parcel_level_plural + '.json')
			self.parcel_data = json.load(f)

			self.population = [self.parcel_data['nodes'][i]['P0010001'] for i in self.G.nodes()]
			for v in self.G.nodes():
				self.G.nodes[v]['TOTPOP'] = self.population[v]
			self.voting_age_population = [self.parcel_data['nodes'][i]['P0030001'] for i in self.G.nodes()]
			for v in self.G.nodes():
				self.G.nodes[v]['VAP'] = self.voting_age_population[v]
			self.minority_population = [self.parcel_data['nodes'][i]['P0030004'] for i in self.G.nodes()]
			for v in self.G.nodes():
				self.G.nodes[v]['MVAP'] = self.minority_population[v]

			self.L = math.ceil((1 - self.population_deviation / 2) * sum(self.population) / self.k)
			self.M = sum(self.population) / self.k
			self.U = math.floor((1 + self.population_deviation / 2) * sum(self.population) / self.k)

			for i in self.G.nodes():
				if self.population[i] > math.floor((1 + self.population_deviation / 2) * sum(self.population) / self.k):
					self.error_message += ' Parcel too large'
			f.close()
			if not nx.is_connected(self.G):
				self.error_message += ' NOT CONNECTED'
				print(self.error_message)
		else:
			# no JSON data
			self.error_message += ' No json data'
			print(self.error_message)

		print('Attempting to read in shape data')
		if os.path.exists('./raw_data/' + parcel_level + '/shape/' + state + '_' + parcel_level_plural + '.shp'):
			self.shape_info = gpd.read_file('./raw_data/' + parcel_level + '/shape/' + state + '_' + parcel_level_plural + '.shp')
			self.plot = False
		else:
			self.error_message += ' no SHAPE data'
			print(self.error_message)

		if self.error_message == '':

			
			if os.path.exists('results/upper_bound_minority_districts/upper_bound_minority_districts_' + self.state + '_' + self.parcel_level + '.txt'):
				with open('results/upper_bound_minority_districts/upper_bound_minority_districts_' + self.state + '_' + self.parcel_level + '.txt', 'r') as file:
					for line in file:
						self.k_minority = int(line)
						print('Upper bound on minority districts read in as :', str(self.k_minority))

			else:
				with open('results/upper_bound_minority_districts/upper_bound_minority_districts_' + self.state + '_' + self.parcel_level + '.txt', 'w') as doc:
					self.minority_district_range = range(self.k)
					self.k_minority = self.k
					self.find_upper_bound_minority_districts()
					doc.write(str(self.k_minority))
					print('Upper bound on minority districts comupted and stored as :', str(self.k_minority))

			
			self.minority_district_range = range(self.k_minority)
			self.majority_district_range = range(self.k_minority, self.k)


			self.find_upper_bound_s()

			if os.path.exists('results/lower_bound_s/lower_bound_s_' + self.state + '_' + self.parcel_level + '.txt'):
				with open('results/lower_bound_s/lower_bound_s_' + self.state + '_' + self.parcel_level + '.txt', 'r') as file:
					for line in file:
						self.lower_bound_s = int(line)

				if not s:
					self.s = self.lower_bound_s
				else:
					self.s = s
				
				print('Read in lower_bound_s: ', str(self.lower_bound_s))
				print('Setting s to: ', str(self.s))

			else:
				print('No lower bound on s to read calculating now')
				self.s_method = 'independent_set'
				self.lower_bound_s = 1
				self.find_lower_bound_s()
				print('success, found lower_bound_s: ', str(self.lower_bound_s))
				
	def get_current_plan_stat(self):
		
		district_plan = gpd.read_file('./raw_data/current_plan/' + self.state + '.shp')
		parcels = self.shape_info

		plt.style.use('classic')
		parcels = parcels.set_crs(epsg=4269)
		
		
		district_plan = district_plan.to_crs(parcels.crs)
		plot = False
		if plot:
			f, ax = plt.subplots()
			parcels.plot(ax=ax, cmap = 'Paired', edgecolor='white', figsize = (10,12)).get_figure()
			district_plan.plot(ax=ax, categorical=True, cmap = 'Paired', edgecolor='black', facecolor='none', linewidth=3, figsize = (10,12)).get_figure()
			plt.show()

		districts = []

		for district_num, district in district_plan.iterrows():
			current_district = []
			district_polygon = district['geometry']
			
			for parcel_num, parcel in parcels.iterrows():
				parcel_polygon = parcel['geometry']

				
				
				if parcel_polygon.intersection(district_polygon).area > 0.0000001:
					current_district.append(parcel_num)

			districts.append(current_district)

		split_number = 0
		for parcel in self.G.nodes():
			counter = 0
			for district in districts:
				if parcel in district:
					counter += 1

			if counter > 1:
				
				split_number += 1

		current_s = 0
		for district in districts:
			subgraph = self.G.subgraph(district)
			if nx.diameter(subgraph) > current_s:
				current_s = nx.diameter(subgraph)
		
		current_number_minority_districts = 0

		for district in districts:
			minority_population = 0
			total_population = 0
			for parcel in district:
				minority_population += self.minority_population[parcel]
				total_population += self.population[parcel]		
		
		return current_s, split_number
	
	def compute_lower_bound_s(self):

		start_time = time.time()
		self.num_s_iter = 0
		s = 0


		upper = int(self.upper_bound_s)
		stability_number = self.k + 1

		while self.lower_bound_s != upper:

			self.num_s_iter += 1


			s = (upper + self.lower_bound_s) // 2
			print('upper: ', upper)
			print('lower: ', self.lower_bound_s)
			print('s: ', s)
		
			stability_number, max_independent_set = compuete_stability_number('compuete', self.G, self.state, self.k, s)
			

			if stability_number == 'timed_out':
				self.s_time_out = True
				break

			if stability_number > self.k:
				self.lower_bound_s = s + 1
			else:
				upper = s

		print('The smallest s is :' + str(self.lower_bound_s))

		end_time = time.time()

		self.time_for_s_lower_bound = end_time - start_time

	def find_lower_bound_s(self):

		self.compute_lower_bound_s()
		with open('results/lower_bound_s/lower_bound_s_' + self.state + '_' + self.parcel_level + '.txt', 'w') as file:
			file.write(str(self.lower_bound_s))

	def find_upper_bound_s(self):

		if os.path.exists('results/upper_bound_s/upper_bound_s_' + self.state + '_' + self.parcel_level + '.txt'):
			with open('results/upper_bound_s/upper_bound_s_' + self.state + '_' + self.parcel_level + '.txt', 'r') as file:
				for line in file:
					self.upper_bound_s = int(float(line))
					print("Upper bound on s read as: ", self.upper_bound_s)
		else:
			self.upper_bound_s = nx.diameter(self.G)

			with open('results/upper_bound_s/upper_bound_s_' + self.state + '_' + self.parcel_level + '.txt', 'w') as file:
				file.write(str(self.upper_bound_s))

			print("Upper bound on s computed as: ", self.upper_bound_s)

	def find_upper_bound_minority_districts(self, continuity=False, relax=False):

		self.minority_district_ONLY_model(0, relax, continuity)
		self.k_minority = self.num_minority_districts

	def find_interesting_s_vals(self):

		output_string = 's_vaule, num_minority_districts, time_out \n'
		s = self.lower_bound_s - 1

		continue_search = True

		if not os.path.exists('./interesting_s_vals.csv'):
			with open('./interesting_s_vals.csv', 'a') as doc:
				doc.write('state, s, number of districts, time out \n')

		while continue_search:
			s += 1

			self.minority_district_ONLY_model(s)
			
			if self.time_out == True:
				continue_search = False

			if type(self.num_minority_districts) != int:
				continue_search = False
			
			if s == self.upper_bound_s:
				continue_search = False

			if self.k_minority == self.num_minority_districts:
				continue_search = False
		
			current_string = str(s) + ', ' + str(self.num_minority_districts) + ',' + str(self.time_out) + '\n'
			output_string += current_string
			#self.break_points.append([s, self.num_minority_districts])

			with open('./interesting_s_vals.csv', 'a') as doc:
				current_string = self.state + ',' + current_string
				doc.write(current_string)
		with open('./results/interesting_s_vals/interesting_s_vals_' + self.state + '_' + self.parcel_level + '.csv', 'w') as doc:
			doc.write(output_string)

	def run_GerryChain_heuristic(self, iterations):

		s = -1
		my_updaters = {'population': updaters.Tally('P0030001', alias='population')}
		start = recursive_tree_part(self.G, range(self.k),sum(self.G.nodes[i]['P0030001'] for i in self.G.nodes())/self.k,'P0030001', self.population_deviation/2,1)
		initial_partition = GeographicPartition(self.G, start, updaters = my_updaters)

		proposal = partial(recom,
						pop_col='P0030001',
						pop_target=sum(self.G.nodes[i]['P0030001'] for i in self.G.nodes())/self.k,
						epsilon=self.population_deviation/2,
						node_repeats=2
						)

		compactness_bound = constraints.UpperBound(
			lambda p: len(p['cut_edges']),
			1.5 * len(initial_partition['cut_edges'])
		)

		pop_constraint = constraints.within_percent_of_ideal_population(initial_partition, self.population_deviation/2)

		my_chain = MarkovChain(
			proposal=proposal,
			constraints=[
				pop_constraint,
				compactness_bound
			],
			accept=accept.always_accept,
			initial_state=initial_partition,
			total_steps=iterations
		)


		min_s = self.upper_bound_s
		print('In GerryChain heuristic, current s and number of minority districts: ', end='')
		max_of_minority_districts = -1
		all_maps = []
		pareto_frontier = []
		obj_vals = []
		for partition in my_chain:
			number_minority_district = 0
			current_s = -1
			for district in range(self.k):
				total_pop_district = 0
				total_pop_minority = 0
				list_of_nodes_in_current_district = []
				for node in partition.graph:
						if partition.assignment[node] == district:
							list_of_nodes_in_current_district.append(node)
							total_pop_district += self.voting_age_population[node]
							total_pop_minority += self.minority_population[node]

				district_subgraph = self.G.subgraph(list_of_nodes_in_current_district)
				district_diamater = nx.diameter(district_subgraph)

				if district_diamater > current_s:
					current_s = district_diamater
				if (total_pop_minority > self.f * total_pop_district):
						number_minority_district += 1

			if current_s < min_s:
				min_s = current_s

			if number_minority_district > max_of_minority_districts:
				max_of_minority_districts = number_minority_district

			print((current_s, number_minority_district),',',sep='',end=' ')
			obj_vals.append([current_s, number_minority_district])
			all_maps.append([partition, current_s, number_minority_district])

		print('Best heuristic solution has the district diameter =', min_s)
		print('Best heuristic solution has # minority districts =', max_of_minority_districts)

		all_maps.sort(key= lambda x: x[1])
		all_maps.sort(key= lambda x: x[2], reverse = True)
		pareto_frontier.append(all_maps[0])
		least_number_of_cut_edges = all_maps[0][1]
		for i in range(1,len(all_maps)):
			if all_maps[i][1] < least_number_of_cut_edges:
				pareto_frontier.append(all_maps[i])
				least_number_of_cut_edges = all_maps[i][1]

		print('Pareto Frontier: ', pareto_frontier)

		optimal_maps = []
		heurisitic_max_s = []
		optimal_minority_districts = []

		i = 0
		for plan in pareto_frontier:
			optimal_maps.append([[i for i in self.G.nodes if plan[0].assignment[i]==j] for j in range(self.k)])
			heurisitic_max_s.append(plan[1])
			optimal_minority_districts.append(plan[2])

		return pareto_frontier

	def fixing_sub_problem(self, option, pass_level, fixed_vertices, skiped_vertices):

		time1 = time.time()
		iterations = 0

		nodes_to_not_be_checked = fixed_vertices + skiped_vertices
		nodes_to_be_checked = list(set(self.G.nodes()) - set(nodes_to_not_be_checked))

		cardinality = len(nodes_to_be_checked)
		if nodes_to_be_checked == []:
			return []
		nodes_to_be_fixed = []
		if nodes_to_be_checked:
			nodes_to_be_checked_is_empty = False
		else:
			nodes_to_be_checkd_is_empty = True
		

		while not nodes_to_be_checked_is_empty:
			time3 = time.time()
			v = nodes_to_be_checked[-1]
			verticies_in_induced_subgraph = list(set(self.G.nodes()) - (set(nodes_to_be_fixed + fixed_vertices)))
			induced_subgraph = self.G.subgraph(verticies_in_induced_subgraph)
			sub_graph = nx.ego_graph(induced_subgraph, v, radius=self.s)

			m = gp.Model()
			m.Params.OutputFlag = 0
			m.Params.TIME_LIMIT = 10
			


			m._t = m.addVars(sub_graph.nodes(), vtype=gp.GRB.BINARY, name='T')

			m._t[v].lb = 1

			if pass_level == 3 or pass_level == 4:
				self.add_lazy_distance_cuts_fixing_problem(m, sub_graph)

			m.setObjective(1, gp.GRB.MAXIMIZE)

			m.addConstr(self.L <= gp.quicksum(self.population[u] * m._t[u] for u in sub_graph.nodes()))

			m.addConstr(self.U >= gp.quicksum(self.population[u] * m._t[u] for u in sub_graph.nodes()))

			###FLOW
			if pass_level == 2 or pass_level == 4:
				sub_DG = nx.DiGraph(sub_graph)
				m._g = m.addVars(sub_DG.edges(), vtype=gp.GRB.CONTINUOUS, name='g')
				m.addConstrs(gp.quicksum(m._g[a, u] for a in sub_graph.neighbors(u)) - gp.quicksum(m._g[u, a] for a in sub_graph.neighbors(u)) == m._t[u] for u in sub_graph.nodes() if u != v)
				m.addConstrs(gp.quicksum(m._g[a, u] for a in sub_graph.neighbors(u)) <= (cardinality - 1) * m._t[u] for u in sub_graph.nodes() if u != v)
				for a in sub_graph.neighbors(v):
					m._g[a, v].ub = 0


			if option == 'minority':
				m.addConstr(gp.quicksum(self.minority_population[u] * m._t[u] for u in sub_graph.nodes()) >= self.f * gp.quicksum(self.voting_age_population[u] * m._t[u] for u in sub_graph.nodes()))
			if option == 'majority':
				m.addConstr(gp.quicksum(self.minority_population[u] * m._t[u] for u in sub_graph.nodes()) <= self.f * gp.quicksum(self.voting_age_population[u] * m._t[u] for u in sub_graph.nodes()))


			m.optimize()

			checked_nodes = []

			if m.status == gp.GRB.TIME_LIMIT:
				print('Time limit reached')
				nodes_to_be_checked.remove(v)
			else:
				if m.status == gp.GRB.OPTIMAL:
					for i in sub_graph.nodes():
						if m._t[i].x > 0.5:
							checked_nodes.append(i)
					


				else:
					checked_nodes.append(v)
					nodes_to_be_fixed.append(v)

				nodes_to_be_checked = [node for node in nodes_to_be_checked if node not in checked_nodes]

			if not nodes_to_be_checked:
				nodes_to_be_checked_is_empty = True

			iterations += 1
			time4 = time.time()
			
		time2 = time.time()

		if option == 'minority':

			self.minority_fixings = list(set(self.minority_fixings + nodes_to_be_fixed))
			print('percent_fixed', str(100 * len(self.minority_fixings) / len(self.G.nodes())), '%')
			print('TOTAL TIME', str(time2 - time1))

		if option == 'majority':

			self.majority_fixings = list(set(self.majority_fixings + nodes_to_be_fixed))
			print('percent_fixed', str(100 * len(self.majority_fixings) / len(self.G.nodes())), '%')
			print('TOTAL TIME', str(time2 - time1))

	def save_fixing_info(self, pass_level, time, draw_map=False):
		output_string = ''

		with open('./results/fixing/txt/' + self.state + '_' + self.parcel_level + '_' + '{:02d}'.format(self.s) + 'pass' + str(pass_level)+ '.txt', 'w') as doc:
			output_string += (str(time) + ' \n')
			for vertex in self.G:
				if vertex in self.minority_fixings:
					output_string += str(vertex) + ' minor\n'
				elif vertex in self.majority_fixings:
					output_string += str(vertex) + ' major\n'
				else:
					output_string += str(vertex) + ' neutr\n'
			doc.write(output_string)
		doc.close()
		if draw_map:
			self.draw_fixing_map()

	def parcel_fixing_iteration(self, pass_level, draw_map=True):

		print('Attempting to read in,', pass_level, 'pass fixing information')
		filename = './results/fixing/txt/'


		filename += self.state + '_' + self.parcel_level + '_' + '{:02d}'.format(self.s)

		if pass_level == 1:
			filename += 'pass1.txt'
		if pass_level == 2:
			filename += 'pass2.txt'
		if pass_level == 3:
			filename += 'pass3.txt'
		if pass_level == 4:
			filename += 'pass4.txt'

		if os.path.exists(filename):
			self.read_fixing_info(filename)
		else:
			time1 = time.time()

			print('No input file for first pass', self.state, self.parcel_level, self.s)

			if pass_level == 4:
				print('Applying minority pass', pass_level, 'fixing subproblem')
				num_of_fixings = len(self.minority_fixings)
				self.fixing_sub_problem('minority', 4, self.minority_fixings, self.majority_fixings)
				num_of_fixings = len(self.minority_fixings) - num_of_fixings
				while num_of_fixings != 0:
					num_of_fixings = len(self.minority_fixings)
					self.fixing_sub_problem('minority', 4, self.minority_fixings, self.majority_fixings)
					num_of_fixings = len(self.minority_fixings) - num_of_fixings


			
			else:
				print('Applying minority pass ', pass_level, ' fixing subproblem')
				self.fixing_sub_problem('minority', pass_level, self.minority_fixings, self.majority_fixings)
			
				time2 = time.time()
				self.save_fixing_info(pass_level, time2 -time1)
			
			
		print('Successfully got', pass_level, 'pass fixing information')
		print('minority', pass_level, 'pass results in ', len(self.minority_fixings), 'fixings out of ', len(self.G.nodes()), 'verticies.')
		
	def apply_parcel_fixing(self, draw_map=False):

		num_fixed = 0

		self.minority_fixings = []
		self.majority_fixings = []
		
		self.parcel_fixing_iteration(1)
		self.parcel_fixing_iteration(2)
		self.parcel_fixing_iteration(3)

		if draw_map:
			self.draw_fixing_map()

		if self.m != False:
			num_fixed = 0
			for vertex in self.minority_fixings:
				for j in self.minority_district_range:
					self.m._X[int(vertex), j].ub = 0
					num_fixed += 1

			print('Applied all minority fixings to fix ', num_fixed, ' out of ', self.m.NumVars, ' variables')

	def read_fixing_info(self, filename):

		fixing_file = open(filename, 'r')
		time = fixing_file.readline()[:-2]
		num_of_fixings = 0


		while True:
			line = fixing_file.readline()


			if not line:
				break

			vertex = int(line[:-7])
			fixing_status = line[-6:-1]

			if vertex in self.minority_fixings:
				continue
			if vertex in self.majority_fixings:
				continue

			if fixing_status == 'minor':
				self.minority_fixings.append(int(vertex))
				num_of_fixings +=1
			if fixing_status == 'major':
				self.majority_fixings.append(int(vertex))

	def create_fixing_csv_file(self):
		
		with open('fixing_info_each_pass.csv', 'w') as doc:
			doc.write('State, parcel_level, s, pass, number of parcels, number of fixings, percent fixed, time, \n')
		
		
		with open('fixing_info.csv', 'w') as doc:
			doc.write('State, parcel level, k, n, m, s, number of fixings, percent fixed, time, \n')


		for file in os.listdir('./results/fixing/txt/'):
			
			
			number_of_parcels = 0
			number_of_minority_fixing_parcels = 0
			number_of_majority_fixing_parcels = 0

			f = open( './results/fixing/txt/' + file, 'r')

			time = f.readline()[:-2]
			state = file[:2]
			pass_level = file[-5:-4]
			tract_level = file[3:9]
			
			if tract_level[-1] == '_':
				tract_level = tract_level[:-1]

			for line in f:
				fixing_status = line[-6:-1]
				number_of_parcels += 1
				if fixing_status == 'minor':
					number_of_minority_fixing_parcels += 1
				if fixing_status == 'major':
					number_of_majority_fixing_parcels += 1

			s = file[-12:-9]

			if s[0] == '_':
				s = s[1:]
			if s[-1] == '_':
				s = s[:-1]


			temp_time = '\"{:,.2f}\"'.format(float(time))
			output_string = state + ',' + tract_level + ',' + s + ',' + pass_level + ',' + '\"{:,}\"'.format(number_of_parcels) + ',' + '\"{:,}\"'.format(number_of_minority_fixing_parcels) + ',' + '{:.2f}'.format(100*number_of_minority_fixing_parcels/number_of_parcels) +  ',' + temp_time + '\n'
			
			with open('fixing_info_each_pass.csv', 'a') as doc:
				doc.write(output_string)		

		with open('fixing_info_each_pass.csv') as file:
			reader = csv.reader(file, delimiter=',')
			
			time = 0
			
			for row in reader:

				if row[3] != '3':
					continue


				state = row[0]
				parcel_level = row[1]
				s = row[2]
				time += float(row[7].replace(',',''))
				number_of_parcels = row[4].replace(',','')
				number_of_parcels_fixed = row[5].replace(',','')
				percent_fixed = row[6]
				

				if row[3] == '3':
					with open('fixing_info.csv', 'a') as doc:
					
						doc.write(state + ',' + parcel_level + ',' + str(congressional_districts[state.upper()]) + ',' + str(len(self.G.nodes)) +',' + str(len(self.G.edges)) + ',' + s  + ',' + number_of_parcels_fixed + ',' + percent_fixed + ',' + str(time) + '\n')

					time = 0

	def add_lazy_distance_cuts(self, model, option, district_labels=True):

		non_fixed_nodes = list(set(self.G.nodes()) - set(self.minority_fixings))

		if non_fixed_nodes:
			induced_graph = self.G.subgraph(non_fixed_nodes)
			power_graph = nx.power(induced_graph, self.s)
			complement_power_graph = nx.complement(power_graph)

			if option != "majority":
				for (u, v) in complement_power_graph.edges():
					for j in self.minority_district_range:
						if district_labels:
							compact_constr = model.addConstr(model._X[u, j] + model._X[v, j] <= model._Z[j])
						else:
							compact_constr = model.addConstr(model._X[u, j] + model._X[v, j] <= 1)

						compact_constr.Lazy = 3

			if option != "minority":
				for (u, v) in complement_power_graph.edges():
					for j in self.majority_district_range:
						if district_labels:
							compact_constr = model.addConstr(model._Y[u, j] + model._Y[v, j] <= model._W[j])
						else:
							compact_constr = model.addConstr(model._Y[u, j] + model._Y[v, j] <= 1)
						compact_constr.Lazy = 3

	def distance_cuts_cliques(self, model, option, district_labels=True):

		non_fixed_nodes = list(set(self.G.nodes()) - set(self.minority_fixings))

		if non_fixed_nodes:
			induced_graph = self.G.subgraph(non_fixed_nodes)
			power_graph = nx.power(induced_graph, self.s)
			complement_power_graph = nx.complement(power_graph)
			cliques = nx.find_cliques(complement_power_graph)

			if option != "majority":
				
				for clique in cliques:
					if district_labels: 
						for j in self.minority_district_range:
							conflict_constraint = self.m.addConstr(gp.quicksum(self.m._X[c, j] for c in clique) <= self.m._Z[j])
					else:
						for j in self.minority_district_range:
							conflict_constraint = self.m.addConstr(gp.quicksum(self.m._X[c, j] for c in clique) <= 1)
				
					conflict_constraint.Lazy = 3
			
			if option != "minority":

				for clique in cliques:
					if district_labels: 
						for j in self.majority_district_range:
							conflict_constraint = self.m.addConstr(gp.quicksum(self.m._Y[c, j] for c in clique) <= self.m._W[j])
					else:
						for j in self.majority_district_range:
							conflict_constraint = self.m.addConstr(gp.quicksum(self.m._Y[c, j] for c in clique) <= 1)

					conflict_constraint.Lazy = 3

	def add_lazy_distance_cuts_fixing_problem(self, model, sub_graph):

		power_graph = nx.power(sub_graph, self.s)

		complement_power_graph = nx.complement(power_graph)

		for (u, v) in complement_power_graph.edges():

			compact_constr = model.addConstr(model._t[u] + model._t[v] <= 1)
			compact_constr.Lazy = 3

	def plotter(self, option, directory):
		matplotlib.use('Agg')
		plt.axis('off')
		plt.grid(False)
		#colors_list = ['w', 'lightcyan','cornflowerblue','lightcoral', 'red', 'firebrick']
		colors_list = ['dimgray', 'tab:cyan', 'lightcoral', 'tab:olive']
		#colors_list = ['white', 'whitesmoke', 'gainsboro', 'silver', 'grey']
		custom_cmap = matplotlib.colors.ListedColormap(colors_list)

		if option == 'minority':
			self.shape_info['minority'] = -1
			color = -1
			for district in self.minority_districts:
				color += 1
				for v in district:
					self.shape_info['minority'][v] = color

			oplot = self.shape_info.plot(column='minority', figsize=(10, 12), legend=False, cmap=custom_cmap).get_figure()
			plt.show()
			oplot.savefig(directory + '_optMinor.png')
			plt.close()

		if option == 'both':
			self.shape_info['district'] = -1
			self.shape_info['minority'] = 3

			color = 0
			for district in self.minority_districts:
				color += 1
				for v in district:
					self.shape_info['minority'][v] = color
				
			
			oplot = self.shape_info.plot(column='minority', figsize = (10,12), legend=False, cmap=custom_cmap).get_figure()
			plt.axis('off')
			oplot.savefig(directory + '_optMinor.png')
			plt.close()

			color = 1
			for district in self.districts:
				if district in self.minority_districts:
					for v in district:
						self.shape_info['district'][v] = 0
					continue
				if 1 in district:
					for v in district:
						self.shape_info['district'][v] = 1
					continue
				color += 1
				for v in district:
					self.shape_info['district'][v] = color
			
				
			oplot = self.shape_info.plot(column='district', figsize = (10,12), legend=False, cmap=custom_cmap).get_figure()
			plt.axis('off')
			oplot.savefig(directory + '_opt.png')
			plt.close()

	def draw_fixing_map(self):

		
		colors_list = ['lightgray','b','c','g']
		custom_cmap = matplotlib.colors.ListedColormap(colors_list)

		minority_map_verticies = []
		local_copy_shape_info_big_map = self.shape_info.copy()
		local_copy_shape_info_big_map['fixing_map'] = 'Vertex never fixed'
		pass_level_to_label_dict = {1:'Fixed in first pass', 2:'Fixed in second pass', 3:'Fixed in third pass'}

		for pass_level in range(3,0,-1):
			if os.path.exists('./results/fixing/txt/' + self.state + '_' + self.parcel_level + '_' + '{:02d}'.format(self.s) + 'pass' + str(pass_level) + '.txt'):
				fixing_file = open('./results/fixing/txt/' + self.state + '_' + self.parcel_level + '_' + '{:02d}'.format(self.s) + 'pass' + str(pass_level) + '.txt', 'r')
				white_nodes = []
				black_nodes = []
				gray_nodes = []

				while True:
					line = fixing_file.readline()

					if not line:
						break

					vertex = line[:-7]
					fixing_status = line[-6:-1]

					if fixing_status == 'minor':
						black_nodes.append(int(vertex))
					if fixing_status == 'major':
						white_nodes.append(int(vertex))
					if fixing_status == 'neutr':
						gray_nodes.append(int(vertex))

				local_copy_shape_info = self.shape_info.copy()
				local_copy_shape_info['fixing_map'] = 'Vertex not fixed'

				my_categories = []

				for v in black_nodes:
					category = 'Parcel must join majority district'
					local_copy_shape_info['fixing_map'][v] = category
					local_copy_shape_info_big_map['fixing_map'][v] = pass_level_to_label_dict[pass_level]
					
					if category not in my_categories:
						my_categories.append(category)

				for v in gray_nodes:

					category = 'Parcel can join either majority or minority district'
					local_copy_shape_info['fixing_map'][v] = category
				
					if category not in my_categories:
						my_categories.append(category)

				for v in white_nodes:
					category = 'Parcel can must join minority district'
					local_copy_shape_info['fixing_map'][v] = category

					if category not in my_categories:
						my_categories.append(category)

				plt.style.use('classic')
				
				cmap = matplotlib.colors.ListedColormap(['0.9', 'gray', 'black', '0.3'])
				
				oplot = local_copy_shape_info.plot(column='fixing_map',categorical=True, cmap = cmap, legend_kwds={'loc': 'best'}, categories=my_categories, legend=True, figsize = (10,12)).get_figure()
				plt.axis('off')
				plt.title('Pass level ' + str(pass_level) + ' for state ' + self.state + ' s: ' + str(self.s))
				plt.grid(False)

				oplot.patch.set_facecolor('xkcd:white')

				oplot.savefig('./results/fixing/img/' + self.state + '_' + self.parcel_level + '_' + str(self.s) + 'pass' + str(pass_level) + '.png')

				plt.style.use('classic')
				oplot = local_copy_shape_info.plot(column='fixing_map',categorical=True, cmap = cmap, categories=my_categories, figsize = (10,12)).get_figure()
				plt.axis('off')
				plt.title('Pass level ' + str(pass_level) + ' for state ' + self.state + ' s: ' + str(self.s))
				plt.grid(False)
				oplot.patch.set_facecolor('xkcd:white')

				oplot.savefig('./results/fixing/img/no_legend_pass_' + str(pass_level) + '_' + self.state + '_' + self.parcel_level + '_' + str(self.s) + '.png')


		plt.style.use('classic')
		oplot = local_copy_shape_info_big_map.plot(column='fixing_map',categorical=True, cmap = custom_cmap, legend_kwds={'loc': 'best'}, categories=['Vertex never fixed', 'Fixed in first pass', 'Fixed in second pass', 'Fixed in third pass'], legend=True, figsize = (10,12)).get_figure()

		plt.axis('off')
		plt.title('State: ' + self.state + ', s: ' + str(self.s) + ', when which parcel got fixed')
		plt.grid(False)
		oplot.patch.set_facecolor('xkcd:white')
		plt.show()
		oplot.savefig('./results/fixing/img/' + self.state + '_' + str(self.s) + '.png')

		plt.style.use('classic')
		oplot = local_copy_shape_info_big_map.plot(column='fixing_map',categorical=True, cmap = 'Set1', categories=['Vertex never fixed', 'Fixed in first pass', 'Fixed in second pass', 'Fixed in third pass'], figsize = (10,12)).get_figure()
		plt.axis('off')
		plt.title('State: ' + self.state + ', s: ' + str(self.s) + ', when which parcel got fixed')
		plt.grid(False)
		oplot.patch.set_facecolor('xkcd:white')
		oplot.savefig('./results/fixing/img/no_legend_' + self.state + '_' + str(self.s) + '.png')

	def draw_symmetry_map(self):
		max_independent_set = compuete_stability_number('read', self.G, self.state, self.k, self.s)
		number_of_symmetry_breaking_allowed = min(len(max_independent_set), len(self.majority_district_range))

		local_copy_shape_info = self.shape_info.copy()
		my_categories = ['parcel in independent_set', 'parcels reachable by parcel in independent_set', 'parcels not reachable in independent set']
		
		for fixing in range(number_of_symmetry_breaking_allowed):
			local_copy_shape_info['fixing_map'] = my_categories[2]
			v = max_independent_set[0]
			max_independent_set.pop(0)
			local_copy_shape_info['fixing_map'][v] = my_categories[0]


			power_graph = nx.power(self.G, self.s)
			for vertex in power_graph.neighbors(v):
				local_copy_shape_info['fixing_map'][vertex] = my_categories[1]


			plt.style.use('classic')
				
			
			cmap = matplotlib.colors.ListedColormap(['0.95', '0.7', '0.25'])
			
			oplot = local_copy_shape_info.plot(column='fixing_map',categorical=True, cmap = cmap, legend_kwds={'loc': 'best'}, categories=my_categories, legend=False, figsize = (10,12)).get_figure()
			plt.axis('off')
			plt.grid(False)
			oplot.patch.set_facecolor('xkcd:white')
			plt.show()

	def add_symmetry_constraints(self):
		max_independent_set = list(compuete_stability_number('read', self.G, self.state, self.k, self.s))
		number_of_symmetry_breaking_allowed = min(len(max_independent_set), len(self.majority_district_range))
		
		i = self.k - 1
		for fixing in range(number_of_symmetry_breaking_allowed):
			v = max_independent_set[0]
			self.m.addConstr(gp.quicksum(self.m._X[v,j] for j in self.minority_district_range) + self.m._Y[v,i] == 1)
			max_independent_set.pop(0)
			i -= 1

	def minority_district_ONLY_model(self, s, relax=False, continuity=False):
		self.s = s
		time1 = time.time()
		self.m = gp.Model()
		self.m.Params.TIME_LIMIT = self.model_time_limit
		
		if relax:
			self.m._X = self.m.addVars(self.G.nodes(), self.minority_district_range, vtype=gp.GRB.CONTINUOUS, name='X')
		else:
			self.m._X = self.m.addVars(self.G.nodes(), self.minority_district_range, vtype=gp.GRB.BINARY, name='X')

		self.m._Z = self.m.addVars(self.minority_district_range, vtype=gp.GRB.CONTINUOUS, ub=1, name='Z')
		
		if s:
			
			self.apply_parcel_fixing()

		#1a
		self.m.setObjective(gp.quicksum(self.m._Z[j] for j in self.minority_district_range), gp.GRB.MAXIMIZE)

		#1b
		self.m.addConstrs(gp.quicksum(self.m._X[v, j] for j in self.minority_district_range) <= 1 for v in self.G.nodes())

		#1c
		self.m.addConstrs(self.L * self.m._Z[j] <= gp.quicksum(self.population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in self.minority_district_range)
		#1c
		self.m.addConstrs(self.U * self.m._Z[j] >= gp.quicksum(self.population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in self.minority_district_range)
		#1d
		self.m.addConstrs(gp.quicksum(self.minority_population[v] * self.m._X[v, j] for v in self.G.nodes()) >= self.f * gp.quicksum(self.voting_age_population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in self.minority_district_range)


		#1f
		self.m.addConstrs(self.m._X[v, j] <= self.m._Z[j] for v in self.G.nodes() for j in self.minority_district_range)

		#g
		self.m.addConstrs(self.m._Z[i] >= self.m._Z[i + 1] for i in range(self.k_minority - 1))
		self.m.addConstr(self.m._Z[0] <= 1)

		if continuity == 'callback':
			self.m.Params.lazyConstraints = 1
			self.m._option = 'minority'
			self.m._district_vars_exist = True
			self.m._G = self.G
			self.m._k = self.k
			self.m._minority_district_range = self.minority_district_range
			self.m.optimize(continuous_ONLY_callback)

		elif continuity == 'flow':

			M = most_possible_nodes_in_one_district(self.population, self.U) - 1
			DG = nx.to_directed(self.G)

			if relax:
				self.m._R = self.m.addVars(self.G.nodes(), self.minority_district_range, vtype=GRB.CONTINUOUS)
			else:
				self.m._R = self.m.addVars(self.G.nodes(), self.minority_district_range, vtype=GRB.BINARY)

			# t[i,j] = amount of flow generated at node i of type j
			t = self.m.addVars(DG.nodes, self.minority_district_range, vtype=GRB.CONTINUOUS)

			# g[j,u,v] = amount of flow sent across arc uv of type j
			g = self.m.addVars(self.minority_district_range, DG.edges, vtype=GRB.CONTINUOUS)

			# compute big-M


			# the following constraints are weaker than some in the orbitope EF
			#if symmetry != 'orbitope':
			self.m.addConstrs( gp.quicksum(self.m._R[i,j] for i in DG.nodes)==1 for j in self.minority_district_range )

			self.m.addConstrs( self.m._R[i,j] <= self.m._X[i,j] for i in DG.nodes for j in self.minority_district_range )

			# flow can only be generated at roots
			self.m.addConstrs( t[i,j] <= (M+1)*self.m._R[i,j] for i in DG.nodes for j in self.minority_district_range )

			# flow balance
			self.m.addConstrs( t[i,j] - self.m._X[i,j] == gp.quicksum(g[j,i,u]-g[j,u,i] for u in DG.neighbors(i)) for i in DG.nodes for j in self.minority_district_range )

			# flow type j can enter vertex i only if (i is assigned to district j) and (i is not root of j)
			self.m.addConstrs( gp.quicksum(g[j,u,i] for u in DG.neighbors(i)) <= M*(self.m._X[i,j]-self.m._R[i,j]) for i in DG.nodes for j in self.minority_district_range )

			self.m.optimize()
		else:
			self.m.optimize()

		if self.m.status == gp.GRB.OPTIMAL:

			self.minority_districts = [[] for i in range(int(self.m.getAttr(gp.GRB.Attr.ObjVal)))]
			self.num_minority_districts = len(self.minority_districts)



			for district_number in range(self.num_minority_districts):
				for vertex in self.G.nodes():
					if self.m._X[vertex, district_number].x > 0.5:
						self.minority_districts[district_number].append(vertex)


			#For root relax
			self.BB_node = self.m.NodeCount
			self.output = int(self.m.getAttr(gp.GRB.Attr.ObjVal))
			self.time_out = False

		elif self.m.status == gp.GRB.TIME_LIMIT:
			self.num_minority_districts = self.m.getAttr(gp.GRB.Attr.ObjVal)
			self.time_out = True
		else:
			self.num_minority_districts = 'Model Status: ' + str(self.m.status)
			self.time_out = False
		time2 = time.time()
		
		self.minority_district_ONLY_model_time = time2 - time1

	def labeling_model(self, reduced_model, parcel_fixing, symmetry_breaking_constraint):

			s = self.s
			self.m = gp.Model()
			self.m.Params.TIME_LIMIT = self.model_time_limit

			
			if reduced_model:
				self.m._X = self.m.addVars(self.G.nodes(), self.minority_district_range, vtype=gp.GRB.BINARY, name='X')
				self.m._Z = self.m.addVars(self.minority_district_range, vtype=gp.GRB.BINARY, name='Z')
			else:
				self.m._X = self.m.addVars(self.G.nodes(), range(self.k), vtype=gp.GRB.BINARY, name='X')
				self.m._Z = self.m.addVars(range(self.k), vtype=gp.GRB.BINARY, name='Z')

			self.m._Y = self.m.addVars(self.G.nodes(), range(self.k), vtype=gp.GRB.BINARY, name='Y')
			self.m._W = self.m.addVars(range(self.k), vtype=gp.GRB.BINARY, name='W')

			
			if reduced_model:
				for i in self.majority_district_range:
					self.m._W[i].lb = 1

			if symmetry_breaking_constraint:
				self.add_symmetry_constraints()
				
			
			if parcel_fixing:
				self.apply_parcel_fixing()

			if s:
				distance_constraints = 'clique'
				#distance_constraints = 'pairwise_lazy'
				
				if distance_constraints == 'pairwise_lazy':
					self.add_lazy_distance_cuts(self.m, 'both')
				elif distance_constraints == 'clique':
					self.distance_cuts_cliques(self.m, 'both')

			# Base Constraints
			if reduced_model:
				self.m.setObjective(gp.quicksum(self.m._Z[i] for i in self.minority_district_range), gp.GRB.MAXIMIZE)

				self.m.addConstrs(gp.quicksum(self.m._X[v, j] for j in self.minority_district_range) + gp.quicksum(self.m._Y[v, j] for j in range(self.k)) == 1 for v in self.G.nodes())

				self.m.addConstrs(self.m._Z[j] + self.m._W[j] == 1 for j in self.minority_district_range)

				self.m.addConstrs(self.m._X[v, j] <= self.m._Z[j] for v in self.G.nodes() for j in self.minority_district_range)
			else:
				self.m.setObjective(gp.quicksum(self.m._Z[i] for i in range(self.k)), gp.GRB.MAXIMIZE)

				self.m.addConstrs(gp.quicksum(self.m._X[v, j] for j in range(self.k)) + gp.quicksum(self.m._Y[v, j] for j in range(self.k)) == 1 for v in self.G.nodes())

				self.m.addConstrs(self.m._Z[j] + self.m._W[j] == 1 for j in range(self.k))

				self.m.addConstrs(self.m._X[v, j] <= self.m._Z[j] for v in self.G.nodes() for j in range(self.k))

			self.m.addConstrs(self.m._Y[v, j] <= self.m._W[j] for v in self.G.nodes() for j in range(self.k))

			# Population Constraints
			if reduced_model:
				self.m.addConstrs(self.L * self.m._Z[j] <= gp.quicksum(self.population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in self.minority_district_range)

				self.m.addConstrs(self.U * self.m._Z[j] >= gp.quicksum(self.population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in self.minority_district_range)

			else:
				self.m.addConstrs(self.L * self.m._Z[j] <= gp.quicksum(self.population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in range(self.k))

				self.m.addConstrs(self.U * self.m._Z[j] >= gp.quicksum(self.population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in range(self.k))


			self.m.addConstrs(self.L * self.m._W[j] <= gp.quicksum(self.population[v] * self.m._Y[v, j] for v in self.G.nodes()) for j in range(self.k))

			self.m.addConstrs(self.U * self.m._W[j] >= gp.quicksum(self.population[v] * self.m._Y[v, j] for v in self.G.nodes()) for j in range(self.k))

			if reduced_model:
				self.m.addConstrs(gp.quicksum(self.minority_population[v] * self.m._X[v, j] for v in self.G.nodes()) >= self.f * gp.quicksum(self.voting_age_population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in self.minority_district_range)
			else:
				self.m.addConstrs(gp.quicksum(self.minority_population[v] * self.m._X[v, j] for v in self.G.nodes()) >= self.f * gp.quicksum(self.voting_age_population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in range(self.k))
			
			self.m.addConstrs(gp.quicksum(self.voting_age_population[v] - self.minority_population[v] * self.m._Y[v, j] for v in self.G.nodes()) >= (1 - self.f) * gp.quicksum(self.voting_age_population[v] * self.m._Y[v, j] for v in self.G.nodes()) for j in range(self.k))

			if reduced_model:
				self.m.addConstrs(self.m._Z[i] >= self.m._Z[i + 1] for i in range(self.k_minority - 1))
			else:
				self.m.addConstrs(self.m._Z[i] >= self.m._Z[i + 1] for i in range(self.k - 1))
			
			self.m._k_minority = self.k_minority
			self.m._G = self.G
			self.m._k = self.k
			self.m._option = 'both'
			self.m._s = s


			self.m.Params.lazyConstraints = 1
			self.m._minority_district_range = self.minority_district_range
			self.m._majority_district_range = range(self.k)
			
			self.m._district_vars_exist = True
			self.m._L = self.L
			self.m._population = self.population
			self.m._n = len(self.G.nodes())

			if False:
				counter = 0
				for minority_district in minority_districts:
					for node in minority_district:
						self.m._X[node, counter].start = 1
					counter += 1
				for majority_district in majority_districts:
					for node in majority_district:
						self.m._Y[node, counter].start = 1
					counter += 1

			#self.m.optimize()
			self.m.optimize(compactness_callback)
			
			self.model_runtime = self.m.runtime


			if self.m.status == gp.GRB.OPTIMAL:

				self.minority_districts = [[] for i in range(int(self.m.getAttr(gp.GRB.Attr.ObjVal)))]
				self.districts = [[] for i in range(self.k)]

				self.num_minority_districts = len(self.minority_districts)
				for district_number in range(self.num_minority_districts):
					for vertex in self.G.nodes():
						if self.m._X[vertex, district_number].x > 0.5:
							self.minority_districts[district_number].append(vertex)
							self.districts[district_number].append(vertex)
				for district_number in range(self.num_minority_districts, self.k):
					for vertex in self.G.nodes():
						if self.m._Y[vertex, district_number].x > 0.5:
							self.districts[district_number].append(vertex)
			elif self.m.status == 3:
				self.num_minority_districts = 'Infeasible'
			elif self.m.status == 9:
				self.num_minority_districts = 'Timeout'
			else:
				self.num_minority_districts = 'Model Status: ' + str(self.m.status)

	def prescribed_labeling_model(self, warm_start=False, symmetry_breaking_constraint=True, include_majority_variables=True):

			self.m = gp.Model()
			self.m.Params.TIME_LIMIT = self.model_time_limit

			self.m._X = self.m.addVars(self.G.nodes(), self.minority_district_range, vtype=gp.GRB.BINARY, name='X')
			if include_majority_variables:
				self.m._Y = self.m.addVars(self.G.nodes(), self.majority_district_range, vtype=gp.GRB.BINARY, name='Y')


			if warm_start:
				counter = 0
				for minority_district in warm_start[0]:
					for node in minority_district:
						self.m._X[node, counter].start = 1
					counter += 1
				for majority_district in warm_start[1]:
					for node in majority_district:
						if include_majority_variables:
							self.m._Y[node, counter].start = 1
					counter += 1
				#self.m.write('out.mst')

			if symmetry_breaking_constraint:
				if include_majority_variables:
					self.add_symmetry_constraints()

			self.apply_parcel_fixing()
			if include_majority_variables:
				self.add_lazy_distance_cuts(self.m, 'both', False)
			else:
				self.add_lazy_distance_cuts(self.m, 'minority', False)
			
			#max_independent_set = compuete_stability_number('read', self.G, self.state, self.k, self.s)
			#number_of_symmetry_breaking_allowed = min(len(max_independent_set), len(self.majority_district_range))
			
			#i = self.k - 1
			#for fixing in range(number_of_symmetry_breaking_allowed):
			#	v = max_independent_set[0]
			#	self.m.addConstr(gp.quicksum(self.m._X[v,j] for j in self.minority_district_range) + self.m._Y[v,i] == 1)
			#	max_independent_set.pop(0)
			#	i -= 1

			# Base Constraints
			self.m.setObjective(self.k_minority, gp.GRB.MAXIMIZE)
			#self.m.setObjective(gp.quicksum(self.m._X))

			if include_majority_variables:
				self.m.addConstrs(gp.quicksum(self.m._X[v, j] for j in self.minority_district_range) + gp.quicksum(self.m._Y[v, j] for j in self.majority_district_range) == 1 for v in self.G.nodes())



			# Population Constraints
			self.m.addConstrs(self.L <= gp.quicksum(self.population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in self.minority_district_range)

			self.m.addConstrs(self.U >= gp.quicksum(self.population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in self.minority_district_range)

			if include_majority_variables:
				self.m.addConstrs(self.L <= gp.quicksum(self.population[v] * self.m._Y[v, j] for v in self.G.nodes()) for j in self.majority_district_range)

			if include_majority_variables:
				self.m.addConstrs(self.U >= gp.quicksum(self.population[v] * self.m._Y[v, j] for v in self.G.nodes()) for j in self.majority_district_range)

			self.m.addConstrs(gp.quicksum(self.minority_population[v] * self.m._X[v, j] for v in self.G.nodes()) >= self.f * gp.quicksum(self.voting_age_population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in self.minority_district_range)

			if include_majority_variables:
				self.m.addConstrs(gp.quicksum(self.population[v] - self.minority_population[v] * self.m._Y[v, j] for v in self.G.nodes()) >= (1 - self.f) * gp.quicksum(self.voting_age_population[v] * self.m._Y[v, j] for v in self.G.nodes()) for j in self.majority_district_range)


			self.m._k_minority = self.k_minority
			self.m._G = self.G
			self.m._k = self.k
			if include_majority_variables:
				self.m._option = 'both'
			else:
				self.m_option = 'minority'
			self.m._s = self.s


			self.m.Params.lazyConstraints = 1
			self.m._minority_district_range = self.minority_district_range
			self.m._majority_district_range = self.majority_district_range

			self.m._district_vars_exist = False
			self.m._L = self.L
			self.m._population = self.population
			self.m._n = len(self.G.nodes())

			self.m.optimize(compactness_callback)


			self.model_runtime = self.m.runtime

			if self.m.status == gp.GRB.OPTIMAL:
				
				self.minority_districts = [[] for i in self.minority_district_range]
				self.districts = [[] for i in range(self.k)]

				for district_number in self.minority_district_range:
					for vertex in self.G.nodes():
						if self.m._X[vertex, district_number].x > 0.5:
							self.minority_districts[district_number].append(vertex)
							self.districts[district_number].append(vertex)
				
				if include_majority_variables:
					for district_number in self.majority_district_range:
						for vertex in self.G.nodes():
							if self.m._Y[vertex, district_number].x > 0.5:
								self.districts[district_number].append(vertex)
				self.num_minority_districts = len(self.minority_districts)
			
			else:
				self.num_minority_districts = 'Model Status: ' + str(self.m.status)

	def generate_gerrychain_run(self, number_of_gerrychain_iter):

		
		time1 = time.time()
		gerrychain_results = instance.run_GerryChain_heuristic(number_of_gerrychain_iter)
		time2 = time.time()
		gerrychain_time = '{0:.2f}'.format(time2-time1)
		gerrychain_output_string = instance.state + ',' + instance.parcel_level + ',' + str(number_of_gerrychain_iter) + ',' + gerrychain_time + ','


		for result in gerrychain_results:

			partition = result[0]
			gerrychain_s = result[1]
			gerrychain_num_minority_district = result[2]

			gerrychain_output_string += '(' + str(gerrychain_s) + '|' + str(gerrychain_num_minority_district) + ') '


			minority_districts_from_gerrychain = []
			majority_districts_from_gerrychain = []


			for district in range(instance.k):
				total_pop_district = 0
				total_pop_minority = 0
				list_of_nodes_in_current_district = []
				for node in partition.graph:
						if partition.assignment[node] == district:
							list_of_nodes_in_current_district.append(node)

							total_pop_district += instance.voting_age_population[node]
							total_pop_minority += instance.minority_population[node]

				district_subgraph = instance.G.subgraph(list_of_nodes_in_current_district)
				district_diamater = nx.diameter(district_subgraph)
				if total_pop_minority > self.f * total_pop_district:
					minority_districts_from_gerrychain.append(list_of_nodes_in_current_district)
				else:
					majority_districts_from_gerrychain.append(list_of_nodes_in_current_district)

			gerrychain_output = [minority_districts_from_gerrychain, majority_districts_from_gerrychain]

		
			with open('./results/gerrychain/' + self.state + '_' + self.parcel_level + '_(' + '{:02d}'.format(gerrychain_s) + ',' + '{:02d}'.format(gerrychain_num_minority_district) + ').pckl', 'wb') as doc:
				pickle.dump(gerrychain_output, doc)
		gerrychain_output_string += '\n'
		
		if not os.path.exists('gerrychain.csv'):
			with open('gerrychain.csv', 'w') as doc:
				doc.write('State, parcel level, gerrychain iterations, time, pareto frontier \n')

		with open('gerrychain.csv', 'a') as doc:
			doc.write(gerrychain_output_string)
		
	def make_sure_preprocess_exists(self):

		result = compuete_stability_number('read', self.G, self.state, self.k, self.s)
		
		self.minority_fixings = []
		self.majority_fixings = []

		self.parcel_fixing_iteration(1)
		self.parcel_fixing_iteration(2)
		self.parcel_fixing_iteration(3)

	def generate_instance_info(self):
		if not os.path.exists('instance_info.csv'):
			with open('instance_info.csv', 'w') as doc:
				doc.write('State, parcel, k, n, m, diam, % Black pop, current s, split parcels, current # minority districts')
		
		with open('instance_info.csv', 'a') as doc:
			current_s, split_number = self.get_current_plan_stat()
			string = '\n'
			string += self.state + ', ' + self.parcel_level + ', ' + str(self.k) + ', '+ str(len(self.G.nodes())) + ', ' + str(len(self.G.edges())) + ', ' + str(nx.diameter(self.G)) + ', ' + format(100 * sum(self.minority_population)/sum(self.voting_age_population), '.2f') + ',' + str(current_s) + ',' +  str(split_number) + self.error_message + ', '
			
			doc.write(string)
			doc.close()

	def generate_upper_bound_minority_table(self):
		
		if not os.path.exists('upper_bound_minority_table.csv'):
			with open('upper_bound_minority_table.csv', 'w') as doc:
				doc.write('State, parcel, k, n, m, objective, time (s), x & z fixed (%)')

		instance.k_minority = instance.k
		instance.minority_district_range = range(instance.k)
		instance.find_upper_bound_minority_districts(continuity=False, relax=False)
		instance.table_3_first_pass = instance.output
		instance.table_3_first_pass_time = instance.minority_district_ONLY_model_time
		instance.minority_district_range = range(instance.k)

		with open('upper_bound_minority_table.csv', 'a') as doc:
				string = '\n'
				string += self.state + ',' 
				string += self.parcel_level + ',' 
				string += str(self.k) + ',' 
				string += str(len(self.G.nodes)) + ',' 
				string += str(len(self.G.edges)) + ',' 
				string += format(self.table_3_first_pass, '.2f') + ','
				string += format(self.table_3_first_pass_time, '.2f') + ','
				string += format(self.table_3_first_pass/self.k * 100, '.2f') + ',' 

				doc.write(string)
				doc.close()

	def generate_lower_bound_s_table(self):
		if not os.path.exists('lower_bound_s_table.csv'):
			with open('lower_bound_s_table.csv', 'w') as doc:
				doc.write('State, land parcel, k, n, m, l_s, Iter, time (s)')
		
		instance.compute_lower_bound_s()
		
		

		with open('lower_bound_s_table.csv', 'a') as doc:
			string = '\n'
			string += self.state + ','
			string += self.parcel_level  + ','
			string += str(self.k) + ','
			string += str(len(self.G.nodes())) + ',' 
			string += str(len(self.G.edges())) + ',' 
			string += str(self.lower_bound_s) + ','
			string += str(self.num_s_iter) + ',' 
			string += format(self.time_for_s_lower_bound, '.2f')

			doc.write(string)
			doc.close()

	def generate_fixing_info(self):
		
		condition = True
		self.s = self.lower_bound_s - 1
		while condition:
			self.s += 1
			print('Current s: ', str(self.s))
			self.make_sure_preprocess_exists()
			if self.minority_fixings == []:
				condition = False
			if self.s == self.upper_bound_s:
				condition = False

	def generate_lower_bound_s_table(self):
		
		if not os.path.exists('lower_bound_s_table.csv'):
			with open('lower_bound_s_table.csv', 'w') as doc:
				doc.write('State, parcel, lower_bound_s, iterations, time')

		instance.compute_lower_bound_s()
		
		with open('lower_bound_s_table.csv', 'a') as doc:
			string = '\n'
			string += self.state + ', ' + self.parcel_level + ', ' + str(self.lower_bound_s) + ', ' + str(self.num_s_iter) + ', ' + format(self.time_for_s_lower_bound, '.2f') + ',' + self.error_message + ', '
			doc.write(string)
			doc.close()
	'''
	def generate_table_7(self):
		
		if not os.path.exists('table_7.csv'):
			with open('table_7.csv', 'w') as doc:
				doc.write('State, parcel, n, m, s, reduced_obj, reduced_time, prescribed_obj, prescribed time \n')

		try_higher_s = True

		while try_higher_s:
			

			self.make_sure_preprocess_exists()
			#
			#self.labeling_model(False, False, False)
			#base_model_objective = self.num_minority_districts
			#base_model_time = '{0:.2f}'.format(self.model_runtime)
			#
			#self.labeling_model(True, False, False)
			#reduced_model_objective_plain = self.num_minority_districts
			#reduced_model_time_plain = '{0:.2f}'.format(self.model_runtime)
			#
			self.labeling_model(True, True, True)
			reduced_model_objective = self.num_minority_districts
			reduced_model_time = '{0:.2f}'.format(self.model_runtime)
			#
			#self.labeling_model(True, True, True)
			#reduced_model_objective_with_symmetry_parcel_fixing = self.num_minority_districts
			#reduced_model_time_with_symmetry_parcel_fixing = '{0:.2f}'.format(self.model_runtime)
			#
			#self.labeling_model(True, True, False)
			#reduced_model_no_sym = self.num_minority_districts
			#reduced_model_time_no_sym = '{0:.2f}'.format(self.model_runtime)
			#			
			prescribed_model_objective, prescribed_model_time = self.optimize_prescribed_labeling_model([])
			prescribed_model_time = '{0:.2f}'.format(prescribed_model_time)



			if prescribed_model_objective == 'Timeout':
				try_higher_s = False
			if prescribed_model_objective == self.k_minority:
				try_higher_s = False
	

			with open('table_7.csv', 'a') as doc:
				doc.write(self.state + ',' + self.parcel_level + ',' + str(len(self.G.nodes)) + ',' + str(len(self.G.edges)) + ',' + str(self.s) + ','  + str(reduced_model_objective) + ',' + reduced_model_time + ',' + str(prescribed_model_objective) + ',' + prescribed_model_time + '\n')
								
			self.s += 1
	'''					
	def generate_point_check_table(self):
		if not os.path.exists('point_check.csv'):
			with open('point_check.csv', 'w') as doc:
				doc.write('State, parcel level, s, number of minority districts, model result, model time \n')
			
		values_of_k_minority_to_investigate = range(self.k_minority, 0, -1)

		for i in values_of_k_minority_to_investigate:
			if i == 0:
				continue

			self.k_minority = i
			self.minority_district_range = range(self.k_minority)
			self.majority_district_range = range(self.k_minority, self.k)
			condition = True

			while condition == True:

				time1 = time.time()
				self.prescribed_labeling_model()
				time2 = time.time()
				MIP_time = '{0:.2f}'.format(time2-time1)


				MIP_output_string = self.state + ',' + self.parcel_level + ',' + str(self.s) + ',' + str(self.k_minority) + ','

				if self.m.status == 2:
					model_result = 'Feasible'
					condition = False
				elif self.m.status == 3:
					model_result = 'Infeasible'
					condition = True
					self.s += 1
				elif self.m.status == 9:
					model_result = 'Timeout'
					condition = False
				else:
					model_result = 'unanticpated termination gurobicode: ' + str(self.m.status)
					condition = False

				MIP_output_string += model_result + ',' + MIP_time + ',' + '\n'
				with open('point_check.csv', 'a') as doc:
					doc.write(MIP_output_string)

				if self.s >= 10:
					condition=False

	def evaluate_gerrychain(self):

		if not os.path.exists('table_8.csv'):
			with open('table_8.csv', 'w') as doc:
				doc.write('State, parcel level, s, number of minority districts, model result, model time \n')
		
		values_of_k_minority_to_investigate = range(self.k_minority)

		
		
		for file in os.listdir('results/gerrychain'):
			gerrychain_state = file[:2]
			if file[3:9] == 'tract_':
				gerrychain_parcel_level = 'tract'
			if file[3:9] == 'county':
				gerrychain_parcel_level = 'county'

			gerrychain_s = int(file[-11:-9])
			gerrychain_num_minority_district = int(file[-8:-6])

			if self.state == gerrychain_state and self.parcel_level == gerrychain_parcel_level:

				with open('results/gerrychain/' + file, 'rb') as f:
					gerrychain_data = pickle.load(f)
					minority_districts_from_gerrychain = gerrychain_data[0]
					majority_districts_from_gerrychain = gerrychain_data[1]

				self.k_minority = gerrychain_num_minority_district
				self.s = gerrychain_s
				
				next_s_to_run = gerrychain_s
				try_to_improve_s = False
				try_to_improve_k = True

				while try_to_improve_k:
					self.k_minority += 1

					if self.k_minority not in values_of_k_minority_to_investigate:
						print('here3')
						try_to_improve_k = False
						continue

					self.minority_district_range = range(self.k_minority)
					self.majority_district_range = range(self.k_minority, self.k)

					time1 = time.time()
					self.prescribed_labeling_model([minority_districts_from_gerrychain, []])
					time2 = time.time()

					MIP_time = '{0:.2f}'.format(time2-time1)

					MIP_output_string = self.state + ',' + self.parcel_level + ',' + str(self.s) + ',' + str(self.k_minority) + ','

					if self.m.status == 2:
						model_result = 'Feasible'
						try_to_improve_k = True
						
					elif self.m.status == 3:
						model_result = 'Infeasible'
						try_to_improve_k = False

					elif self.m.status == 9:
						model_result = 'Timeout'
						try_to_improve_k = False

					else:
						model_result = 'unanticpated termination gurobicode: ' + str(self.m.status)
						try_to_improve_k = False


					MIP_output_string += model_result + ',' + MIP_time + ',' + '\n'
					with open('table_8.csv', 'a') as doc:
						doc.write(MIP_output_string)

				while try_to_improve_s:
		
					self.s = next_s_to_run

					time1 = time.time()
					self.prescribed_labeling_model([minority_districts_from_gerrychain, majority_districts_from_gerrychain])
					time2 = time.time()

					MIP_time = '{0:.2f}'.format(time2-time1)
					MIP_output_string = self.state + ',' + self.parcel_level + ',' + str(self.s) + ',' + str(self.k_minority) + ','

					if self.m.status == 2:
						model_result = 'Feasible'
						condition = True
						next_s_to_run -= 1
					elif self.m.status == 3:
						model_result = 'Infeasible'
						condition = False
					elif self.m.status == 9:
						model_result = 'Timeout'
						condition = False
					else:
						model_result = 'unanticpated termination gurobicode: ' + str(self.m.status)
						condition = False


					MIP_output_string += model_result + ',' + MIP_time + ',' + '\n'
					with open('table_8.csv', 'a') as doc:
						doc.write(MIP_output_string)

	def optimize_prescribed_labeling_model(self, warm_start):
		k_minority = self.k_minority
		optimization_time = 0
		while True:
								
			self.make_sure_preprocess_exists()
			self.prescribed_labeling_model()
			optimization_time += self.model_runtime

			
			if self.m.status == 2:
				
				model_objective = self.num_minority_districts
				self.k_minority = k_minority
				self.minority_district_range = range(self.k_minority)
				self.majority_district_range = range(self.k_minority, self.k)
				return model_objective, optimization_time

			elif self.m.status == 3:
				
				self.k_minority -= 1
				self.minority_district_range = range(self.k_minority)
				self.majority_district_range = range(self.k_minority, self.k)
				if self.k_minority < 0:
					self.k_minority = k_minority
					self.minority_district_range = range(self.k_minority)
					self.majority_district_range = range(self.k_minority, self.k)
					return 'Infeasible', optimization_time

			elif self.m.status == 9:

				model_objective = 'Timeout'
				self.k_minority = k_minority
				self.minority_district_range = range(self.k_minority)
				self.majority_district_range = range(self.k_minority, self.k)
				return model_objective, optimization_time

	def black_population_map(self):
		
		#colors_list = ['w', 'lightcyan','cornflowerblue','lightcoral', 'red', 'firebrick']
		#colors_list = ['w', 'lightcoral','indianred', 'red', 'firebrick', 'darkred']
		colors_list = ['white', 'whitesmoke', 'gainsboro', 'silver', 'grey']
		custom_cmap = matplotlib.colors.ListedColormap(colors_list)

		

		self.shape_info['black_population_plot'] = -1
		color = 0
		for v in self.G.nodes():
			print(self.parcel_data['nodes'][v]['BASENAME'])
			
			pop_percent = self.minority_population[v] / self.voting_age_population[v]
			print(pop_percent)

			if pop_percent < 0.2:
				color = '0% - 20%'

			elif pop_percent < .3:
				color = '20% - 30%'
				
			elif pop_percent < .4:
				color = '30% - 40%'
				
			elif pop_percent < .5:
				color = '40% - 50%'
				
			else:
				color = '50% - 100%'
					
			self.shape_info['black_population_plot'][v] = color

		fig = plt.figure(1, figsize=(9,15), dpi = 80) 
		ax = fig.add_subplot()
		
		oplot = self.shape_info.plot(ax=ax,column='black_population_plot', categorical=True, cmap = custom_cmap, edgecolor='black', legend_kwds={'loc':(0,.05)}, categories=['0% - 20%', '20% - 30%','30% - 40%', '40% - 50%', '50% - 100%'], legend=True, figsize = (10,12)).get_figure()
		ax.get_legend().set_title('% BVAP population of each county in Mississippi')
		
		for v in self.G.nodes():
			name = self.parcel_data['nodes'][v]['BASENAME']
			x = float(self.parcel_data['nodes'][v]['INTPTLON']) -.155
			y = float(self.parcel_data['nodes'][v]['INTPTLAT']) -.05

			if name == 'Lee':
				x += .05
				#y += .05
			elif name == 'Washington':
				x += .015
				y += -.05
			elif name == 'Tunica':
				x += .02
				y += -.01
			elif name == 'Benton':
				x += .04
				y += .1
			elif name == 'Sunflower':
				x += -.01
			elif name == 'Leflore':
				x += .02
			elif name == 'Montgomery':
				x += .01
				y += -.03
			elif name == 'Washington':
				x += .01
			elif name == 'Humphreys':
				x += .02
				y += .05
			elif name == 'Sharkey':
				x += .065
				y += .05
			elif name == 'Warren':
				x += .14
				y += .07
			elif name == 'Smith':
				x += .03
			elif name == 'Lawrence':
				x += .02
			elif name == 'Pike':
				x += .02
			elif name == 'Walthall':
				x += .02
			elif name == 'Lamar':
				x += .02
				#y += .05
			elif name == 'Forrest':
				x += .02
				y += .23
			elif name == 'Perry':
				x += .02
			elif name == 'Tishomingo':
				x += .02
				y += .07
			elif name == 'Issaquena':
				x += .02
				y += .05
			elif name == 'Dolliver':
				x += .3
				
			ax.text(x, y, name, fontsize=7, color='black')
			
		plt.axis('off')
		plt.show()
		plt.close()

	def upper_bound_s_model(self, v):
		
		DG = nx.to_directed(self.G)
		self.m = gp.Model()
		#self.m.Params.OutputFlag = 0
		self.m.Params.TIME_LIMIT = 3550

		self.m._t = self.m.addVars(self.G.nodes(), ub=1, name='t')
		self.m._g = self.m.addVars(DG.edges(), vtype=gp.GRB.BINARY, name='g')


		verticies_except_v = list(set(self.G.nodes()) - {v})
		self.m.setObjective(gp.quicksum(self.m._t[u] for u in verticies_except_v), gp.GRB.MAXIMIZE)

		#3b
		#self.m.addConstr(self.m._t[v] == 1)
		self.m._t[v].lb = 1
		#3c
		self.m.addConstr(gp.quicksum(self.m._g[v,j] for j in DG.neighbors(v)) == 1)
		#3d
		self.m.addConstr(gp.quicksum(self.m._g[j,v] for j in DG.neighbors(v)) == 0)
		#3?
		self.m.addConstrs(gp.quicksum(self.m._g[j,u] for j in DG.neighbors(u)) == self.m._t[u] for u in verticies_except_v)
		#3e1
		self.m.addConstrs(gp.quicksum(self.m._g[u,j] for j in DG.neighbors(u)) <= self.m._t[u] for u in verticies_except_v)
		#3e2
		#self.m.addConstrs(self.m._t[u] <= 1 for u in verticies_except_v)
		#3f1
		self.m.addConstr(self.L <= gp.quicksum(self.population[u] * self.m._t[u] for u in verticies_except_v))
		#3f2
		self.m.addConstr(self.U >= gp.quicksum(self.population[u] * self.m._t[u] for u in verticies_except_v))
		#####
		minority = False
		if minority:
			self.m.addConstr(gp.quicksum(self.minority_population[v] * self.m._t[u] for u in self.G.nodes()) >= self.f * gp.quicksum(self.voting_age_population[u] * self.m._t[u] for u in self.G.nodes()))

		
		self.m.addConstrs(self.m._t[i] + self.m._t[j] <= 1 + self.m._g[i,j] + self.m._g[j,i] for (i,j) in self.G.edges())
		self.m._DG = DG
		self.m._v = v

		self.m.Params.lazyConstraints = 1
		
	def upper_bound_s_algo(self):
		plot = False
		upper_bound_s = nx.diameter(self.G)
		time1 = time.time()
		number_of_iterations = 0
		for v in self.G.nodes():
			self.upper_bound_s_model(v)
			self.m.addConstr(gp.quicksum(self.m._t[u] for u in list(set(self.G.nodes()) - {v})) >= upper_bound_s + 1)
			self.m.optimize(upper_bound_s_callback)
			if self.m.status == 2:
				print("FEASIBLE")
				objective_value = self.m.getObjective().getValue()
				if objective_value > upper_bound_s:
					upper_bound_s = objective_value
					
					if plot:
						self.minority_districts = [[]]
						edge_list = []
						for i in self.m.getVars():
						
							if i.VarName[0] == 't':
								if i.X > 0.5:
									
									self.minority_districts[0].append(int(i.VarName[2:-1]))

							if i.VarName[0] == 'g':
									if i.X > 0.5:
										
										edge_list.append(tuple(map(int, i.VarName[2:-1].split(','))))						
						
			number_of_iterations += 1

			if time.time() - time1 > 3600:
				break

		time2 = time.time()
		upper_bound_s_time = time2 - time1

		#self.plotter('minority', 'AL bad district')
		if plot:
			self.shape_info['minority'] = -1
			color = 0
			for district in self.minority_districts:
				color += 1
				for v in district:
					self.shape_info['minority'][v] = color

			#oplot.grid = off
			oplot = self.shape_info.plot(column='minority', figsize=(10, 12), legend=True).get_figure()
			#pos = nx.spring_layout(G, seed=1)  # positions for all nodes
			#nx.draw(G, pos, with_labels=True)
			solution_subgraph = self.G.edge_subgraph(edge_list)
			#nx.draw_networkx_edges(solution_subgraph, width=5.0, alpha=.9, edge_color="tab:red")
			#nx.draw(self.G, pos=nx.get_node_attributes(self.G, 'pos'), node_size = ns, width = w, node_color=c,edge_color=ec)
			centroids = self.shape_info.centroid
			c_x = centroids.x
			c_y = centroids.y
			pos = {node:(c_x[node],c_y[node]) for node in self.G.nodes}
			for node in self.G.nodes():
				self.G.nodes[node]["pos"] = np.array(pos[node])
			nx.draw(solution_subgraph, pos=nx.get_node_attributes(self.G, 'pos'), node_size = 0, edge_color='red')
			#plt.show()
			oplot.savefig('Bad_'+ self.state + '_' + self.parcel_level + '_Map.png')
			plt.close()

		return upper_bound_s_time, upper_bound_s, number_of_iterations

	def create_upper_bound_s_table(self):
	
		if not os.path.exists('upper_bound_s_table.csv'):
			with open('upper_bound_s_table.csv', 'w') as doc:
				doc.write('State, parcel, k, m, n, upper bound s time, upper bound s, number of iterations, diameter \n')

		upper_bound_s_time, upper_bound_s, number_of_iterations = self.upper_bound_s_algo()
		
		
		with open('results/upper_bound_s/upper_bound_s_' + self.state + '_' + self.parcel_level + '.txt', 'w') as file:
			file.write(str(upper_bound_s))

		with open('upper_bound_s_table.csv', 'a') as doc:	
			print("Upper bound on s computed as: ", upper_bound_s)
			string = self.state + ',' + self.parcel_level + ',' + str(self.k) + ',' + str(len(self.G.edges())) + ',' + str(len(self.G.nodes())) + ',' + '{0:.2f}'.format(upper_bound_s_time) + ',' + str(upper_bound_s) + ',' + str(number_of_iterations) + ',' + str(nx.diameter(self.G)) + '\n'
			
			doc.write(string)

	def DFJ_spanning_tree(self, longest_path=False, biggest_district=False, starting_tree=False, gerrychain_s_val=False):
		m = gp.Model()
		m.Params.lazyConstraints = 1
		m.Params.TIME_LIMIT = 3600

		centroids = self.shape_info.centroid
		c_x = centroids.x
		c_y = centroids.y
		pos = {node:(c_x[node],c_y[node]) for node in self.G.nodes}
		for node in self.G.nodes():
			self.G.nodes[node]["pos"] = np.array(pos[node])
		
		DG = nx.DiGraph(self.G)

		
		if longest_path:
			source = longest_path[0]
			sink = longest_path[-1]
			r_q = source
		else:
			r_q = 0

		#

		#4a
		m._q = m.addVars(DG.edges, vtype=GRB.BINARY, name='q' )
		# Constraints: each node (besides r) should have one incoming arc
		m.addConstrs(gp.quicksum(m._q[j,i] for j in DG.neighbors(i)) == 1 for i in self.G.nodes if i != r_q)
		m.addConstr(gp.quicksum(m._q[j,r_q] for j in DG.neighbors(r_q)) == 0)

		m._DG = DG
		m._G = self.G
		m._r_q = r_q

		if starting_tree:
			for edge in starting_tree:
				m._q[edge].start = 1

		#4b-4e
		
		m._b = m.addVars(DG.edges, vtype=GRB.BINARY, name='b')
		m._r = m.addVars(DG.nodes, vtype=GRB.BINARY, name='r')

		m.addConstrs(m._b[i,j] + m._b[j,i] <= m._q[i,j] + m._q[j,i] for (i,j) in self.G.edges())
		m.addConstrs(gp.quicksum(m._b[u,v] for u in DG.neighbors(v)) + m._r[v] <= 1 for v in DG.nodes)
		m.addConstr(gp.quicksum(m._r[v] for v in DG.nodes) == 1)
		
		#POPULATION FLOW 4f-4i
		
		m._f = m.addVars(DG.edges, name='f', lb=0)
		m._g = m.addVars(DG.nodes, name='g')

		m.addConstrs(self.L * m._r[v] <= m._g[v] for v in DG.nodes)
		m.addConstrs(self.U * m._r[v] >= m._g[v] for v in DG.nodes)

		m.addConstrs(m._g[v] + gp.quicksum(m._f[u,v] for u in DG.neighbors(v)) == gp.quicksum(m._f[v,u] for u in DG.neighbors(v)) + self.population[v] * (gp.quicksum(m._b[u,v] for u in DG.neighbors(v)) + m._r[v]) for v in DG.nodes)

		m.addConstrs(m._f[u,v] <= (self.U - self.population[u]) * m._b[u,v] for (u,v) in DG.edges)
		
		#LONGEST SHORTEST PATH 4j-4p

		
		m._t = m.addVars(DG.nodes, vtype=GRB.BINARY, name='t')
		m._c = m.addVars(DG.nodes, vtype=GRB.BINARY, name='c', ub=1)
		m._h = m.addVars(DG.edges, vtype=GRB.BINARY, name='h', lb=0)

		m.addConstr(gp.quicksum(m._t[v] for v in DG.nodes) == 1)
		
		m.addConstrs(gp.quicksum(m._h[u,v] for u in DG.neighbors(v)) + m._r[v] == m._c[v] for v in DG.nodes)

		m.addConstrs(m._r[v] + m._t[v] <= m._c[v] for v in DG.nodes)

		m.addConstrs(gp.quicksum(m._h[v,u] for u in DG.neighbors(v)) - gp.quicksum(m._h[u,v] for u in DG.neighbors(v)) == m._r[v] - m._t[v] for v in DG.nodes)

		m.addConstrs(m._h[u,v] <= m._b[u,v] for (u,v) in DG.edges)

		m.addConstrs(m._c[u] + m._c[v] <= 1 + m._h[u,v] + m._h[v,u] for (u,v) in self.G.edges)
		

		if longest_path:
			for v in longest_path:
				print(longest_path)
				m._c[v].start = 1
			#m._c[source].start = 1
			#m._c[sink].start = 1
			m._t[sink].start = 1


		m.setObjective(gp.quicksum(m._c[v] for v in DG.nodes) - 1, GRB.MAXIMIZE)
	
		
		time1 = time.time()
		m.optimize(DFJ_callback)
		time2 = time.time()




		plot = False
		if plot: 
			tree_edges_q = [(i,j) for i,j in DG.edges if m._q[i,j].x > 0.5]	
			tree_edges_b = [(i,j) for i,j in DG.edges if m._b[i,j].x > 0.5]	
			tree_edges_f = [(i,j) for i,j in DG.edges if m._f[i,j].x > 0.5]	
			tree_edges_h = [(i,j) for i,j in DG.edges if m._h[i,j].x > 0.5]

			district = [j for i,j in DG.edges if m._f[i,j].x >0.5]
			source = [i for i in DG.nodes if m._r[i].x > 0.5]
			terminal = [i for i in DG.nodes if m._t[i].x > 0.5]

			district.append(source[0])
			pop_counter = 0
			for parcel in district:
				
				pop_counter += self.population[parcel]
			
			
			node_color_list = []
			for node in DG.nodes:
				if node == r_q:
					node_color_list.append('yellow')
				elif node in terminal:
					node_color_list.append('blue')
				elif node in source:
					node_color_list.append('red')
				elif node in district:
					node_color_list.append('green')
				else:
					node_color_list.append('white')
			#node_color_list = ["green" if i in district  else "black" for i in DG.nodes]
			

			plotting_options_dict = {"node_size":100, 'edge_width':3, 'node_color':node_color_list, 'edge_color':'red'}
			
			oplot = self.shape_info.plot(figsize=(10, 12), legend=True,).get_figure()
			tree_graph = nx.DiGraph()
			tree_graph.add_nodes_from(DG)
			tree_graph.add_edges_from(tree_edges_q)

			nx.draw(tree_graph, pos=nx.get_node_attributes(self.G, 'pos'), node_size=plotting_options_dict['node_size'], width=plotting_options_dict['edge_width'], node_color=plotting_options_dict['node_color'],edge_color='k', with_labels=True)
								

			oplot = self.shape_info.plot(figsize=(10, 12), legend=True,).get_figure()
			tree_graph = nx.DiGraph()
			tree_graph.add_nodes_from(DG)
			tree_graph.add_edges_from(tree_edges_f)
			
			nx.draw(tree_graph, pos=nx.get_node_attributes(self.G, 'pos'), node_size=plotting_options_dict['node_size'], width=plotting_options_dict['edge_width'], node_color=plotting_options_dict['node_color'],edge_color='b', with_labels=True)
			
			tree_graph = nx.DiGraph()
			tree_graph.add_nodes_from(DG)
			tree_graph.add_edges_from(tree_edges_b)
			
			nx.draw(tree_graph, pos=nx.get_node_attributes(self.G, 'pos'), node_size=plotting_options_dict['node_size'], width=plotting_options_dict['edge_width'], node_color=plotting_options_dict['node_color'],edge_color=plotting_options_dict['edge_color'])

			tree_graph = nx.DiGraph()
			tree_graph.add_nodes_from(DG)
			tree_graph.add_edges_from(tree_edges_h)
			
			nx.draw(tree_graph, pos=nx.get_node_attributes(self.G, 'pos'), node_size=plotting_options_dict['node_size'], width=plotting_options_dict['edge_width'], node_color=plotting_options_dict['node_color'],edge_color='blue')
			

			plt.show()
		

		if not os.path.exists('./upper_bound_s_table.csv'):
			with open('upper_bound_s_table.csv', 'w') as doc:
				doc.write('State, Parcel Level, k, n, m, LB, UB, comp time, diameter, gerrychain_value \n')

		with open('upper_bound_s_table.csv', 'a') as doc:
			doc.write(self.state + ',' + self.parcel_level + ',' + str(self.k) + ',' + str(len(self.G.nodes)) + ',' + str(len(self.G.edges)) + ',' + str(m.ObjVal) + ',' + str(m.ObjBound) + ',' + '{0:.2f}'.format(time2-time1) + ',' + str(nx.diameter(self.G)) + ',' + str(gerrychain_s_val) + '\n')

	def get_solution_spanning_tree_with_gerrychain(self):
		iterations = 50
		s = -1
		my_updaters = {'population': updaters.Tally('P0010001', alias='population')}
		start = recursive_tree_part(self.G, range(self.k),sum(self.G.nodes[i]['P0010001'] for i in self.G.nodes())/self.k,'P0010001', self.population_deviation/2,1)
		initial_partition = GeographicPartition(self.G, start, updaters = my_updaters)

		proposal = partial(recom,
						pop_col='P0010001',
						pop_target=sum(self.G.nodes[i]['P0010001'] for i in self.G.nodes())/self.k,
						epsilon=self.population_deviation/2,
						node_repeats=2
						)

		compactness_bound = constraints.UpperBound(
			lambda p: len(p['cut_edges']),
			1000 * len(initial_partition['cut_edges'])
		)

		pop_constraint = constraints.within_percent_of_ideal_population(initial_partition, self.population_deviation/2)

		my_chain = MarkovChain(
			proposal=proposal,
			constraints=[
				pop_constraint,
				compactness_bound
			],
			accept=accept.always_accept,
			initial_state=initial_partition,
			total_steps=iterations
		)


		max_s = 0
		print('In GerryChain heuristic, current s ', end='')

		all_maps = []
		pareto_frontier = []
		obj_vals = []
		for partition in my_chain:
			number_minority_district = 0
			current_s = -1
			for district in range(self.k):
				
				list_of_nodes_in_current_district = []
				for node in partition.graph:
						if partition.assignment[node] == district:
							list_of_nodes_in_current_district.append(node)

				district_subgraph = self.G.subgraph(list_of_nodes_in_current_district)
				district_diamater = nx.diameter(district_subgraph)

				if district_diamater > current_s:
					current_s = district_diamater

				if current_s > max_s:
					max_s = current_s
					biggest_district = district_subgraph


			print((current_s),',',sep='',end=' ')
			obj_vals.append([current_s])
			all_maps.append([partition, current_s])

		print('Best heuristic solution has the district diameter =', max_s)

		all_shortest_paths = nx.all_pairs_shortest_path(biggest_district)
		longest_path = []
		

		for i in all_shortest_paths:
			vertex = i[0]
			dict_from_i = i[1]
			longest_path_from_current_vertex = max(dict_from_i.values(), key=len)
			if len(longest_path_from_current_vertex) > len(longest_path):
				longest_path = longest_path_from_current_vertex

		source = longest_path[0]
		sink = longest_path[-1]
		
		plot = False
		if plot:
			self.shape_info['district'] = -1

			for parcel in biggest_district:
				self.shape_info['district'][parcel] = 0
			for parcel in longest_path:
				self.shape_info['district'][parcel] = 1
			self.shape_info['district'][source] = 2
			self.shape_info['district'][sink] = 3

			oplot = self.shape_info.plot(column='district', figsize=(10, 12), legend=True).get_figure()
			plt.show()

		rest_of_path = longest_path[1:-1]

		return longest_path, biggest_district, max_s

	def spanning_tree_generator(self, path, district):
		model_one = gp.Model()
		model_one.Params.lazyConstraints = 1
		model_one.Params.TIME_LIMIT = 3600

		
		centroids = self.shape_info.centroid
		c_x = centroids.x
		c_y = centroids.y
		pos = {node:(c_x[node],c_y[node]) for node in self.G.nodes}
		for node in self.G.nodes():
			self.G.nodes[node]["pos"] = np.array(pos[node])
		
		
		
		DG = nx.DiGraph(district)
	

		#4a
		model_one._q = model_one.addVars(DG.edges, vtype=GRB.BINARY, name='q' )
		r_q = path[0]

		
		last_vertex = r_q

		for i in path[1:]:
			model_one._q[(last_vertex, i)].lb = 1
			last_vertex = i
		# Constraints: each node (besides r) should have one incoming arc
		model_one.addConstrs(gp.quicksum(model_one._q[j,i] for j in DG.neighbors(i)) == 1 for i in district.nodes if i != r_q)
		model_one.addConstr(gp.quicksum(model_one._q[j,r_q] for j in DG.neighbors(r_q)) == 0)

		model_one._DG = DG
		model_one._G = district
		model_one._r_q = r_q

		model_one.optimize(DFJ_callback)

		edges_to_fix = []
		for v in model_one.getVars():
			if v.X > 0.5:
				#print(f"{v.VarName} = {v.X}")
				vertex_one = int(v.VarName.split(',')[0][2:])
				vertex_two = int(v.VarName.split(',')[1][:-1])
		
				edges_to_fix.append([vertex_one, vertex_two])

		#2nd model
		model_two = gp.Model()
		model_two.Params.lazyConstraints = 1
		model_two.Params.TIME_LIMIT = 3600

		DG = nx.DiGraph(self.G)
		

		#4a
		model_two._q = model_two.addVars(DG.edges, vtype=GRB.BINARY, name='q' )
		r_q = path[0]
		

		for edge in edges_to_fix:
			model_two._q[edge[0], edge[1]].lb = 1
		# Constraints: each node (besides r) should have one incoming arc
		model_two.addConstrs(gp.quicksum(model_two._q[j,i] for j in DG.neighbors(i)) == 1 for i in DG.nodes if i != r_q)
		model_two.addConstr(gp.quicksum(model_two._q[j,r_q] for j in DG.neighbors(r_q)) == 0)

		model_two._DG = DG
		model_two._G = self.G
		model_two._r_q = r_q

		model_two.optimize(DFJ_callback)
		tree_edges_q = [(i,j) for i,j in DG.edges if model_two._q[i,j].x > 0.5]	
		
		plot = False
		if plot: 
			
			
			node_color_list = ["green" if i in district  else "black" for i in DG.nodes]
			

			plotting_options_dict = {"node_size":20, 'edge_width':3, 'edge_color':'red'}
			
			oplot = self.shape_info.plot(figsize=(10, 12), legend=True,).get_figure()
			tree_graph = nx.DiGraph()
			tree_graph.add_nodes_from(DG)
			tree_graph.add_edges_from(tree_edges_q)

			#nx.draw(tree_graph, pos=nx.get_node_attributes(self.G, 'pos'), node_size=plotting_options_dict['node_size'], width=plotting_options_dict['edge_width'],edge_color='r', with_labels=False)
			nx.draw(tree_graph, pos=nx.get_node_attributes(self.G, 'pos'), node_size=plotting_options_dict['node_size'], width=plotting_options_dict['edge_width'], node_color=node_color_list,edge_color='k', with_labels=False)
			
			plt.show()

		return tree_edges_q

	def generate_symmetry_table(self):

		if not os.path.exists('symmetry_table.csv'):
			with open('symmetry_table.csv', 'w') as doc:
				doc.write('State, parcel level, k, n, m, s, obj_wo_symmetry, time_wo_symmetry, obj_w_symmetry, time_w_symmetry, \n')

		string = self.state + ',' + self.parcel_level + ',' + str(self.k) + ',' + str(len(self.G.nodes)) + ',' + str(len(self.G.edges)) + ','  + str(self.s) + ','

		time1 = time.time()
		self.labeling_model(True, True, False)
		time2 = time.time()

		
		no_symmetry_time = format(time2 - time1, '.2f')
		no_symmetry_obj = str(self.num_minority_districts)
		

		string += no_symmetry_obj + ',' + no_symmetry_time + ','  
		
		
		time1 = time.time()
		self.labeling_model(True, True, True)
		time2 = time.time()

		symmetry_time = format(time2 - time1, '.2f')
		symmetry_obj = str(self.num_minority_districts)

		string += symmetry_obj + ',' + symmetry_time + ', \n'

	

		with open('symmetry_table.csv', 'a') as doc:
			doc.write(string)

	def check_points_with_prescribed_model(self):

		if not os.path.exists('point_check.csv'):
			with open('point_check.csv', 'w') as doc:
				doc.write('State, parcel, n, m, num_minority_districts, s, feasible/infeasible/unknown, time \n')

		

		s_save = self.s
		
		interesting_k = range(self.k_minority, 0, -1)
		
		
		for k in interesting_k:
			string = self.state + ',' + self.parcel_level + ',' + str(len(self.G.nodes())) + ',' + str(len(self.G.edges())) + ',' + str(k) + ','
			try_higher_s = True
			self.s = s_save - 1
			self.k_minority = k
			time_out_in_a_row = 0
			while try_higher_s:

				self.s += 1
				string += str(self.s) + ','

				self.minority_district_range = range(self.k_minority)
				self.majority_district_range = range(self.k_minority, self.k)
			
				self.make_sure_preprocess_exists()
				print("MODEL BEING SOLVED WITH:")
				print("S:", str(self.s))
				print("k_minority:", str(self.k_minority))
				self.prescribed_labeling_model()
					
				if self.m.status == 2:
					#optimal
					string += 'Feasible' + ',' + '{0:.2f}'.format(self.model_runtime) + '\n'
					try_higher_s = False

				if self.m.status == 3:
					#infeasible
					string += 'Infeasible' + ',' + '{0:.2f}'.format(self.model_runtime) + '\n' + ',' + ','  + ',' + ',' + str(k) + ','
					


				if self.m.status == 9:
					#timeout
					

					time_out_in_a_row += 1
					if time_out_in_a_row == 2:
						try_higher_s = False
						string += 'Unknown' + ',' + 'timeout \n'
					else:
						string += 'Unknown' + ',' + 'timeout \n' + ',' + ','  + ',' + ',' + str(k) + ','

		


			with open('point_check.csv', 'a') as doc:
				doc.write(string)

	def short_bursts(self):
		
		minority = 'Black'
		deviation = 0.01      # 0.01 means 1% means +/-0.5%
		total_steps = 10000 # SHORT BURST STEPS
		total_time = 3600
		#total_steps = 100
		burst_length = 10
		num_bursts = round( total_steps / burst_length )

		random.seed(2024)

		ideal_population = self.M
		DG = nx.DiGraph(self.G)
		DG._k = self.k
		DG._L = self.L
		DG._U = self.U
		for v in DG.nodes():
			DG.nodes[v]['TOTPOP'] = self.population[v]
		time1 =time.time()
		mip_districts = recursive_bipartition_heuristic(DG)
		time2 = time.time()
		mip_labeling = { i : j for j in range(self.k) for i in mip_districts[j] }

		warm_start_s = -1
		for district in mip_districts:
			district_subgraph = self.G.subgraph(district)
			diameter_of_district = nx.diameter(district_subgraph)
			if diameter_of_district > warm_start_s:
				warm_start_s = diameter_of_district

		
		if not os.path.exists('short_bursts_first_half.csv'):
			with open('short_bursts_first_half.csv', 'w') as doc:
				doc.write('State, parcel, n, m, num_minority_districts, s, iterations, region_aware, total time \n')

		print("seed_plan_"+self.state+"_"+self.parcel_level+" =",mip_districts)
		report_metrics(self.G, mip_districts, minority)
		
		with open('short_bursts_first_half.csv', 'a') as doc:
			num_minority = number_of_gingles_districts(self.G, mip_districts, minority)
			doc.write(self.state + ',' + self.parcel_level + ',' + str(len(self.G.nodes)) + ',' + str(len(self.G.edges)) + ',' + str(num_minority) + ',' + str(warm_start_s) + ',' +  'NA' + ',' + 'Warm Start' + ',' + '{0:.2f}'.format(time2-time1) + '\n')
		
		

		for region_aware in {False}:
			time1 = time.time()
			print("\n*****************************************")
			################print("Starting short bursts for",self.state,self.parcel_level,level)
			print("Region-aware:",region_aware)
			print("*****************************************\n")
			
			# GerryChain/ShortBursts/Gingleator setup
			chain_updaters = {
			"population": updaters.Tally("TOTPOP", alias="population"),
			"VAP": updaters.Tally("VAP"),
			"MVAP": updaters.Tally("MVAP")
			}
			
			initial_partition = Partition(graph=self.G,
				assignment=mip_labeling,
			updaters=chain_updaters
			)

			my_surcharge = 0.5 if region_aware else 0.0
			proposal = partial(
				proposals.recom,
				pop_col="TOTPOP",
				pop_target=ideal_population,
				epsilon=deviation/2,
				node_repeats=1,
				region_surcharge={"COUNTY": my_surcharge}
			)
			
			constraints = constraints_class.within_percent_of_ideal_population(initial_partition, deviation/2)
				
			gingles = Gingleator(
				proposal,
				constraints,
				initial_partition,
				minority_pop_col="MVAP",
				total_pop_col="VAP",
				score_function=Gingleator.reward_partial_dist
			)
	
			# run short bursts
			max_scores_sb = np.zeros(total_steps)
			scores_sb = np.zeros(total_steps)
			pp = np.zeros(total_steps)
			s1 = np.zeros(total_steps)
			s2 = np.zeros(total_steps)
			s3 = np.zeros(total_steps)
			
			incumbent_plan = mip_districts.copy()
			incumbent_s = -1
			incumbent_i = -1
			for i, part in enumerate(gingles.short_bursts(burst_length, num_bursts, with_progress_bar=True)):
				max_scores_sb[i] = gingles.best_score
				scores_sb[i] = gingles.score(part)
				districts = [ list() for j in range(self.k) ]
				for v in self.G.nodes:
					j = part.assignment[v]
					districts[j].append(v)
				s1[i] = number_of_counties_split(self.G, districts)
				s2[i] = number_of_county_splits(self.G, districts)
				current_s = -1
				for district in districts:
					district_subgraph = self.G.subgraph(district)
					diameter_of_district = nx.diameter(district_subgraph)
					if diameter_of_district > current_s:
						current_s = diameter_of_district
				s3[i] = current_s
				pp[i] = average_polsby_popper(self.G, districts)
				
				num_minority = number_of_gingles_districts(self.G, districts, minority)
				with open(self.state + '_SAM RECORD.txt', 'a') as doc:
					doc.write(str(s3[i]) + ', ' + str(num_minority) + '\n')
				# update incumbent?
				
				case1 = math.floor( scores_sb[i] ) > math.floor( scores_sb[incumbent_i] )
				case2 = math.floor( scores_sb[i] ) == math.floor( scores_sb[incumbent_i] )
				#case3 = pp[i] > pp[incumbent_i] 
				case3 = s3[i] < s3[incumbent_i] 
				if case1 or (case2 and case3):
					incumbent_i = i
					incumbent_s = current_s
					incumbent_plan = districts.copy()
				current_time = time.time()
				
				if current_time - time1 > total_time:
					num_iterations = i
					#break
			# reporting
			print("i gingles pp s1 s2 s3")
			for i in range(total_steps):
				if i%1000==0:
					print(i,round(scores_sb[i], 4), round(pp[i],4), round(s1[i]), round(s2[i]), s3[i])

			self.districts = incumbent_plan
			self.minority_districts = []
			for district in incumbent_plan:
				total_population = 0
				minority_population = 0
				for v in district:
					minority_population += self.minority_population[v]
					total_population += self.voting_age_population[v]		
				if minority_population >= self.f * total_population:
					self.minority_districts.append(district)
				print(district)
				print(minority_population)
				print(total_population)

			print(self.minority_districts)

			
			print("incumbent_plan_"+self.state+"_"+self.parcel_level+" =",incumbent_plan)
	
			report_metrics(self.G, incumbent_plan, minority)

			time2 = time.time()
			


			
			with open('short_bursts_first_half.csv', 'a') as doc:
	
				num_minority = number_of_gingles_districts(self.G, incumbent_plan, minority)
				doc.write(self.state + ',' + self.parcel_level + ',' + str(len(self.G.nodes)) + ',' + str(len(self.G.edges)) + ',' + str(num_minority) + ',' + str(incumbent_s) + ',' +  str(i) + ',' + str(region_aware) + ',' + '{0:.2f}'.format(time2-time1) + '\n')
			

	def lower_bound_s_minority_dist_only(self):

		original_s = self.s
		if not os.path.exists('lower_bound_s_new_idea.csv'):
			with open('lower_bound_s_new_idea.csv', 'w') as doc:
				doc.write('State, parcel, n, m, num_minority_districts, s, model status, total time \n')

		loop = True
		while loop:
			
			self.prescribed_labeling_model(include_majority_variables=False)

			if self.m.status == 2:
				model_result = 'Feasible'

			elif self.m.status == 3:
				model_result = 'Infeasible'
				

			elif self.m.status == 9:
				model_result = 'Timeout'

			else:
				model_result = 'unanticpated termination gurobicode: ' + str(self.m.status)



			with open('lower_bound_s_new_idea.csv', 'a') as doc:

				doc.write(self.state + ',' + self.parcel_level + ',' + str(len(self.G.nodes)) + ',' + str(len(self.G.edges)) + ',' + str(self.k_minority) + ',' + str(self.s) + ',' + model_result + ',' + '{0:.2f}'.format(self.model_runtime) + '\n')

			if model_result == 'Infeasible':
				self.s += 1
			elif model_result == 'Feasible':
				self.s = original_s
				self.k_minority -= 1
				self.minority_district_range = range(self.k_minority)
				self.majority_district_range = range(self.k_minority, self.k)

				if self.k_minority == 0:
					loop = False				
			else:
				loop = False

		


if __name__ == '__main__':

	state_codes = get_state_codes()
	congressional_districts = get_congressional_codes()

	if os.path.exists('results.csv'):
		os.remove('results.csv')

	file = open('data.json', 'r')
	data = json.load(file)
	

	for request in data['first_half']:
		print('State: ', request['state'])
		print('Parcel level: ', request['parcel_level'])

		instance = problem_instance(request['state'], request['parcel_level'], request['s'])

 
		if instance.error_message != '':
			print('ERROR: ', instance.error_message)
			break
		
	

		#instance.black_population_map()
		#instance.draw_symmetry_map()

		#instance.generate_gerrychain_run(10)
		#instance.make_sure_preprocess_exists()
		#instance.generate_instance_info()
		#instance.generate_upper_bound_minority_table()
		#instance.generate_lower_bound_s_table()
		#instance.find_interesting_s_vals()
		#instance.generate_fixing_info()
		
		#instance.check_points_with_prescribed_model()
		
		#instance.short_bursts()
		#instance.lower_bound_s_minority_dist_only()
		
		#instance.draw_fixing_map()
		#instance.black_population_map()

	file.close()