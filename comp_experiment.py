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
try:
	from gerrychain import MarkovChain
	from gerrychain import (GeographicPartition, Graph, updaters, constraints, accept, proposals)
	from gerrychain import Partition
	from gerrychain import constraints as constraints_class
	from gerrychain.optimization import Gingleator
	from gerrychain.tree import recursive_tree_part
	from gerrychain.proposals import recom
	_HAS_GERRYCHAIN = True
except ImportError:
	# gerrychain is only used by the GerryChain heuristic methods and for graph
	# loading. The MIP formulations and the A1 benchmark do not need it.
	_HAS_GERRYCHAIN = False
	Graph = None
from functools import partial
from matplotlib.colors import LinearSegmentedColormap
from aux_functions import *


def load_graph_json(path):
	"""Load a contiguity graph from a JSON adjacency file as a networkx graph.
	Uses gerrychain's Graph.from_json when available (original behavior, so the
	GerryChain methods keep working); otherwise falls back to plain networkx so
	the MIP formulations and the A1 benchmark run without gerrychain installed."""
	if _HAS_GERRYCHAIN:
		return Graph.from_json(path)
	with open(path) as f:
		data = json.load(f)
	return nx.adjacency_graph(data, multigraph=False)


# Module-level lookups used by problem_instance. These were previously defined
# only inside the __main__ block, so importing this module (e.g., from the A1
# benchmark driver) left them undefined; define them here so import works.
state_codes = get_state_codes()
congressional_districts = get_congressional_codes()
pandas.options.mode.chained_assignment = None


class problem_instance(object):

	def __init__(self, state, parcel_level, s, group):

		self.state = state
		self.parcel_level = parcel_level
		self.model_time_limit = 3600
		self.group = group

		self.error_message = ''
		# Set up-front so that a labeling_model(..., parcel_fixing=False, ...) run
		# does not dereference an unset attribute in add_lazy_distance_cuts.
		self.minority_fixings = []
		self.f = 0.5
		
		self.population_deviation = .01

		if parcel_level == 'county':
			parcel_level_plural = 'counties'

		if parcel_level == 'tract':
			parcel_level_plural = 'tracts'

		self.k = congressional_districts[state.upper()]

		print('Attempting to read in json data')
		if os.path.exists('./raw_data/' + parcel_level + '/json/' + state + '_' + parcel_level_plural + '.json'):
			self.G = load_graph_json('./raw_data/' + parcel_level + '/json/' + state + '_' + parcel_level_plural + '.json')

			f = open('./raw_data/' + parcel_level + '/json/' + state + '_' + parcel_level_plural + '.json')
			self.parcel_data = json.load(f)

			self.population = [self.parcel_data['nodes'][i]['P0010001'] for i in self.G.nodes()]
			for v in self.G.nodes():
				self.G.nodes[v]['TOTPOP'] = self.population[v]
			self.voting_age_population = [self.parcel_data['nodes'][i]['P0030001'] for i in self.G.nodes()]
			for v in self.G.nodes():
				self.G.nodes[v]['VAP'] = self.voting_age_population[v]
			


			if self.group == 'black':
				# ANY-PART Black VAP (alone OR in combination) -- the VRA/Gingles
				# standard, matching Belotti et al.'s carve definition so the
				# majority-minority comparison is apples-to-apples. (Previously
				# this used P0030004 = Black alone.)
				_black_anypart = [
					'P0030004', 'P0030011', 'P0030016', 'P0030017', 'P0030018',
					'P0030019', 'P0030027', 'P0030028', 'P0030029', 'P0030030',
					'P0030037', 'P0030038', 'P0030039', 'P0030040', 'P0030041',
					'P0030042', 'P0030048', 'P0030049', 'P0030050', 'P0030051',
					'P0030052', 'P0030053', 'P0030058', 'P0030059', 'P0030060',
					'P0030061', 'P0030064', 'P0030065', 'P0030066', 'P0030067',
					'P0030069', 'P0030071']
				self.minority_population = [
					sum(self.parcel_data['nodes'][i][c] for c in _black_anypart)
					for i in self.G.nodes()]
			elif self.group == 'hispanic':
				# Hispanic any-part == Hispanic VAP (single code), unchanged.
				self.minority_population = [self.parcel_data['nodes'][i]['P0040002'] for i in self.G.nodes()]


			for v in self.G.nodes():
				if self.group != 'none':
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
			if self.group != 'none':
					
				
				if os.path.exists('results_' + self.group + '/upper_bound_minority_districts/upper_bound_minority_districts_' + self.state + '_' + self.parcel_level + '.txt'):
					with open('results_' + self.group + '/upper_bound_minority_districts/upper_bound_minority_districts_' + self.state + '_' + self.parcel_level + '.txt', 'r') as file:
						for line in file:
							# int(line) fails when a previous run cached a float,
							# e.g. '33.0' from a time-limited solve that stored a
							# Gurobi bound rather than a proven integer optimum.
							# Round through float so an old cache cannot poison
							# every later step with
							#   TypeError: 'float' object cannot be interpreted
							#   as an integer     (from range(self.k_minority))
							self.k_minority = int(round(float(line.strip())))
							print('Upper bound on minority districts read in as :', str(self.k_minority))

				else:
					with open('results_' + self.group + '/upper_bound_minority_districts/upper_bound_minority_districts_' + self.state + '_' + self.parcel_level + '.txt', 'w') as doc:
						self.minority_district_range = range(self.k)
						self.k_minority = self.k
						self.find_upper_bound_minority_districts()
						self.k_minority = int(round(float(self.k_minority)))
						doc.write(str(self.k_minority))
						print('Upper bound on minority districts comupted and stored as :', str(self.k_minority))

				self.minority_district_range = range(self.k_minority)
				self.majority_district_range = range(self.k_minority, self.k)


			self.find_upper_bound_s()

			if os.path.exists('results_' + self.group + '/lower_bound_s/lower_bound_s_' + self.state + '_' + self.parcel_level + '.txt'):
				with open('results_' + self.group + '/lower_bound_s/lower_bound_s_' + self.state + '_' + self.parcel_level + '.txt', 'r') as file:
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
		
		# count majority-minority districts in the enacted plan: a district is MM
		# when its minority voting-age population is at least half of its total
		# voting-age population (same any-part VAP definition used everywhere else).
		current_number_minority_districts = 0
		for district in districts:
			minority_vap = 0
			total_vap = 0
			for parcel in district:
				minority_vap += self.minority_population[parcel]
				total_vap += self.voting_age_population[parcel]
			if total_vap > 0 and minority_vap >= 0.5 * total_vap:
				current_number_minority_districts += 1

		return current_s, split_number, current_number_minority_districts
	
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
		
			stability_number, max_independent_set = compute_stability_number('compute', self.G, self.state, self.k, s, self.group)
			

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
		with open('results_' + self.group + '/lower_bound_s/lower_bound_s_' + self.state + '_' + self.parcel_level + '.txt', 'w') as file:
			file.write(str(self.lower_bound_s))

	def find_upper_bound_s(self):

		if os.path.exists('results_' + self.group + '/upper_bound_s/upper_bound_s_' + self.state + '_' + self.parcel_level + '.txt'):
			with open('results_' + self.group + '/upper_bound_s/upper_bound_s_' + self.state + '_' + self.parcel_level + '.txt', 'r') as file:
				for line in file:
					self.upper_bound_s = int(float(line))
					print("Upper bound on s read as: ", self.upper_bound_s)
		else:
			self.upper_bound_s = nx.diameter(self.G)

			with open('results_' + self.group + '/upper_bound_s/upper_bound_s_' + self.state + '_' + self.parcel_level + '.txt', 'w') as file:
				file.write(str(self.upper_bound_s))

			print("Upper bound on s computed as: ", self.upper_bound_s)

	def find_upper_bound_minority_districts(self, continuity=False, relax=False):

		self.minority_district_ONLY_model(0, relax, continuity)
		self.k_minority = self.num_minority_districts
		print(self.k_minority)
		print(self.num_minority_districts)

	def find_interesting_s_vals(self):


		if not os.path.exists('results_' + self.group + '/upper_bound_minority_table.csv'):
			with open('results_' + self.group + '/upper_bound_minority_table.csv', 'w') as doc:
				doc.write('State, parcel, k, n, m, objective, time (s), x & z fixed (%)')


		output_string = 's_vaule, num_minority_districts, time_out \n'
		s = self.lower_bound_s - 1

		continue_search = True

		if not os.path.exists('results_' + self.group + '/interesting_s_vals.csv'):
			with open('results_' + self.group + '/interesting_s_vals.csv', 'a') as doc:
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

			with open('./results_' + self.group + '/interesting_s_vals.csv', 'a') as doc:
				current_string = self.state + ',' + current_string
				doc.write(current_string)
		with open('./results_' + self.group + '/interesting_s_vals/interesting_s_vals_' + self.state + '_' + self.parcel_level + '.csv', 'w') as doc:
			doc.write(output_string)

	def run_GerryChain_heuristic(self, iterations):

		# Imported lazily: gerrychain is only needed for the heuristic, so the other
		# tables keep working without it. (These names were previously undefined --
		# neither this module nor aux_functions imported them.)
		from gerrychain import GeographicPartition, MarkovChain, updaters, constraints, accept
		from gerrychain.tree import recursive_tree_part, bipartition_tree
		from gerrychain.proposals import recom
		import random as _random

		# --- seeding robustness -------------------------------------------------
		# recursive_tree_part draws a random spanning tree and looks for a
		# population-balanced edge to cut. At +/-0.5% deviation with many districts
		# this can fail, and gerrychain's default budget is only
		#     method = partial(bipartition_tree, max_attempts=10000)
		# with pair reselection DISABLED. When the budget runs out it raises
		#     RuntimeError: Could not find a possible cut after 10000 attempts
		# and the whole instance dies before the chain even starts. That is what
		# killed FL on the first sweep.
		#
		# NOTE this is a *seeding* budget, not the number of chain iterations. No
		# increase in `iterations` can rescue it, because the failure happens
		# before iteration 1.
		#
		# Fix: raise the per-attempt budget, enable pair reselection (which lets
		# the bipartition step pick a different district pair instead of giving
		# up), and retry the whole seed from a fresh random seed if it still
		# fails. The retries matter because the failure is stochastic: the same
		# instance and settings can succeed or fail depending on the draw.
		_seed_attempts = int(os.environ.get("GC_SEED_MAX_ATTEMPTS", 200000))
		_seed_retries = int(os.environ.get("GC_SEED_RETRIES", 25))
		_seed_method = partial(bipartition_tree,
							   max_attempts=_seed_attempts,
							   allow_pair_reselection=True)
		# the chain's proposal gets the same hardening so a long run cannot die
		# halfway through on an unlucky recombination step
		_recom_method = partial(bipartition_tree,
								max_attempts=_seed_attempts,
								allow_pair_reselection=True)

		s = -1
		# Balance TOTAL population (P0010001) -- the quantity the MIP constrains via
		# L <= sum_v p_v x_vj <= U, where self.population is P0010001. This routine
		# previously balanced P0030001 (voting-age population) instead. Because VAP
		# is only ~77% of total population and that share varies across districts, a
		# plan balanced to +/-0.5% on VAP is generally OUTSIDE +/-0.5% on total
		# population, so its districts violate [L,U] and the plan is unusable both as
		# a feasibility certificate and as a MIP warm start. (The GerryChain routine
		# used for the tradeoff curve, below, already used P0010001 correctly.)
		my_updaters = {'population': updaters.Tally('P0010001', alias='population')}
		_pop_target = sum(self.G.nodes[i]['P0010001'] for i in self.G.nodes()) / self.k

		# Time-boxed, retrying seed. See gc_seed.py for why this is a subprocess
		# with a wall-clock budget rather than a bigger max_attempts: the seeding
		# runtime is heavy-tailed, so restarting beats waiting.
		from gc_seed import seed_partition
		_json = ('./raw_data/' + self.parcel_level + '/json/' + self.state + '_'
				 + ('counties' if self.parcel_level == 'county' else 'tracts') + '.json')
		start = seed_partition(_json, self.k, 'P0010001',
							   self.population_deviation / 2,
							   log=lambda m: print(m, flush=True))
		if start is None:
			raise RuntimeError(
				'GerryChain could not seed %s %s at +/-%.2f%% deviation within '
				'the seeding budget. No balanced recursive bipartition was '
				'found; this is a property of the instance, not a crash.'
				% (self.state, self.parcel_level,
				   100 * self.population_deviation / 2))

		initial_partition = GeographicPartition(self.G, start, updaters = my_updaters)

		proposal = partial(recom,
						pop_col='P0010001',
						pop_target=_pop_target,
						epsilon=self.population_deviation/2,
						node_repeats=2,
						method=_recom_method
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
		prefilter_hits = 0

		# CSR for the cheap (no-MIP) combinatorial fixing pre-screen
		if not hasattr(self, '_csr_fix'):
			self._csr_fix = CSRGraph(self.G.number_of_nodes(), self.G.edges())

		nodes_to_not_be_checked = fixed_vertices + skiped_vertices
		nodes_to_be_checked = list(set(self.G.nodes()) - set(nodes_to_not_be_checked))

		cardinality = len(nodes_to_be_checked)
		if nodes_to_be_checked == []:
			return []
		nodes_to_be_fixed = []

		# Worklist as a LIFO stack (matching the original's nodes_to_be_checked[-1]
		# order) with a 'done' set for O(1) lazy removal, plus an incrementally
		# maintained 'allowed' membership mask for the induced subgraph
		# (= every vertex except those already fixed). Previously the induced-vertex
		# set and the [False]*n mask were rebuilt from scratch every iteration --
		# O(n) each, i.e. O(n^2) total at tract scale -- as were two list-comprehension
		# removals. The fix/skip decisions and their order are unchanged.
		n_nodes = self.G.number_of_nodes()
		fixed_set = set(fixed_vertices)
		allowed = [False] * n_nodes
		for u in self.G.nodes():
			if u not in fixed_set:
				allowed[u] = True
		stack = list(nodes_to_be_checked)
		done = set()

		# Phase A -- iterate the cheap combinatorial pre-screen to a fixpoint
		# (Proposition "Validity of iterative fixing" in the Online Supplement).
		# Each fix removes a parcel from `allowed`, shrinking the radius-s balls of
		# the remaining parcels and possibly exposing further fixable parcels, so we
		# re-scan until a full sweep fixes nothing new. This is sound for any order,
		# fixes the maximal set the pre-screen certifies far more cheaply than the
		# per-parcel solve, and shrinks the workload of the CP-SAT phase below.
		if option == 'minority':
			prescreen_changed = True
			while prescreen_changed:
				prescreen_changed = False
				for v in nodes_to_be_checked:
					if v in done:
						continue
					if can_fix_minority(self._csr_fix, v, self.s, allowed,
					                    self.population, self.L, self.U,
					                    self.minority_population,
					                    self.voting_age_population, self.f):
						nodes_to_be_fixed.append(v)
						allowed[v] = False
						done.add(v)
						prefilter_hits += 1
						iterations += 1
						prescreen_changed = True

		# Per-pass budget on the number of CP-SAT feasibility solves (pass 3). The
		# cheap pre-screen -- which does the bulk of the fixing -- runs unlimited;
		# once this many borderline parcels have been decided by CP-SAT we stop the
		# pass, leaving the rest unfixed. This hard-bounds the cost of the expensive
		# passes (large tracts at mid/high s, where most parcels are feasible and
		# each solve also rebuilds an O(n^2) distance structure). Sound -- we only
		# ever fix provably-infeasible parcels -- and low-s passes fix almost
		# everything through the pre-screen without reaching this budget, so
		# ell_s^fix is unaffected.
		MAX_CPSAT_SOLVES_PER_PASS = 150
		cpsat_solves = 0

		while stack:
			v = stack.pop()
			if v in done:
				continue

			# -------- cheap combinatorial pre-screen (no MIP) --------
			# If v provably cannot lie in any feasible minority-majority district
			# of diameter <= s, fix it immediately and skip the MIP. Sound:
			# can_fix_minority only returns True when the radius-s ball around v
			# cannot supply population L and minority share f (validated in
			# test_prefilter.py). The MIP would fix exactly these vertices too.
			if option == 'minority' and can_fix_minority(self._csr_fix, v, self.s, allowed,
			                    self.population, self.L, self.U,
			                    self.minority_population,
			                    self.voting_age_population, self.f):
				nodes_to_be_fixed.append(v)
				allowed[v] = False
				done.add(v)
				prefilter_hits += 1
				iterations += 1
				continue
			# ---------------------------------------------------------

			# radius-s ball of v within the induced subgraph via array BFS
			# (replaces the slow per-node nx.ego_graph). bfs1 restricts traversal
			# to `allowed` nodes and caps depth at s; reached nodes are q1[:t1].
			self._csr_fix.bfs1(v, allowed=allowed, cap=self.s)
			ball = [self._csr_fix._q1[i] for i in range(self._csr_fix._t1)]

			# Pass 3 (distance cuts / diameter <= s): solve the per-node feasibility
			# with CP-SAT rather than the MIP. The far-pair constraints are binary
			# clauses that CP-SAT's clause propagation resolves in ~1s even on the
			# large high-s balls, whereas Gurobi routinely hits the 10s limit there
			# (leaving parcels unfixed). CP-SAT thus both runs far faster and, by
			# actually deciding those subproblems, fixes at least as many parcels.
			if pass_level == 3:
				# Budget exhausted: skip the (expensive) CP-SAT solve for v and
				# leave it unfixed, but DO NOT break -- keep looping so the cheap
				# pre-screen still fixes the remaining parcels it can. (Breaking
				# here also cut off pre-screen fixings, badly undercounting the
				# high-s passes of the large instances.)
				if cpsat_solves >= MAX_CPSAT_SOLVES_PER_PASS:
					done.add(v)
					iterations += 1
					continue
				cpsat_solves += 1
				status, witness = self._fixing_feasible_cpsat(v, ball)
				if status == 'timeout':
					print('Time limit reached')
					done.add(v)
				elif status == 'feasible':
					for i in witness:
						done.add(i)
				else:   # infeasible: v cannot lie in any such district -> fix it
					nodes_to_be_fixed.append(v)
					allowed[v] = False
					done.add(v)
				iterations += 1
				continue

			sub_graph = self.G.subgraph(ball)

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

			if m.status == gp.GRB.TIME_LIMIT:
				print('Time limit reached')
				done.add(v)   # cannot conclude; stop checking v (left unfixed)
			elif m.status == gp.GRB.OPTIMAL:
				# v lies in a feasible population-balanced minority district of
				# diameter <= s, so none of the witnessed nodes (t > 0.5, which
				# include v since t[v].lb = 1) can be fixed.
				for i in sub_graph.nodes():
					if m._t[i].x > 0.5:
						done.add(i)
			else:
				# infeasible: v cannot lie in any such district -> fix it.
				nodes_to_be_fixed.append(v)
				allowed[v] = False
				done.add(v)

			iterations += 1

		time2 = time.time()

		print('combinatorial pre-screen fixed', prefilter_hits, 'of', len(nodes_to_be_fixed), 'vertices without a MIP')

		if option == 'minority':

			self.minority_fixings = list(set(self.minority_fixings + nodes_to_be_fixed))
			print('percent_fixed', str(100 * len(self.minority_fixings) / len(self.G.nodes())), '%')
			print('TOTAL TIME', str(time2 - time1))

		if option == 'majority':

			self.majority_fixings = list(set(self.majority_fixings + nodes_to_be_fixed))
			print('percent_fixed', str(100 * len(self.majority_fixings) / len(self.G.nodes())), '%')
			print('TOTAL TIME', str(time2 - time1))

	def _fixing_feasible_cpsat(self, v, ball):
		"""Pass-3 per-node feasibility, solved with CP-SAT instead of the MIP.

		Question: is there a set T of parcels inside v's radius-s ball, containing
		v, with population in [L, U], a minority voting-age majority, and every pair
		within graph-distance s (i.e. diameter <= s)? The diameter is enforced by
		binary 'far-pair' clauses (not both endpoints of a >s-apart pair), which
		CP-SAT resolves far faster than the equivalent MIP -- typically ~1s even on
		the large high-s balls where Gurobi hits the time limit.

		Returns ('feasible', witness_nodes) / ('infeasible', None) / ('timeout', None).
		A minority-majority uses f = 1/2, encoded exactly as 2*mvap - vap >= 0.
		"""
		from ortools.sat.python import cp_model

		B = ball
		Bset = set(B)
		relabel = {u: i for i, u in enumerate(B)}
		edges = [(relabel[u], relabel[w]) for u in B for w in self.G.neighbors(u)
		         if w in Bset and relabel[u] < relabel[w]]
		csr = CSRGraph(len(B), edges)
		far = far_pairs(csr, self.s, range(len(B)))

		model = cp_model.CpModel()
		x = {u: model.NewBoolVar('x%d' % u) for u in B}
		model.Add(x[v] == 1)
		model.Add(sum(self.population[u] * x[u] for u in B) >= self.L)
		model.Add(sum(self.population[u] * x[u] for u in B) <= self.U)
		model.Add(sum((2 * self.minority_population[u] - self.voting_age_population[u]) * x[u] for u in B) >= 0)
		for (a, b) in far:
			model.AddBoolOr([x[B[a]].Not(), x[B[b]].Not()])

		solver = cp_model.CpSolver()
		# 3s per-parcel cap: parcels that cannot be decided this fast almost never
		# resolve within 10s either, so a tighter cap makes the unavoidable
		# time-outs ~3x cheaper while losing almost no fixings (a parcel that
		# times out is left unfixed in both cases -- sound, just conservative).
		solver.parameters.max_time_in_seconds = 3
		solver.parameters.num_search_workers = 8
		result = solver.Solve(model)
		if result in (cp_model.OPTIMAL, cp_model.FEASIBLE):
			return 'feasible', [u for u in B if solver.Value(x[u]) == 1]
		if result == cp_model.INFEASIBLE:
			return 'infeasible', None
		return 'timeout', None

	def save_fixing_info(self, pass_level, time, draw_map=False):
		output_string = ''

		# ensure the output directory exists (it may have been deleted to force a
		# clean recompute of the fixing cache)
		os.makedirs('./results_' + self.group + '/fixing/txt/', exist_ok=True)

		with open('./results_' + self.group + '/fixing/txt/' + self.state + '_' + self.parcel_level + '_' + '{:02d}'.format(self.s) + 'pass' + str(pass_level)+ '.txt', 'w') as doc:
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
		filename = './results_' + self.group + '/fixing/txt/'


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

			print('No cached file for pass', pass_level, '-', self.state, self.parcel_level, 's =', self.s)

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
		# Passes 2 (flow) and 3 (distance-cut) solve a per-node MIP (10s limit
		# each) for EVERY unfixed vertex. At tract scale and/or large s the
		# radius-s ball covers most of the graph, so these passes fix ~0 vertices
		# while costing thousands of timed-out MIPs (hours). Only run them on small
		# instances where they are cheap and productive. Skipping them is sound:
		# fixing only removes provably-zero variables, so fewer fixings just leaves
		# a larger (still exact) model.
		if self.G.number_of_nodes() <= 300:
			self.parcel_fixing_iteration(2)
			self.parcel_fixing_iteration(3)
		else:
			print('Skipping fixing passes 2-3 (n=%d > 300): per-node MIPs are '
			      'unproductive and slow at this scale.' % self.G.number_of_nodes())

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

		with open(filename, 'r') as fixing_file:
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
		
		with open('results_' + self.group + '/fixing_info_each_pass.csv', 'w') as doc:
			doc.write('State, parcel_level, s, pass, number of parcels, number of fixings, percent fixed, time, \n')
		
		
		with open('results_' + self.group + '/fixing_info.csv', 'w') as doc:
			doc.write('State, parcel level, k, n, m, s, number of fixings, percent fixed, time, \n')


		for file in os.listdir('./results_' + self.group + '/fixing/txt/'):
		
			if file[-4:] != '.txt':
				continue
			print(file)			
			number_of_parcels = 0
			number_of_minority_fixing_parcels = 0
			number_of_majority_fixing_parcels = 0

			f = open( './results_' + self.group +'/fixing/txt/' + file, 'r')
			print('./results_' + self.group + '/fixing/txt/' + file)
			#print(f.readline())
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
			f.close()
			with open('results_' + self.group + '/fixing_info_each_pass.csv', 'a') as doc:
				doc.write(output_string)		

		with open('results_' + self.group + '/fixing_info_each_pass.csv') as file:
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
					with open('results_' + self.group + '/fixing_info.csv', 'a') as doc:
					
						doc.write(state + ',' + parcel_level + ',' + str(congressional_districts[state.upper()]) + ',' + str(len(self.G.nodes)) +',' + str(len(self.G.edges)) + ',' + s  + ',' + number_of_parcels_fixed + ',' + percent_fixed + ',' + str(time) + '\n')

					time = 0

	# Memory budget for the OPTIONAL up-front pairwise far-pair distance cuts.
	# The lazy separator callback (compactness_callback) enforces the full
	# diameter <= s bound on its own, so these pairwise conflict cuts only tighten
	# the root relaxation -- they are NOT needed for correctness. Materializing all
	# of them is O(|far pairs| x #district labels); Gurobi stores every row even
	# when it is marked Lazy=3, which OOMs on large tracts (e.g. TX at the tract
	# level has ~24 minority labels x millions of far pairs => hundreds of millions
	# of rows). Past this budget we skip the up-front family and defer entirely to
	# the parsimonious lazy separation of minimal length-s a,b-separators (cf.
	# Salemi & Buchanan, "Parsimonious formulations for low-diameter clusters",
	# Math. Prog. Comp. 2020). The optimum is unchanged either way.
	MAX_UPFRONT_DISTANCE_CUTS = 5_000_000

	def add_lazy_distance_cuts(self, model, option, district_labels=True):

		non_fixed_nodes = list(set(self.G.nodes()) - set(self.minority_fixings))

		if non_fixed_nodes and option != "majority" and len(self.minority_district_range):
			induced_graph = self.G.subgraph(non_fixed_nodes)
			# conflict pairs = complement(G^s) on induced_graph (dist > s), via
			# fast BFS instead of nx.power + nx.complement (validated equal).
			_nodes = list(induced_graph.nodes())
			_relabel = {u: i for i, u in enumerate(_nodes)}
			_edges = [(_relabel[u], _relabel[v]) for u, v in induced_graph.edges()]
			_csr = CSRGraph(len(_nodes), _edges)
			conflict_pairs = [(_nodes[a], _nodes[b])
			                  for (a, b) in far_pairs(_csr, self.s, range(len(_nodes)))]

			n_cuts = len(conflict_pairs) * len(self.minority_district_range)
			if n_cuts > self.MAX_UPFRONT_DISTANCE_CUTS:
				print("  [distance cuts] deferring %d up-front minority cuts to the "
				      "lazy separator callback (budget %d)"
				      % (n_cuts, self.MAX_UPFRONT_DISTANCE_CUTS))
			else:
				for (u, v) in conflict_pairs:
					for j in self.minority_district_range:
						if district_labels:
							compact_constr = model.addConstr(model._X[u, j] + model._X[v, j] <= model._Z[j])
						else:
							compact_constr = model.addConstr(model._X[u, j] + model._X[v, j] <= 1)

						compact_constr.Lazy = 3

		if non_fixed_nodes and option != "minority" and len(self.majority_district_range):
			# Majority districts MAY contain minority-fixed nodes, so their
			# within-district distances must be measured on the FULL graph, not
			# the minority-fixing-induced subgraph. Using the induced subgraph
			# here wrongly forbids majority configurations that route through a
			# fixed node (full distance <= s but induced distance > s), which
			# cuts off FEASIBLE plans -- e.g. a warm-start majority district that
			# uses a fixed node as a connector ("MIP start violates Rxxxx").
			_full_csr = CSRGraph(self.G.number_of_nodes(), list(self.G.edges()))
			majority_conflict_pairs = [(u, v) for (u, v)
			                           in far_pairs(_full_csr, self.s, list(self.G.nodes()))]

			n_cuts = len(majority_conflict_pairs) * len(self.majority_district_range)
			if n_cuts > self.MAX_UPFRONT_DISTANCE_CUTS:
				print("  [distance cuts] deferring %d up-front majority cuts to the "
				      "lazy separator callback (budget %d)"
				      % (n_cuts, self.MAX_UPFRONT_DISTANCE_CUTS))
			else:
				for (u, v) in majority_conflict_pairs:
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
			# complement of the s-th power of the induced graph, built WITHOUT
			# materializing the dense power graph (nx.power blows up memory at
			# large s). far_pairs returns exactly the pairs at induced distance
			# > s (or disconnected) -- i.e. the edges of complement(G_induced^s).
			_nodes_list = list(induced_graph.nodes())
			_idx = {v: i for i, v in enumerate(_nodes_list)}
			_csr = CSRGraph(len(_nodes_list),
			                [(_idx[u], _idx[v]) for u, v in induced_graph.edges()])
			_far = far_pairs(_csr, self.s, range(len(_nodes_list)))
			complement_power_graph = nx.Graph()
			complement_power_graph.add_nodes_from(_nodes_list)
			complement_power_graph.add_edges_from(
				(_nodes_list[a], _nodes_list[b]) for a, b in _far)
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

		# conflict pairs = complement(G^s) within sub_graph, i.e. pairs with
		# dist_{sub_graph}(u,v) > s. Computed via fast BFS (validated equal to
		# nx.complement(nx.power(sub_graph, s)) in test_power.py) instead of
		# materializing the dense power graph and its complement.
		nodes = list(sub_graph.nodes())
		relabel = {u: i for i, u in enumerate(nodes)}
		edges = [(relabel[u], relabel[v]) for u, v in sub_graph.edges()]
		csr = CSRGraph(len(nodes), edges)

		for (a, b) in far_pairs(csr, self.s, range(len(nodes))):
			u = nodes[a]
			v = nodes[b]
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
			if os.path.exists('./results_' + self.group + '/fixing/txt/' + self.state + '_' + self.parcel_level + '_' + '{:02d}'.format(self.s) + 'pass' + str(pass_level) + '.txt'):
				fixing_file = open('./results_' + self.group + '/fixing/txt/' + self.state + '_' + self.parcel_level + '_' + '{:02d}'.format(self.s) + 'pass' + str(pass_level) + '.txt', 'r')
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

				oplot.savefig('./results_' + self.group + '/fixing/img/' + self.state + '_' + self.parcel_level + '_' + str(self.s) + 'pass' + str(pass_level) + '.png')

				plt.style.use('classic')
				oplot = local_copy_shape_info.plot(column='fixing_map',categorical=True, cmap = cmap, categories=my_categories, figsize = (10,12)).get_figure()
				plt.axis('off')
				plt.title('Pass level ' + str(pass_level) + ' for state ' + self.state + ' s: ' + str(self.s))
				plt.grid(False)
				oplot.patch.set_facecolor('xkcd:white')

				oplot.savefig('./results_' + self.group + '/fixing/img/no_legend_pass_' + str(pass_level) + '_' + self.state + '_' + self.parcel_level + '_' + str(self.s) + '.png')


		plt.style.use('classic')
		oplot = local_copy_shape_info_big_map.plot(column='fixing_map',categorical=True, cmap = custom_cmap, legend_kwds={'loc': 'best'}, categories=['Vertex never fixed', 'Fixed in first pass', 'Fixed in second pass', 'Fixed in third pass'], legend=True, figsize = (10,12)).get_figure()

		plt.axis('off')
		plt.title('State: ' + self.state + ', s: ' + str(self.s) + ', when which parcel got fixed')
		plt.grid(False)
		oplot.patch.set_facecolor('xkcd:white')
		plt.show()
		oplot.savefig('./results_' + self.group + '/fixing/img/' + self.state + '_' + str(self.s) + '.png')

		plt.style.use('classic')
		oplot = local_copy_shape_info_big_map.plot(column='fixing_map',categorical=True, cmap = 'Set1', categories=['Vertex never fixed', 'Fixed in first pass', 'Fixed in second pass', 'Fixed in third pass'], figsize = (10,12)).get_figure()
		plt.axis('off')
		plt.title('State: ' + self.state + ', s: ' + str(self.s) + ', when which parcel got fixed')
		plt.grid(False)
		oplot.patch.set_facecolor('xkcd:white')
		oplot.savefig('./results_' + self.group + '/fixing/img/no_legend_' + self.state + '_' + str(self.s) + '.png')

	def draw_symmetry_map(self):
		max_independent_set = compute_stability_number('read', self.G, self.state, self.k, self.s)
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
		# Independent-set symmetry-breaking constraints (Online Supplement,
		# Proposition 1): for the i-th independent-set vertex v (i = 1,...,q),
		#     y_{v, k-i+1} + sum_{j=1}^{min(i, k_m)} x_{v, j} = 1.
		# Together with the assignment constraints these fix a subset of the x, y
		# variables to zero, hence they carve out a FACE of P_0 and preserve its
		# integrality (Lemma in the main text). They are super valid, so the
		# optimum is preserved. Here i (paper, 1-indexed) = fixing + 1, the non-MM
		# target k-i+1 is the 0-indexed label self.k-1-fixing, and the MM prefix
		# min(i, k_m) = min(fixing+1, self.k_minority).
		max_independent_set = list(compute_stability_number('read', self.G, self.state, self.k, self.s, self.group))
		number_of_symmetry_breaking_allowed = min(len(max_independent_set), len(self.majority_district_range))

		i = self.k - 1
		for fixing in range(number_of_symmetry_breaking_allowed):
			v = max_independent_set[0]
			prefix = min(fixing + 1, self.k_minority)  # min(i, k_m) with i = fixing+1
			self.m.addConstr(gp.quicksum(self.m._X[v, j] for j in range(prefix)) + self.m._Y[v, i] == 1)
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

		#1c population balance: build each district's population expression ONCE
		# and reuse it for both the lower (L) and upper (U) bound (previously this
		# linear expression was constructed twice).
		pop_expr = {j: gp.quicksum(self.population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in self.minority_district_range}
		self.m.addConstrs(self.L * self.m._Z[j] <= pop_expr[j] for j in self.minority_district_range)
		self.m.addConstrs(pop_expr[j] <= self.U * self.m._Z[j] for j in self.minority_district_range)
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
			# This model MAXIMIZES the number of majority-minority districts and
			# its optimum is the upper bound k_m. On a timeout the VALID upper
			# bound is the solver's bound, not the incumbent ObjVal: the
			# incumbent is merely a feasible value, so using it would understate
			# k_m and could fix away optimal solutions. Floor to an integer,
			# since the objective counts districts.
			_bnd = self.m.getAttr(gp.GRB.Attr.ObjBound)
			self.num_minority_districts = int(math.floor(_bnd + 1e-6))
			# self.output and self.BB_node were left UNSET on this branch, so
			# generate_upper_bound_minority_table died with
			#   AttributeError: 'problem_instance' object has no attribute 'output'
			# the first time an instance hit the limit (CA majority-Latino).
			self.output = self.num_minority_districts
			self.BB_node = self.m.NodeCount
			self.time_out = True
		else:
			self.num_minority_districts = 'Model Status: ' + str(self.m.status)
			self.time_out = False
		time2 = time.time()
		
		self.minority_district_ONLY_model_time = time2 - time1

	def no_minority_district_model(self, s, relax=False, continuity=False):
		self.s = s
		time1 = time.time()
		self.m = gp.Model()
		self.m.Params.TIME_LIMIT = self.model_time_limit
		
		self.m._X = self.m.addVars(self.G.nodes(), range(self.k), vtype=gp.GRB.BINARY, name='X')

		self.m.addConstrs(gp.quicksum(self.m._X[v, j] for j in range(self.k)) == 1 for v in self.G.nodes())

		#1c
		self.m.addConstrs(self.L <= gp.quicksum(self.population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in range(self.k))
		#1c
		self.m.addConstrs(self.U >= gp.quicksum(self.population[v] * self.m._X[v, j] for v in self.G.nodes()) for j in range(self.k))
		#1d
	
	
		self.m._k_minority = self.k
		self.m._G = self.G
		self.m._k = self.k
		self.m._option = 'minority'
		self.m._s = s


		self.m.Params.lazyConstraints = 1
		self.m._minority_district_range = range(self.k)
		
		self.m._district_vars_exist = False
		self.m._L = self.L
		self.m._population = self.population
		self.m._n = len(self.G.nodes())
		self.m.optimize(compactness_callback)
	

		if self.m.status == gp.GRB.OPTIMAL:

			#self.minority_districts = [[] for i in range(int(self.m.getAttr(gp.GRB.Attr.ObjVal)))]
			districts = [[] for i in range(self.k)]
			#self.num_minority_districts = self.k




			for district_number in range(self.k):
				for vertex in self.G.nodes():
					if self.m._X[vertex, district_number].x > 0.5:
						districts[district_number].append(vertex)

			print(districts)
			self.minority_districts = districts
			self.num_minority_districts = len(districts)
			self.plotter('minority', 'results_none/images/' + str(self.state) + str(self.s))
			self.time_out = False

		elif self.m.status == gp.GRB.TIME_LIMIT:
			self.num_minority_districts = self.m.getAttr(gp.GRB.Attr.ObjVal)
			self.time_out = True
		else:
			self.num_minority_districts = 'Model Status: ' + str(self.m.status)
			self.time_out = False
		time2 = time.time()
		
		self.no_minority_district_model_time = '{0:.2f}'.format(time2-time1)

	def find_interesting_s_vals_no_minority(self):

		s = self.lower_bound_s - 1

		continue_search = True

		if not os.path.exists('results_' + self.group + '/interesting_s_vals.csv'):
			with open('results_' + self.group + '/interesting_s_vals.csv', 'a') as doc:
				doc.write('state, s, result, time out, time \n')

		output_string = ''

		while continue_search:
			s += 1

			self.no_minority_district_model(s)
			
			if self.time_out == True:
				continue_search = False

			if type(self.num_minority_districts) != int:
				continue_search = True
			
			if s == self.upper_bound_s:
				continue_search = False

			if self.k == self.num_minority_districts:
				continue_search = False
		
			current_string = str(s) + ', ' + str(self.num_minority_districts) + ',' + str(self.time_out) + ',' + str(self.no_minority_district_model_time) + '\n'
			output_string += current_string
			#self.break_points.append([s, self.num_minority_districts])

			with open('./results_' + self.group + '/interesting_s_vals.csv', 'a') as doc:
				current_string = self.state + ',' + current_string
				doc.write(current_string)
		with open('./results_' + self.group + '/interesting_s_vals/interesting_s_vals_' + self.state + '_' + self.parcel_level + '.csv', 'w') as doc:
			doc.write(output_string)

	def labeling_model(self, reduced_model, parcel_fixing, symmetry_breaking_constraint, bare=False, warm_start=None):
			# warm_start (optional): (minority_districts, majority_districts), each
			# a list of node-lists, with the majority-minority districts FIRST. It
			# is injected as a Gurobi MIP start below. Use it with a heuristic plan
			# (e.g. from GerryChain) whose diameter is <= s, so the start is feasible
			# and the solver begins from a known incumbent instead of from scratch.

			# bare=True: compare the BARE formulation only -- no up-front clique
			# distance cuts (the z-ordering IS kept; it is structurally required
			# for the separator callback to type districts). Compactness
			# (contiguity + diameter <= s) is enforced by the lazy separator
			# callback. parcel_fixing and symmetry_breaking_constraint should also
			# be False for a bare run.
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

			if symmetry_breaking_constraint and warm_start is None:
				# The independent-set symmetry-breaking constraints pin specific
				# vertices to specific district INDICES. An externally supplied warm
				# start (GerryChain / CCC) labels districts arbitrarily, so it lands
				# in a different symmetry class and violates these equalities
				# ("User MIP start violates constraint R0"). Skip them when warm
				# starting -- the model stays EXACT (symmetry breaking is only a
				# speed-up), and the warm start is then feasible.
				self.add_symmetry_constraints()


			if parcel_fixing:
				self.apply_parcel_fixing()

			if s and not bare:
				# Use the pairwise far-pair distance cuts (array BFS via far_pairs,
				# added as a Lazy=3 pool that Gurobi pulls on demand). This avoids
				# nx.find_cliques, whose maximal-clique enumeration blows up on the
				# dense complement-power graph at tract scale / large s. The lazy
				# separator callback still enforces the full diameter bound, so the
				# model remains exact; this only changes the up-front cut family.
				distance_constraints = 'pairwise_lazy'

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

			# Keep the z-ordering even when bare: it is structurally required for
			# the separate x/y formulation's separator callback to target the
			# correct district indices (it pins majority-minority districts to the
			# low indices). It is a valid symmetry breaker, not a performance trick.
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

			if warm_start is not None:
				# Inject the heuristic plan as a COMPLETE Gurobi MIP start: set every
				# x/y/z/w start explicitly (0 everywhere, then the plan's 1s), so
				# Gurobi does not have to complete a partial start (which it failed to
				# do, discarding the start). Majority-minority districts take the low
				# indices (x/z); the rest take the remaining indices (y/w).
				ws_minority, ws_majority = warm_start
				for v in self.G.nodes():
					for j in self.minority_district_range:
						self.m._X[v, j].Start = 0
					for j in range(self.k):
						self.m._Y[v, j].Start = 0
				for j in self.minority_district_range:
					self.m._Z[j].Start = 0
				for j in range(self.k):
					self.m._W[j].Start = 0
				counter = 0
				for minority_district in ws_minority:
					for node in minority_district:
						self.m._X[node, counter].Start = 1
					self.m._Z[counter].Start = 1
					counter += 1
				for majority_district in ws_majority:
					for node in majority_district:
						self.m._Y[node, counter].Start = 1
					self.m._W[counter].Start = 1
					counter += 1

			# retain all feasible solutions found during the search so the
			# benchmark can report every incumbent + its map (mode 1 keeps the
			# best PoolSolutions found without extra search effort).
			self.m.Params.PoolSearchMode = 1
			self.m.Params.PoolSolutions = 100

			#self.m.optimize()
			try:
				self.m.optimize(compactness_callback)
			except gp.GurobiError as e:
				# A solve-time Gurobi failure (most commonly 'Out of memory' on the
				# largest tracts, even after deferring the up-front distance cuts)
				# must not abort the whole sweep. Record a marker -- keeping the
				# best bound if one was established -- so the caller still writes a
				# row and the driver moves on to the next instance.
				self.model_runtime = getattr(self.m, 'runtime', 0.0)
				self.solve_status = 'error'
				self.incumbent_obj = None
				self.mip_gap = None
				try:
					self.best_bound = self.m.ObjBound
				except Exception:
					self.best_bound = None
				msg = str(e)
				self.num_minority_districts = 'OOM' if 'memory' in msg.lower() else ('Error: ' + msg)
				print("  [solve] GurobiError (%s); recording marker and continuing" % msg)
				return

			self.model_runtime = self.m.runtime

			# Capture the bounds and MIP gap for EVERY terminal status (additive:
			# does not change num_minority_districts below). When the 1-hour limit
			# is hit these let callers report the incumbent (lower bound for this
			# maximization), Gurobi's best bound (upper bound), and the gap, rather
			# than only a bare 'Timeout'.
			self.solve_status = self.m.status
			self.incumbent_obj = None
			self.best_bound = None
			self.mip_gap = None
			try:
				if self.m.SolCount > 0:
					self.incumbent_obj = self.m.ObjVal      # best feasible (LB)
					self.mip_gap = self.m.MIPGap            # relative gap
				self.best_bound = self.m.ObjBound           # best bound (UB)
			except Exception:
				pass

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
			
			#max_independent_set = compute_stability_number('read', self.G, self.state, self.k, self.s)
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

			try:
				self.m.optimize(compactness_callback)
			except gp.GurobiError as e:
				# Same guard as labeling_model: a solve-time failure (typically
				# 'Out of memory' on the largest tracts) must not abort the sweep.
				# Leave the status non-optimal so the caller records it and moves on.
				self.model_runtime = getattr(self.m, 'runtime', 0.0)
				self.gurobi_error = str(e)
				print("  [solve] GurobiError (%s); continuing" % e)
				return

			self.gurobi_error = None
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

	def generate_gerrychain_run(self, number_of_gerrychain_iter=10000):

		# Allow the chain length to be raised without editing code. More
		# iterations give the short-burst search more chances to find a plan with
		# additional majority-minority districts, so this is the knob that
		# improves RESULT QUALITY. (It is NOT the knob that fixes a seeding
		# failure; see run_GerryChain_heuristic for that.)
		number_of_gerrychain_iter = int(
			os.environ.get("GC_ITERS", number_of_gerrychain_iter))
		print('  GerryChain iterations: %d' % number_of_gerrychain_iter, flush=True)

		os.makedirs('./results_' + self.group + '/gerrychain/', exist_ok=True)

		time1 = time.time()
		gerrychain_results = self.run_GerryChain_heuristic(number_of_gerrychain_iter)
		time2 = time.time()
		gerrychain_time = '{0:.2f}'.format(time2-time1)
		gerrychain_output_string = self.state + ',' + self.parcel_level + ',' + str(number_of_gerrychain_iter) + ',' + gerrychain_time + ','


		for result in gerrychain_results:

			partition = result[0]
			gerrychain_s = result[1]
			gerrychain_num_minority_district = result[2]

			gerrychain_output_string += '(' + str(gerrychain_s) + '|' + str(gerrychain_num_minority_district) + ') '


			minority_districts_from_gerrychain = []
			majority_districts_from_gerrychain = []


			for district in range(self.k):
				total_pop_district = 0
				total_pop_minority = 0
				list_of_nodes_in_current_district = []
				for node in partition.graph:
						if partition.assignment[node] == district:
							list_of_nodes_in_current_district.append(node)

							total_pop_district += self.voting_age_population[node]
							total_pop_minority += self.minority_population[node]

				if total_pop_minority > self.f * total_pop_district:
					minority_districts_from_gerrychain.append(list_of_nodes_in_current_district)
				else:
					majority_districts_from_gerrychain.append(list_of_nodes_in_current_district)

			gerrychain_output = [minority_districts_from_gerrychain, majority_districts_from_gerrychain]

		
			with open('./results_' + self.group + '/gerrychain/' + self.state + '_' + self.parcel_level + '_(' + '{:02d}'.format(gerrychain_s) + ',' + '{:02d}'.format(gerrychain_num_minority_district) + ').pckl', 'wb') as doc:
				pickle.dump(gerrychain_output, doc)
		gerrychain_output_string += '\n'
		
		if not os.path.exists('gerrychain.csv'):
			with open('gerrychain.csv', 'w') as doc:
				doc.write('State, parcel level, gerrychain iterations, time, pareto frontier \n')

		with open('gerrychain.csv', 'a') as doc:
			doc.write(gerrychain_output_string)
		
	def make_sure_preprocess_exists(self, match_model_fixing=False):
		"""Precompute the fixing files for the current s.

		match_model_fixing=True computes ONLY what apply_parcel_fixing will actually
		load. That method skips passes 2-3 when n > 300, so computing pass 3 for a
		large instance costs 300-650 s of per-node CP-SAT and is then discarded --
		which dominated the frontier sweep, where the MIP itself often resolves in
		1-2 s."""

		result = compute_stability_number('read', self.G, self.state, self.k, self.s, self.group)
		
		self.minority_fixings = []
		self.majority_fixings = []

		# Pass 1 = population + minority (cheap pre-screen + small MIP).
		# Pass 3 = distance cuts (diameter <= s), a strong and fast certificate.
		# Pass 2 (single-commodity flow for contiguity) is intentionally skipped:
		# its per-node feasibility MIPs are by far the slowest part and, on hard
		# parcels, they hit the time limit without fixing anything anyway (so they
		# contribute almost nothing while dominating the runtime). Dropping it keeps
		# the procedure sound -- distance cuts are a valid relaxation, so a parcel is
		# only ever fixed when provably impossible -- at the cost of a few
		# contiguity-only fixings. Final fixing counts come from the pass-3 file.
		self.parcel_fixing_iteration(1)
		if match_model_fixing and self.G.number_of_nodes() > 300:
			print('  (skipping pass-3 precompute: apply_parcel_fixing discards it at n=%d > 300)'
			      % self.G.number_of_nodes())
		else:
			self.parcel_fixing_iteration(3)

	def generate_instance_info(self):
		if not os.path.exists('results_' + self.group + '/instance_info.csv'):
			with open('results_' + self.group + '/instance_info.csv', 'w') as doc:
				doc.write('State, parcel, k, n, m, diam, % minority VAP, current s, split parcels, current # minority districts')

		with open('results_' + self.group + '/instance_info.csv', 'a') as doc:
			# the enacted-plan statistics need raw_data/current_plan/<state>.shp,
			# which is not present for every state (e.g. AZ); fall back to 'NA' for
			# those two columns rather than aborting the whole instance.
			try:
				current_s, split_number, current_minority = self.get_current_plan_stat()
			except Exception as _e:
				print('  (no enacted-plan shapefile for %s; current s / splits / #MM = NA)' % self.state)
				current_s, split_number, current_minority = 'NA', 'NA', 'NA'
			string = '\n'
			string += self.state + ', ' + self.parcel_level + ', ' + str(self.k) + ', '+ str(len(self.G.nodes())) + ', ' + str(len(self.G.edges())) + ', ' + str(nx.diameter(self.G)) + ', ' + format(100 * sum(self.minority_population)/sum(self.voting_age_population), '.2f') + ',' + str(current_s) + ',' +  str(split_number) + self.error_message + ', ' + str(current_minority)
			
			doc.write(string)
			doc.close()

	def generate_upper_bound_minority_table(self):
		
		if not os.path.exists('results_' + self.group + '/upper_bound_minority_table.csv'):
			with open('results_' + self.group + '/upper_bound_minority_table.csv', 'w') as doc:
				doc.write('State, parcel, k, n, m, obj. (k_m), time (s), x & z fixed (%), w fixed (%)')

		self.k_minority = self.k
		self.minority_district_range = range(self.k)
		self.find_upper_bound_minority_districts(continuity=False, relax=False)
		self.table_3_first_pass = self.output
		self.table_3_first_pass_time = self.minority_district_ONLY_model_time
		self.minority_district_range = range(self.k)

		# Once the upper bound k_m on the number of majority-minority districts is
		# known, the remaining (k - k_m) districts cannot be majority-minority, so
		# their z and x variables are zero-fixed and their w variables one-fixed.
		# The fixed fraction is therefore (k - k_m)/k (identical for x&z and w).
		fixed_pct = (self.k - self.table_3_first_pass) / self.k * 100

		with open('results_' + self.group + '/upper_bound_minority_table.csv', 'a') as doc:
				string = '\n'
				string += self.state + ','
				string += self.parcel_level + ','
				string += str(self.k) + ','
				string += str(len(self.G.nodes)) + ','
				string += str(len(self.G.edges)) + ','
				string += str(int(self.table_3_first_pass)) + ','
				string += format(self.table_3_first_pass_time, '.2f') + ','
				string += format(fixed_pct, '.2f') + ','
				string += format(fixed_pct, '.2f')

				doc.write(string)
				doc.close()

	def generate_fixing_info(self):

		# Cap on how far above the lower bound we sweep s. The fixing procedure is
		# most valuable near ell_s -- it certifies ell_s^fix (the largest s that is
		# entirely fixed) and tightens the hardest, smallest-diameter models. At
		# large s almost every parcel is feasible, so each pass performs a per-node
		# feasibility solve for most parcels while fixing very little; those passes
		# dominate the runtime and add negligible value. We therefore sweep s only
		# up to ell_s + FIXING_SPAN (or until a pass fixes nothing, or s reaches the
		# graph diameter). Increase FIXING_SPAN for more of the decaying tail at the
		# cost of markedly longer runs on the large tract instances.
		FIXING_SPAN = 5

		condition = True
		self.s = self.lower_bound_s - 1
		s_cap = min(int(self.upper_bound_s), self.lower_bound_s + FIXING_SPAN)
		n_parcels = self.G.number_of_nodes()
		while condition:
			self.s += 1
			print('Current s: ', str(self.s))
			self.make_sure_preprocess_exists()
			num_fixed = len(self.minority_fixings)
			if num_fixed == 0:
				condition = False
			# Stop at the span cap, EXCEPT keep sweeping while every parcel is still
			# fixed (100%). That regime is cheap -- the pre-screen fixes all parcels
			# at once -- and its end is exactly ell_s^fix, the smallest s admitting a
			# majority-minority district; we therefore only stop past the cap once the
			# fixing first drops below 100%. (Matters for very sparse-minority states
			# such as MO, whose 100% region extends well past ell_s + FIXING_SPAN.)
			elif self.s >= s_cap and num_fixed < n_parcels:
				condition = False
			if self.s >= int(self.upper_bound_s):
				condition = False

	def generate_lower_bound_s_table(self):
		
		if not os.path.exists('results_' + self.group + '/lower_bound_s_table.csv'):
			with open('results_' + self.group + '/lower_bound_s_table.csv', 'w') as doc:
				doc.write('State, parcel, lower_bound_s, iterations, time')

		self.compute_lower_bound_s()

		with open('results_' + self.group + '/lower_bound_s_table.csv', 'a') as doc:
			string = '\n'
			string += self.state + ', ' + self.parcel_level + ', ' + str(self.lower_bound_s) + ', ' + str(self.num_s_iter) + ', ' + format(self.time_for_s_lower_bound, '.2f') + ',' + self.error_message + ', '
			doc.write(string)
			doc.close()
					
	def generate_point_check_table(self):
		if not os.path.exists('results_' + self.group + '/point_check.csv'):
			with open('results_' + self.group + '/point_check.csv', 'w') as doc:
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
				with open('results_' + self.group + '/point_check.csv', 'a') as doc:
					doc.write(MIP_output_string)

				if self.s >= 10:
					condition=False

	def evaluate_gerrychain(self):

		if not os.path.exists('table_8.csv'):
			with open('table_8.csv', 'w') as doc:
				doc.write('State, parcel level, s, number of minority districts, model result, model time \n')
		
		values_of_k_minority_to_investigate = range(self.k_minority)

		
		
		for file in os.listdir('results_' + self.group + '/gerrychain'):
			gerrychain_state = file[:2]
			if file[3:9] == 'tract_':
				gerrychain_parcel_level = 'tract'
			if file[3:9] == 'county':
				gerrychain_parcel_level = 'county'

			gerrychain_s = int(file[-11:-9])
			gerrychain_num_minority_district = int(file[-8:-6])

			if self.state == gerrychain_state and self.parcel_level == gerrychain_parcel_level:

				with open('results_' + self.group + '/gerrychain/' + file, 'rb') as f:
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
	
		if not os.path.exists('results_' + self.group + '/upper_bound_s_table.csv'):
			with open('results_' + self.group + '/upper_bound_s_table.csv', 'w') as doc:
				doc.write('State, parcel, k, m, n, upper bound s time, upper bound s, number of iterations, diameter \n')

		upper_bound_s_time, upper_bound_s, number_of_iterations = self.upper_bound_s_algo()
		
		
		with open('results_' + self.group + '/upper_bound_s/upper_bound_s_' + self.state + '_' + self.parcel_level + '.txt', 'w') as file:
			file.write(str(upper_bound_s))

		with open('results_' + self.group + '/upper_bound_s_table.csv', 'a') as doc:	
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

	def _fmt_symmetry_result(self):
		# Read the bound/incumbent/gap captured by labeling_model after the most
		# recent solve. Returns (obj_str, bound_str, gap_str, gap_value) where
		# obj is the incumbent (lower bound for this maximization), bound is
		# Gurobi's best bound (upper bound), gap is the relative MIP gap in %,
		# and gap_value is the numeric fractional gap (or None) used to compute
		# the with/without-symmetry gap improvement.
		if getattr(self, 'incumbent_obj', None) is not None:
			obj_str = str(int(round(self.incumbent_obj)))
		elif isinstance(self.num_minority_districts, str) and self.num_minority_districts != 'Timeout':
			# pass through terminal markers: 'Infeasible', 'OOM', 'Error: ...'
			obj_str = self.num_minority_districts
		else:
			obj_str = 'no incumbent'

		if getattr(self, 'best_bound', None) is not None:
			bound_str = '{0:.2f}'.format(self.best_bound)
		else:
			bound_str = 'NA'

		gap_value = getattr(self, 'mip_gap', None)
		gap_str = '{0:.2f}'.format(gap_value * 100.0) if gap_value is not None else 'NA'

		return obj_str, bound_str, gap_str, gap_value

	def generate_symmetry_table(self):
		"""Panel A -- INFEASIBLE regime: run at s = ell_s (set by the constructor
		when data.json leaves s unset). At ell_s the tract instances have no valid
		districting, so this measures how fast infeasibility can be PROVEN with vs
		without the symmetry-breaking constraints."""
		self._symmetry_run('results_' + self.group + '/symmetry_table.csv')

	def generate_symmetry_feasible_table(self):
		"""Panel B -- FEASIBLE regime: run at an s where a plan is known to exist
		(supplied per instance via data.json). Both solves then have a feasible
		region, so incumbents/gaps exist and gap_improvement is meaningful."""
		self._symmetry_run('results_' + self.group + '/symmetry_table_feasible.csv')

	def _symmetry_run(self, csv_path):
		# NOTE: parcel fixing is ON in both solves, matching the original symmetry
		# experiment (and Panel A's recorded rows). The with/without-symmetry arms
		# are therefore apples-to-apples; only the symmetry constraints differ.
		if not os.path.exists(csv_path):
			with open(csv_path, 'w') as doc:
				doc.write('State, parcel level, k, n, m, s, '
					+ 'obj_wo_symmetry, bound_wo_symmetry, gap_wo_symmetry(%), time_wo_symmetry, '
					+ 'obj_w_symmetry, bound_w_symmetry, gap_w_symmetry(%), time_w_symmetry, '
					+ 'gap_improvement(pp), \n')
		else:
			# RESUME: skip only if this (state, parcel level, s) triple is already
			# recorded. Keying on s as well lets the same instance appear at its
			# infeasible s (Panel A) and its feasible s (Panel B) without one
			# masking the other. Each solve pair costs up to ~2 hours, so this
			# makes long runs robust to restarts.
			with open(csv_path) as doc:
				existing = doc.read().splitlines()
			for line in existing[1:]:
				f = [c.strip() for c in line.split(',')]
				if len(f) > 5 and f[0] == self.state and f[1] == self.parcel_level and f[5] == str(self.s):
					print('  (symmetry: %s %s at s=%s already recorded -- skipping)'
					      % (self.state, self.parcel_level, self.s))
					return

		string = self.state + ',' + self.parcel_level + ',' + str(self.k) + ',' + str(len(self.G.nodes)) + ',' + str(len(self.G.edges)) + ','  + str(self.s) + ','

		# ---- without symmetry-breaking constraints ----
		time1 = time.time()
		self.labeling_model(True, True, False)
		time2 = time.time()
		no_symmetry_time = format(time2 - time1, '.2f')
		obj_wo, bound_wo, gap_wo, gap_wo_val = self._fmt_symmetry_result()

		string += obj_wo + ',' + bound_wo + ',' + gap_wo + ',' + no_symmetry_time + ','

		# ---- with symmetry-breaking constraints ----
		time1 = time.time()
		self.labeling_model(True, True, True)
		time2 = time.time()
		symmetry_time = format(time2 - time1, '.2f')
		obj_w, bound_w, gap_w, gap_w_val = self._fmt_symmetry_result()

		# gap improvement in percentage points (positive => symmetry shrank the gap)
		if gap_wo_val is not None and gap_w_val is not None:
			gap_improvement = '{0:.2f}'.format((gap_wo_val - gap_w_val) * 100.0)
		else:
			gap_improvement = 'NA'

		string += obj_w + ',' + bound_w + ',' + gap_w + ',' + symmetry_time + ',' + gap_improvement + ', \n'

		with open(csv_path, 'a') as doc:
			doc.write(string)
 	
	def _district_diameter_ok(self, nodes, s):
		"""True iff G[nodes] is connected with diameter <= s."""
		sub = self.G.subgraph(nodes)
		if sub.number_of_nodes() == 0:
			return False
		for v in nodes:
			reach = nx.single_source_shortest_path_length(sub, v, cutoff=s)
			if len(reach) < len(nodes):     # some node is >s away (or unreachable)
				return False
		return True

	def _plan_certifies(self, districts, s, k_m):
		"""True iff `districts` is a valid plan for the prescribed model at (k_m, s):
		an exact k-partition, every district population-balanced with diameter <= s,
		and at least k_m districts majority-minority. Such a plan is a CONSTRUCTIVE
		feasibility certificate, so no MIP is needed for that frontier point.
		(At least k_m suffices: the non-majority-minority districts carry no minority
		constraint, so surplus majority-minority districts may occupy the y-slots.)"""
		if len(districts) != self.k:
			return False
		seen = set()
		for d in districts:
			for v in d:
				if v in seen:
					return False
				seen.add(v)
		if len(seen) != self.G.number_of_nodes():
			return False

		mm = 0
		for d in districts:
			pop = sum(self.population[v] for v in d)
			if pop < self.L or pop > self.U:
				return False
			if not self._district_diameter_ok(d, s):
				return False
			vap = sum(self.voting_age_population[v] for v in d)
			mvap = sum(self.minority_population[v] for v in d)
			if mvap > self.f * vap:
				mm += 1
		return mm >= k_m

	def _plan_max_diameter(self, districts):
		"""Max district diameter of a plan, or None if some district is disconnected."""
		worst = 0
		for d in districts:
			sub = self.G.subgraph(d)
			if sub.number_of_nodes() == 0 or not nx.is_connected(sub):
				return None
			worst = max(worst, nx.diameter(sub))
		return worst

	def _heuristic_ub(self, plans, k_m):
		"""Smallest s at which a stored heuristic plan CERTIFIES k_m majority-minority
		districts, i.e. min over valid plans with >= k_m MM districts of the plan's
		max district diameter (a plan of diameter d certifies every s >= d). Returns
		None if no stored plan achieves k_m. This is an upper bound on s*(k_m) and
		costs no MIP."""
		best = None
		for p in plans:
			if not self._plan_certifies(p, self.G.number_of_nodes(), k_m):
				continue          # invalid partition / population / MM count
			d = self._plan_max_diameter(p)
			if d is not None and (best is None or d < best):
				best = d
		return best

	def _load_heuristic_plans(self):
		"""All stored GerryChain plans for this instance, as lists of k node-lists."""
		d = './results_' + self.group + '/gerrychain/'
		prefix = self.state + '_' + self.parcel_level + '_('
		plans = []
		if not os.path.isdir(d):
			return plans
		for fn in sorted(os.listdir(d)):
			if not (fn.startswith(prefix) and fn.endswith('.pckl')):
				continue
			try:
				with open(d + fn, 'rb') as doc:
					minority_districts, majority_districts = pickle.load(doc)
				plans.append(list(minority_districts) + list(majority_districts))
			except Exception as e:
				print('  (skipping unreadable plan %s: %s)' % (fn, e))
		return plans

	def check_points_with_prescribed_model(self):
		"""Table 6 -- the (s, #MM) tradeoff frontier, s*(k_m) = the smallest s at which
		a plan with k_m majority-minority districts exists.

		Uses all three accelerations together:
		  * heuristic (GerryChain)  -- a verified plan CERTIFIES feasibility, no MIP;
		  * variable fixing         -- shrinks the prescribed feasibility MIP;
		  * symmetry breaking       -- what makes the INFEASIBLE s values provable.
		Exploits monotonicity: a plan with k_m majority-minority districts is also a
		plan for k_m - 1, so s*(k_m) is non-decreasing in k_m. Sweeping k_m upward and
		carrying s forward therefore skips every s already known to be infeasible."""
		csv_path = 'results_' + self.group + '/point_check.csv'
		if not os.path.exists(csv_path):
			with open(csv_path, 'w') as doc:
				doc.write('State, parcel, n, m, k_m, s*, status, time, certified_by \n')

		k_m_max = self.k_minority
		ell_s = self.s
		u_s = getattr(self, 'upper_bound_s', None) or nx.diameter(self.G)
		plans = self._load_heuristic_plans()
		print('  [frontier] %d heuristic plan(s) available; ell_s=%s, u_s=%s'
		      % (len(plans), ell_s, u_s))

		# RESUME: skip (state, level, k_m) triples already recorded, and restart the
		# monotone carry from the largest feasible s* already established.
		done_k_m, s_resume = set(), None
		with open(csv_path) as doc:
			for line in doc.read().splitlines()[1:]:
				f = [c.strip() for c in line.split(',')]
				if len(f) > 6 and f[0] == self.state and f[1] == self.parcel_level:
					try:
						done_k_m.add(int(f[4]))
					except ValueError:
						continue
					# carry forward the established LOWER bound on s*; f[5] is either
					# an exact "16" or a bracket "[13,16]".
					tok = f[5].lstrip('[').split(',')[0].lstrip('>=')
					try:
						s_resume = max(s_resume or 0, int(tok))
					except ValueError:
						pass
		if done_k_m:
			print('  [frontier] resuming: k_m %s already recorded' % sorted(done_k_m))

		s_lo = s_resume if s_resume else ell_s   # monotone: s*(k_m) >= s*(k_m - 1)
		timeout_budget = 2                       # give up on a k_m after this many timeouts
		rows = []

		for k_m in range(1, k_m_max + 1):
			if k_m in done_k_m:
				continue
			self.k_minority = k_m
			self.minority_district_range = range(k_m)
			self.majority_district_range = range(k_m, self.k)

			elapsed, certified_by = 0.0, ''

			def _solve_at(s_val):
				"""Feasibility of the prescribed model at s_val: True / False / None
				(None = time limit or solver error). Only precomputes the fixing the
				model will actually use."""
				self.s = s_val
				self.make_sure_preprocess_exists(match_model_fixing=True)
				print('  [frontier] MIP: k_m=%d  s=%d' % (k_m, s_val))
				self.prescribed_labeling_model()
				if getattr(self, 'gurobi_error', None):
					return None
				if self.m.status == 2:
					return True
				if self.m.status == 3:
					return False
				return None

			# (1) FREE upper bound from the heuristic: smallest diameter at which a
			# stored plan already achieves k_m majority-minority districts.
			s_ub = self._heuristic_ub(plans, k_m)
			if s_ub is not None:
				certified_by = 'heuristic'
				print('  [frontier] k_m=%d: heuristic certifies s*<=%d' % (k_m, s_ub))

			s_lb = s_lo
			status, s_star = 'Unknown', '>=' + str(s_lb)

			# (2) EXISTENCE FIRST. At s = u_s >= diam(G) the diameter bound is
			# vacuous, so this single solve decides whether ANY plan with k_m
			# majority-minority districts exists. If not, we are done in one solve
			# instead of scanning every s from ell_s to u_s.
			if s_ub is None:
				r = _solve_at(u_s)
				elapsed += self.model_runtime
				if r is False:
					status, s_star, certified_by = 'Infeasible', '>' + str(u_s), 'MIP (vacuous s)'
				elif r is True:
					s_ub, certified_by = u_s, 'MIP'
				else:
					status, s_star, certified_by = 'Unknown', '>=' + str(s_lb), 'timeout at u_s'

			# (3) BINARY SEARCH for the smallest feasible s in (s_lb, s_ub], using
			# monotonicity: the feasible set only grows with s.
			if s_ub is not None and status != 'Infeasible':
				while s_lb < s_ub:
					mid = (s_lb + s_ub) // 2
					r = _solve_at(mid)
					elapsed += self.model_runtime
					if r is True:
						s_ub = mid
						certified_by = 'MIP'
					elif r is False:
						s_lb = mid + 1
						certified_by = 'MIP'
					else:
						break                     # timeout: cannot narrow further
				if s_lb >= s_ub:
					status, s_star = 'Feasible', str(s_ub)
				else:
					status, s_star = 'Bracketed', '[%d,%d]' % (s_lb, s_ub)

			rows.append('%s,%s,%d,%d,%d,%s,%s,%.2f,%s\n'
			            % (self.state, self.parcel_level, self.G.number_of_nodes(),
			               self.G.number_of_edges(), k_m, s_star, status, elapsed, certified_by))
			with open(csv_path, 'a') as doc:
				doc.write(rows[-1])

			# Monotonicity: s*(k_m) is non-decreasing in k_m, so the next k_m can
			# start from this k_m's established LOWER bound (every s below it is
			# already proven infeasible, and infeasibility only worsens as k_m grows).
			s_lo = max(s_lo, s_lb)
			if status == 'Unknown' or 'error' in certified_by:
				break                    # k_m unresolved => larger k_m is no easier

		# restore instance state
		self.s = ell_s
		self.k_minority = k_m_max
		self.minority_district_range = range(k_m_max)
		self.majority_district_range = range(k_m_max, self.k)

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
	
	group = 'hispanic'
	#group = 'black'
	#group = 'none'
	#dataset = 'hispanic_paper'
	#dataset = 'paper'
	#dataset = 'paper'
	#dataset = 'hispanic_only'
	dataset = 'single'

	for request in data[dataset]:
		
		print('State: ', request['state'])
		print('Parcel level: ', request['parcel_level'])

		#if request['state'] == 'TX':
		#	continue
		
		instance = problem_instance(request['state'], request['parcel_level'], request['s'], group)

 
		if instance.error_message != '':
			print('ERROR: ', instance.error_message)
			break
		
		
		#instance.create_upper_bound_s_table()
		#instance.find_upper_bound_s()
		#instance.generate_lower_bound_s_table()
		#instance.find_interesting_s_vals_no_minority()
		#instance.create_fixing_csv_file()
		#instance.black_population_map()
		#instance.draw_symmetry_map()
		#instance.generate_gerrychain_run(10)
		#instance.make_sure_preprocess_exists()
		#instance.generate_instance_info()
		#instance.generate_instance_info()
		#instance.generate_upper_bound_minority_table()
		#instance.generate_lower_bound_s_table()
		#instance.generate_point_check_table()
		#instance.make_sure_preprocess_exists()
		#instance.check_points_with_prescribed_model()
		#instance.generate_lower_bound_s_table()
		#instance.generate_instance_info()
		#instance.create_upper_bound_s_table()
		#instance.generate_symmetry_table()
		#instance.generate_fixing_info()
		#instance.check_points_with_prescribed_model()
		#instance.short_bursts()
		#instance.lower_bound_s_minority_dist_only()
		#instance.draw_fixing_map()
		#instance.black_population_map()
		#instance.find_upper_bound_minority_districts()
		
	file.close()