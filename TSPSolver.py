#!/usr/bin/python3
from random import randint
from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np
from TSPClasses import *
import heapq
import itertools



class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution, 
		time spent to find solution, number of permutations tried during search, the 
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''
	
	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()

		for c in cities:
			for j in cities:
				if c != j:
					print(str(c._name) + "-" + str(j._name) + ": " + str(c.costTo(j)))

		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results
	
	
	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''
		
	def branchAndBound( self, time_allowance=60.0 ):
		default = self.greedy(time_allowance)
		bssf = default['cost']
		cities = self._scenario.getCities()
		ncities = len(cities)

		##############
		# START TOUR #
		##############
		start_time = time.time()

		# Create initial matrix
		c_m = []
		for i in range(ncities):  # rows -> start city: from
			row = []
			for j in range(ncities):  # columns -> end city: to
				row.append(cities[i].costTo(cities[j]))
			c_m.append(row)

		c_m = np.array(c_m)
		c_m, c = reduce(c_m,0)

		# Create Priority Queue w/ initial state & default bssf
		Q = PQ(c_m, bssf, c, randint(0,ncities-1), ncities)

		while Q.queue:# and time.time()-start_time < time_allowance:
			Q.expand(Q.choose())

		end_time = time.time()

		############
		# END TOUR #
		############

		route = []
		if Q.solution is not None:
			for i in range(len(Q.solution.path)):
				route.append(cities[Q.solution.path[i]])
			route = TSPSolution(route)
		else:
			route = default["soln"]

		result = {
			'cost': Q.bssf if Q.solution is not None else default['cost'],
			'time': end_time - start_time,
			'count': Q.n_sol,
			'soln': route,
			'max': Q.max_len,
			'total': Q.states,
			'pruned': Q.pruned
		}
		return result

	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''
		
	def fancy( self,time_allowance=60.0 ):
		cities = self._scenario.getCities()
		ncities = len(cities)
		route = []
		nsoln = 0

		##############
		# start tour #
		##############
		start_time = time.time()
		F = Fancy(cities, ncities)
		while time.time() - start_time <= time_allowance:
			F.find_tour(start_time, time_allowance)
			nsoln += 1
			if F.route != np.inf:
				route = TSPSolution(F.convert(F.route))
				if route.cost < np.inf:
					break



		end_time = time.time()
		############
		# end tour #
		############

		if end_time - start_time > time_allowance:
			route = None
			nsoln = None

		result = {
			'cost': route.cost if route is not None else math.inf,
			'time': end_time - start_time,
			'count': nsoln,
			'soln': route,
			'max': None,
			'total': None,
			'pruned': None
		}
		return result



	''' <summary>
			This is the entry point for the greedy solver, which you must implement for 
			the group project (but it is probably a good idea to just do it for the branch-and
			bound project as a way to get your feet wet).  Note this could be used to find your
			initial BSSF.
			</summary>
			<returns>results dictionary for GUI that contains three ints: cost of best solution, 
			time spent to find best solution, total number of solutions found, the best
			solution found, and three null values for fields not used for this 
			algorithm</returns> 
		'''

	def greedy(self, time_allowance=60.0):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)

		##############
		# START TOUR #
		##############
		start_time = time.time()

		Q = Queue(cities)

		route = []
		solutions = 0
		foundTour = False
		while not foundTour and time.time()-start_time < time_allowance:
			Q.choose()
			if len(Q.queue) == 1:
				solutions += 1
				Q.add_last()
				route = TSPSolution(Q.route)
				if route.cost != math.inf:
					foundTour = True
				else:
					# reset for new tour
					Q.reset()
		end_time = time.time()

		result = {
			'cost': route.cost,
			'time': end_time - start_time,
			'count': solutions,
			'soln': route,
			'max': None,
			'total': None,
			'pruned': None
		}

		print("done")

		return result

#######################################
# FANCY ALGORITHM CLASSES AND METHODS #
#######################################
class Fancy:
	def __init__(self, c, n):
		self.cities = c
		self.ncities = n

		self.start = 0
		self.route = [self.start]
		self.cost = 0


	def find_tour(self, t_start, t_allowed):
		self.start = randint(0, self.ncities-1)
		self.route = [self.start]
		self.cost = 0

		for i in range(self.ncities - 1):
			if self.route == np.inf or time.time() - t_start >= t_allowed:
				break
			next = self.find_farthest()
			self.insert(next)


	def find_farthest(self):
		city = 0
		max_dist = 0
		for i in self.route:
			for j in range(self.ncities):
				if j not in self.route and j != i:
					dist = self.cities[i].costTo(self.cities[j])
					if max_dist < dist < np.inf:
						max_dist = dist
						city = j
		return city

	def insert(self, city):
		# insert at shortest
		min_dist = math.inf
		min_route = np.inf
		city_index = 0

		for i in range(len(self.route)):
			r = self.route.copy()
			r.insert(i, city)
			route = TSPSolution(self.convert(r))
			if route.cost < min_dist:
				min_dist = route.cost
				min_route = r
		self.route = min_route

	def convert(self, r):
		route = []
		for i in r:
			route.append(self.cities[i])
		return route






########################################
# GREEDY ALGORITHM CLASSES AND METHODS #
########################################
class Queue:
	def __init__(self, c):
		self.queue = []
		self.route = []

		self.cities = c
		self.start = randint(0, len(c) - 1)
		self.location = self.start

		self.make_queue(c)

	def make_queue(self, cities):
		assert (type(cities[0]) == City)
		for i in range(len(cities)):
			self.queue.append(self.Node(cities[i]._index))
		self.queue[self.location].cost = 0

	def choose(self):
		curr = self.queue.pop(self.location)
		source = self.cities[curr.id]
		mei = 0		# minimum edge index
		cost = source.costTo(self.cities[self.queue[mei].id])
		for i in range(1, len(self.queue)):
			dest = self.cities[self.queue[i].id]
			if source.costTo(dest) < cost:
				mei = i
				cost = source.costTo(self.cities[self.queue[mei].id])

		self.location = mei
		self.queue[mei].cost = curr.cost + cost
		self.queue[mei].prev = curr.id

		self.route.append(self.cities[curr.id])

	def add_last(self):
		curr = self.queue.pop(self.location)
		self.route.append(self.cities[curr.id])

	def reset(self):
		self.queue = []
		self.route = []

		self.start = randint(0, len(self.cities) - 1)
		self.location = self.start

		self.make_queue(self.cities)

	def __str__(self):
		out = ""
		for n in self.queue:
			out += str(n)
		return out

	class Node:
		def __init__(self, i):
			self.id = i
			self.cost = math.inf
			self.prev = None

		def __str__(self):
			return "\n{ id: " + str(self.id) + " cost: " + str(self.cost) + " prev: " + str(self.prev) + " }"


##################################################
# BRANCH AND BOUND ALGORITHM CLASSES AND METHODS #
##################################################


####################
# MATRIX REDUCTION #
####################
def reduce(m, c):
	t = m		# adjacency matrix
	cost = c  	# base cost of matrix

	# reduce rows
	for i in range(len(t)):
		min_l = min(t[i])
		if min_l > 0 and min_l != np.inf:
			cost += min_l
			for c in range(len(t[i])):
				if t[i,c] != np.inf:
					t[i,c] -= min_l

	# reduce columns
	for j in range(len(t)):
		min_l = min(t[:,j])
		if min_l > 0 and min_l != np.inf:
			cost += min_l
			for i in range(len(m[j])):
				t[i,j] -= min_l

	return t, cost

########################
# PRIORITY QUEUE CLASS #
########################
class PQ:
	def __init__(self, m, b, c, i, n):
		# properites
		self.queue = []			# array of states
		self.start = i 			# start city
		self.ncities = n  		# number of cities to visit
		self.solution = None  	# array of solutions

		self.states = 1			# number of states created - determines state id
		self.bssf = 0			# best search so far - cost of best solution
		self.pruned = 0			# number of pruned states
		self.max_len = 0  		# max length of the queue
		self.n_sol = 0			# number of solutions created - determines solution id
		self.lb = 0				# lower bound

		# initialization
		self.bssf = b
		S = self.State(self.states,m,c,i,[i])
		self.insert(S)
		self.lb = self.queue[0].cost

	def choose(self):
		return self.delete_min()

	def delete_min(self):
		# pop states that are >= BSSF due to solution updates
		self.pop_obsolete()

		min_i = 0
		for i in range(1, len(self.queue)):
			# drill deeper states w/ longer paths
			if len(self.queue[i].path) > len(self.queue[min_i].path):
				min_i = i
			elif len(self.queue[i].path) == len(self.queue[min_i].path):
				if self.queue[i].cost < self.queue[min_i].cost:
					min_i = i
		return self.queue.pop(min_i)

	def pop_obsolete(self):
		pop = []
		for i in range(len(self.queue)):
			if self.queue[i].cost >= self.bssf:
				pop.append(i)
		for i in range(len(pop), 0,-1):
			self.queue.pop(pop[i - 1])
			self.pruned += 1

	def expand(self, S):  	# S - state
		for j in range(len(S.matrix[S.city])):
			if S.matrix[S.city,j] != np.inf:
				P = self.create_state(S, S.city, j)

				if P.cost < self.bssf:
					if len(P.path) == self.ncities:
						self.update_solution(P)
					else:
						self.insert(P)
				else:
					self.pruned += 1

	def create_state(self, s, f, t):  	# s - state to expand, f - (f)rom start city, t - (t)o end city
		n = np.copy(s.matrix)
		cost = s.cost + n[f,t]

		# set
		for j in range(len(n[f])):
			n[f,j] = np.inf
		for i in range(len(n)):
			n[i,t] = np.inf
		n[t,f] = np.inf

		n,cost = reduce(n,cost)
		path = s.path.copy()
		path.append(t)

		self.states += 1
		return self.State(self.states, n, cost, t, path)

	def insert(self, S):	# S - state
		self.queue.append(S)
		self.max_len = max(self.max_len, len(self.queue))  # Is queue larger?

	def update_solution(self, S):	# S - state
		self.bssf = S.cost
		self.solution = S
		self.n_sol += 1

	###############
	# STATE CLASS #
	###############
	class State:
		def __init__(self, n, m, b, c, p):
			self.s_id = n			# state id
			self.matrix = m			# adjacency matrix of edge costs
			self.cost = b			# cost of tour for this state
			self.city = c			# current city
			self.path = p			# full path from starting city to current city




