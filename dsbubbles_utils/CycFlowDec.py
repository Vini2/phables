import numpy as np
import copy

# Simple cycle flow decomposition class
class CycFlowDec:
	def __init__(self,F,state,tol):
		self.walks = [dict(),dict()]
		self.walks[0][(state,)] = Walk(state)
		self.tick = 0
		self.opp_tick = 1
		self.cycles = dict()
		self.S = np.zeros([F.shape[0],F.shape[0]])
		self.A = []
		for j in range(F.shape[0]):
			s = np.sum(F[:,j])
			if s != 0:
				self.S[:,j] = F[:,j]/s
			adj = []
			for k in range(F.shape[0]):
				if(F[k,j] > 0):
					adj.append(k)
			self.A.append(adj)
		self.F = F
		self.burnin = 0
		self.tol = tol
		self.alpha = -1

	# Run decomposition without cycle logging for burnin steps
	# and then with cycle logging for nstep
	def run(self,burnin,nstep):
		self.cycles.clear()
		self.alpha = -1
		self.burnin = burnin
		j = 0
		ntot = nstep + self.burnin
		while(j < ntot):
			self.step()
			j += 1

	# Perform one percolation step from self.walks[self.tick]
	# to self.walks[self.opp_tick]
	def step(self):
		for path,walk in self.walks[self.opp_tick].items():
			walk.flow = 0
		for path,walk in self.walks[self.tick].items():
			residual_flow = 0
			for state in self.A[path[-1]]:
				if (walk.flow > self.tol or state in walk.visited):
					self.progress_walk(path,walk,state)
				else:
					p = self.S[state,path[-1]]
					residual_flow += walk.flow * p
			if path in self.walks[self.opp_tick]:
				self.walks[self.opp_tick][path].flow += residual_flow
			else:
				next_walk = copy.deepcopy(walk)
				next_walk.flow = residual_flow
				self.walks[self.opp_tick][path] = next_walk
		self.swap_ticks()
		if self.burnin > 0:
			self.burnin -= 1

	# Update self.walk[self.opp_tick] and self.cycles given a state to append
	# to walk with a corresponding path
	def progress_walk(self,path,walk,state):
		p = self.S[state,path[-1]]
		next_flow = walk.flow * p
		if state in walk.visited:
			if self.burnin == 0:
				cycle = walk.get_cycle(path,state)
				if cycle in self.cycles:
					self.cycles[cycle] += next_flow
				else:
					self.cycles[cycle] = next_flow
			next_path = path[:walk.visited[state]+1]
			if next_path in self.walks[self.opp_tick]:
				self.walks[self.opp_tick][next_path].flow += next_flow
			else:
				j_vec = range(walk.visited[state]+1,walk.N)
				next_walk = copy.deepcopy(walk)
				for j in j_vec:
					next_walk.visited.pop(path[j])
					next_walk.N -= 1
				self.walks[self.opp_tick][next_path] = next_walk
		else:
			path = path + (state,)
			if path in self.walks[self.opp_tick]:
				self.walks[self.opp_tick][path].flow += next_flow
			else:
				next_walk = copy.deepcopy(walk)
				next_walk.add(state,next_flow)
				self.walks[self.opp_tick][path] = next_walk

	# Swap self.tick and self.opp_tick
	def swap_ticks(self):
		if self.tick == 0:
			self.tick = 1
			self.opp_tick = 0
		else:
			self.tick = 0
			self.opp_tick = 1

	# Calculate self.alpha for scaling self.cycles to match the
	# self.F
	def calc_alpha(self):
		cycle_sum = 0
		for cycle in self.cycles:
			cycle_sum += len(cycle) * self.cycles[cycle]
		self.alpha = np.sum(self.F)/cycle_sum

	# Scale self.cycles to match self.F
	def scale_cycles(self):
		self.calc_alpha()
		for cycle in self.cycles:
			self.cycles[cycle] *= self.alpha

	# Return the edge-wise mean relative error (MRE) of self.cycles
	# and self.F for edges greater than tol
	def calc_MRE(self,tol):
		SRE = 0
		if self.alpha == -1:
			self.scale_cycles()
		cycle_sum = 0
		N = 0
		for j in range(self.S.shape[0]):
			for k in range(self.S.shape[0]):
				if (self.F[k,j] > tol):
					cycle_sum = 0
					for cycle in self.cycles:
						if self.cycle_has_edge(cycle,(j,k)):
							cycle_sum += self.cycles[cycle]
					SRE += abs(
						self.F[k,j] - 
						cycle_sum
					)/self.F[k,j]
					N += 1
		return SRE/N
	
	# Return True if cycle has edge, False otherwise
	def cycle_has_edge(self,cycle,edge):
		N = len(cycle)
		for j in range(N-1):
			if (cycle[j] == edge[0] and
				cycle[j+1] == edge[1]
				):
				return True
		if (cycle[-1] == edge[0] and
			cycle[0] == edge[1]
			):
			return True
		return False

# Walk stack information class
class Walk:
	def __init__(self,state):
		self.visited = dict() # indices
		self.visited[state] = 0
		self.flow = 1
		self.N = 1
	
	# Add state with flow to self
	def add(self,state,flow):
		self.visited[state] = self.N
		self.flow = flow
		self.N += 1

	# Return the cycle from the previous visitation of
	# state on path
	def get_cycle(self,path,state):
		index = self.visited[state]
		mindex = (
			index +
			np.argmin(path[index:])
		)
		cycle = (
			path[mindex:] + 
			path[index:mindex]
		)
		return cycle
