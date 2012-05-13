import sys, os, math, string, random, unittest

class NotImplementedException(Exception):
	def __init__(self):
		myvar = 1

class SpawnResult:
	def __init__(self):
		self.offspring = None
		self.mutated = False

class Evolvable:
	def __init__(self):
		pass
	def fitness(self):
		raise NotImplementedException, "Must override"
	def spawn(self, mutation_rate):
		raise NotImplementedException, "Must override"
	
def pickIndexByProbability(cum_probs, p):
	assert p>=0.0 and p<=1.0
	i = 0
	while p>cum_probs[i]:
		i+=1
	return i

class Population:
	def __init__(self, num_members, member_type, mutation_rate):
		self.members = []
		for n in xrange(num_members):
			m = member_type()
			self.members.append(m)
		self.mutation_rate = mutation_rate

	def evolve(self, num_generations):
		for n in xrange(num_generations):
			nextgen = []
			total_fitness = 0.0
			for m in self.members:
				m.cache_fitness = m.fitness()
				total_fitness += m.cache_fitness
			self.members.sort(key = lambda x: x.cache_fitness, reverse=True)
			cum_probs = []
			cum_prob = 0.0
			for m in self.members:
				cum_prob += m.cache_fitness/total_fitness
				cum_probs.append(cum_prob)
			for nm in xrange(len(self.members)):
				# Pick parent according to fitness
				parent = self.members[pickIndexByProbability(cum_probs, random.random())]
				# Reproduce with mutation
				spawn_result = parent.spawn(self.mutation_rate)
				# Insert offspring into new gen
				nextgen.append(spawn_result.offspring)
			# Replace population
			self.members = nextgen
			# Increment
				
		
	
		
