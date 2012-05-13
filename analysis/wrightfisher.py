import sys, os, math, string, random, unittest
import stats

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

class EvolvableSequence(Evolvable):
	def __init__(self, sequence, alphabet = 'ATGC'):
		self.sequence = sequence
		self.alphabet = alphabet
		self.mut_choices = dict([(x,alphabet.replace(x,'')) for x in alphabet])

	def fitness(self):
		return 1.0

	def spawn(self, per_site_mutation_rate):
		res = SpawnResult()
		mut_sequence = [x for x in self.sequence]
		for l in range(len(mut_sequence)):
			if random.random() < per_site_mutation_rate:
				mut_sequence[l] = random.choice(self.mut_choices[mut_sequence[l]])
				res.mutated = True
		offs = EvolvableSequence(''.join(mut_sequence), self.alphabet)
		res.offspring = offs
		return res

	def __str__(self):
		return self.sequence

	def eq(self,x):
		return self.sequence == x.sequence
	
def pickIndexByProbability(cum_probs, p):
	assert p>=0.0 and p<=1.0
	i = 0
	while p>cum_probs[i]:
		i+=1
	return i

class LineageTracker:
	def __init__(self):
		self.frequency_table = {}
	def add(self, item, generation):
		pass

class Population:
	def __init__(self, mutation_rate):
		self.members = []
		self.mutation_rate = mutation_rate
		self.tracker = LineageTracker()
		self.generation_count = 0

	def addMember(self, member):
		self.tracker.add(member, self.generations)
		self.members.append(member)

	def makeCumulativeProbabilities(self):
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
		return cum_probs

	def evolve(self, num_generations):
		for n in xrange(num_generations):
			nextgen = []
			cum_probs = self.makeCumulativeProbabilities()
			# Increment generations
			self.generation_count += 1
			# Replicate into the new generation
			for nm in xrange(len(self.members)):
				# Pick parent according to fitness
				parent = self.members[pickIndexByProbability(cum_probs, random.random())]
				# Reproduce with mutation
				spawn_result = parent.spawn(self.mutation_rate)
				# Insert offspring into new gen
				nextgen.append(spawn_result.offspring)
				# Track
				self.tracker.add(spawn_result, self.generations)
			# Replace population
			self.members = nextgen
				
	@property
	def size(self):
		return len(self.members)

	@property
	def generations(self):
		return self.generation_count
		
