import sys, os, math, string, random, unittest
import stats

class NotImplementedException(Exception):
	def __init__(self):
		myvar = 1

class SpawnResult:
	def __init__(self):
		self.offspring = None
		self.mutations = []

	@property
	def mutated(self):
		return len(self.mutations)>0

	@property
	def num_mutations(self):
		return len(self.mutations)

class Evolvable:
	def __init__(self):
		pass
	def fitness(self):
		raise NotImplementedException, "Must override"
	def spawn(self, mutation_rate):
		raise NotImplementedException, "Must override"

class MutationInfo:
	def __init__(self, location, from_base, to_base):
		self.location = location
		self.from_base = from_base
		self.to_base = to_base

class SimpleMutator:
	def __init__(self, per_site_mutation_rate, alphabet = 'ATGC'):
		self.per_site_mutation_rate = per_site_mutation_rate
		self.alphabet = alphabet
		self.mut_choices = dict([(x,alphabet.replace(x,'')) for x in alphabet])

	def mutate(self, sequence):
		mut_sequence = [x for x in sequence]
		mutations = []
		for l in range(len(mut_sequence)):
			if random.random() < self.per_site_mutation_rate:
				from_base = mut_sequence[l]
				mut_sequence[l] = random.choice(self.mut_choices[mut_sequence[l]])
				mutations.append(MutationInfo(l, from_base, mut_sequence[l]))
		return ''.join(mut_sequence), mutations
		

class EvolvableSequence(Evolvable):
	def __init__(self, sequence):
		self.sequence = sequence

	def fitness(self):
		return 1.0

	def spawn(self, mutator): #per_site_mutation_rate):
		res = SpawnResult()
		(mut_sequence, mutations) = mutator.mutate(self.sequence)
		offs = EvolvableSequence(''.join(mut_sequence))
		res.offspring = offs
		res.mutations = mutations
		return res

	@property
	def key(self):
		return self.sequence

	def __str__(self):
		return self.sequence

	def __eq__(self,x):
		return self.sequence == x.sequence
	
def pickIndexByProbability(cum_probs, p):
	assert p>=0.0 and p<=1.0
	i = 0
	while p>cum_probs[i]:
		i+=1
	return i

# Track lineages. When they arise, start keeping a count of their frequencies; when they're extinguished
# This is borrowed heavily from Claus's GeneBank.
class GeneBankEntry:
	def __init__(self, organism, parent, generations):
		self.organism = organism
		self.parent = parent
		self.generations = generations
		self.refcount = 1

	def addReference(self):
		self.refcount += 1
		return self.refcount

	def removeReference(self):
		self.refcount -= 1
		return self.refcount
	
class GeneBank:
	def __init__(self):
		self.table = {}
		
	def add(self, item, generation):
		key = item.key
		try:
			entry = self.table[key]
			entry.addReference()
		except KeyError, ke:
			# DAD: implement
			self.table[key] = 1
			

	def frequency(self, item):
		return self.frequency_table[item.key()]

class GeneBankAnalyzer:
	def __init__(self):
		pass

class Population:
	def __init__(self, mutator):
		self.members = []
		self.mutator = mutator
		self.genebank = GeneBank()
		self.generation_count = 0

	def addMember(self, member):
		self.genebank.add(member, self.generations)
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
				# Pick parent according to fitness: Wright-Fisher sampling
				parent = self.members[pickIndexByProbability(cum_probs, random.random())]
				# Reproduce with mutation
				spawn_result = parent.spawn(self.mutator)
				# Insert offspring into new gen
				nextgen.append(spawn_result.offspring)
				# Track
				#self.genebank.add(spawn_result, self.generations)
			# Replace population
			self.members = nextgen
				
	@property
	def size(self):
		return len(self.members)

	@property
	def generations(self):
		return self.generation_count
		
