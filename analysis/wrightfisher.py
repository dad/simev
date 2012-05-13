import sys, os, math, string, random, unittest
import stats

class NotImplementedException(Exception):
	def __init__(self):
		myvar = 1

class SpawnResult:
	"""Record storing results of replication with mutation."""
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
	"""An interface for a class that can reproduce with mutation and has a fitness function.""" 
	def fitness(self):
		"""Provides a measure of reproductive fitness."""
		raise NotImplementedException, "Must override"
	
	def spawn(self, mutator):
		"""Replication with mutation."""
		raise NotImplementedException, "Must override"

class MutationInfo:
	"""Record storing information about mutations."""
	def __init__(self, location, from_base, to_base):
		self.location = location
		self.from_base = from_base
		self.to_base = to_base

class Mutator:
	"""Interface for a class that can mutate a sequence."""
	def mutate(self, sequence):
		raise NotImplementedException, "Must override"

class SimpleMutator(Mutator):
	"""Unbiased nucleotide mutations."""
	def __init__(self, per_site_mutation_rate, alphabet = 'ATGC'):
		self.per_site_mutation_rate = per_site_mutation_rate
		self.alphabet = alphabet
		self.mut_choices = dict([(x,alphabet.replace(x,'')) for x in alphabet])

	def mutate(self, sequence):
		"""Mutate the sequence with the specified per-site mutation rate."""
		mut_sequence = [x for x in sequence]
		mutations = []
		for l in range(len(mut_sequence)):
			if random.random() < self.per_site_mutation_rate:
				from_base = mut_sequence[l]
				mut_sequence[l] = random.choice(self.mut_choices[mut_sequence[l]])
				mutations.append(MutationInfo(l, from_base, mut_sequence[l]))
		return ''.join(mut_sequence), mutations

class EvolvableSequence(Evolvable):
	"""Base implementation of a sequence class that can evolve."""
	def __init__(self, sequence):
		self.sequence = sequence

	def copy(self):
		e = EvolvableSequence(self.sequence)
		return e

	def fitness(self):
		return 1.0

	def spawn(self, mutator):
		res = SpawnResult()
		(mut_sequence, mutations) = mutator.mutate(self.sequence)
		offs = EvolvableSequence(mut_sequence)
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
	
	def __getitem__(self,i):
		return self.sequence[i]
	
def pickIndexByProbability(cum_probs, p):
	assert p>=0.0 and p<=1.0
	i = 0
	while p>cum_probs[i]:
		i+=1
	return i

# Track lineages. When they arise, start keeping a count of their frequencies; when they're extinguished
# This is borrowed heavily from Claus's GeneBank.
class GeneBankEntry:
	def __init__(self, organism, parent, fitness, birthtime):
		self.organism = organism
		self.parent = parent
		self.fitness = fitness
		self.birthtime = birthtime
		self.refcount = 0
		self.coalescent = False
		self._id = None

	def addReference(self):
		self.refcount += 1
		return self.refcount

	def removeReference(self):
		self.refcount -= 1
		return self.refcount

	def setCoalescent(self, bool):
		self.coalescent = bool

	@property
	def id(self, the_id):
		self._id = the_id
	
	@property
	def id(self):
		return self._id

	def __str__(self):
		return "{} ({})".format(self.organism, self.refcount)

	def __eq__(self, other):
		# DAD: not quite sure how strict to be here
		return self.organism == other.organism
	
class GeneBank:
	def __init__(self):
		# Table maps keys (IDs) to GeneBankEntries.
		self.table = {}
		# An ID counter. Increments with each entry created.
		self.maxID = 0

	# Factory model: GeneBank creates entries which are owned and distributed by it.
	def createEntry(self, organism, parent_entry, fitness, birthtime):
		entry = GeneBankEntry(organism, parent_entry, fitness, birthtime)
		if not parent_entry is None:
			parent_entry.addReference()
		self.table[self.maxID] = entry
		entry.id = self.maxID
		self.maxID += 1
		return entry

	def getEntry(self, id):
		return self.table.get(id)
		
	def addEntry(self, entry):
		entry.addReference()

	def removeEntry(self, entry):
		if not entry is None:
			if entry.removeReference()==0:
				# Refcount is zero: time to remove the organism from our world.
				del self.table[entry.id]
				if not entry.parent is None:
					self.removeEntry(entry.parent)

class Population:
	def __init__(self, population_size, mutator):
		self.population_size = population_size
		# List of members of the population. These are GeneBankEntry objects, which wrap organisms to allow reference-counting.
		self._members = []
		# A buffer for building subsequent generations of the population.
		self._members_buffer = []
		# Class supporting the .mutate() method.
		self.mutator = mutator
		# The gene bank: the set of organisms along the line of descent.
		self.genebank = GeneBank()
		# The number of generations since the beginning of the simulation.
		self.generation_count = 0

	def addMember(self, entry):
		# Add to the population
		self._members.append(entry)

	def createOffspring(self, organism_entry, mutate=True):
		"""Creates but does not add offspring."""
		# Begin assuming no mutation will occur
		new_entry = organism_entry
		if mutate:
			# Reproduce with mutation
			spawnres = organism_entry.organism.spawn(self.mutator)
			if spawnres.mutated:
				# Create an entry for the GeneBank
				new_entry = self.genebank.createEntry(spawnres.offspring, organism_entry, spawnres.offspring.fitness(), self.generation_count)
		# Add to the GeneBank
		self.genebank.addEntry(new_entry)
		return new_entry

	def populate(self, organism):
		"""Add population_size copies of the specified organism"""
		# Make a parent for the whole population
		parent = organism
		parent_entry = self.genebank.createEntry(parent, None, parent.fitness(), self.generation_count)
		parent_entry.setCoalescent(True)
		# All new population
		for i in range(self.population_size):
			self.addMember(self.createOffspring(parent_entry, mutate=False))
		# Remove the parent -- it's no longer in the population
		self.genebank.removeEntry(parent_entry)

	def makeCumulativeProbabilities(self):
		total_fitness = 0.0
		for m in self._members:
			m.cache_fitness = m.organism.fitness()
			total_fitness += m.cache_fitness
		self._members.sort(key = lambda x: x.cache_fitness, reverse=True)
		cum_probs = []
		cum_prob = 0.0
		for m in self._members:
			cum_prob += m.cache_fitness/total_fitness
			cum_probs.append(cum_prob)
		return cum_probs

	def evolve(self, num_generations):
		for n in xrange(num_generations):
			# Clear out the buffer
			self._members_buffer = []
			cum_probs = self.makeCumulativeProbabilities()
			# Increment generations
			self.generation_count += 1
			# Replicate into the new generation
			for nm in xrange(len(self._members)):
				# Pick parent according to fitness: Wright-Fisher sampling
				parent_entry = self._members[pickIndexByProbability(cum_probs, random.random())]
				# Reproduce with mutation
				spawn_entry = self.createOffspring(parent_entry)
				# Insert offspring into new gen
				self._members_buffer.append(spawn_entry)
				# Track
				#self.genebank.add(spawn_result, self.generations)
			# Replace population with the buffer
			self._members = self._members_buffer

	def erase(self):
		"""Get rid of all information in the population."""
		self._members = []
		self.genebank.erase()

	## Analytical methods
	def frequency(self, entry):
		"""Get the frequency of a particular organism."""
		return self.table[item.key()]

	def histogram(self):
		"""Return a sorted count of the number of each type of organism in the population"""
		ids = [e.id for e in self._members]
		hist = []
		for id in set(ids):
			hist.append(ids.count(id))
		return sorted(hist, reverse=True)

	def dominantOrganism(self):
		ids = [e.id for e in self._members]
		id_counts = sorted([(ids.count(id),id) for id in ids], reverse=True)
		return self.genebank.getEntry(id_counts[0][1])
		
	def averageFitness(self):
		total_fitness = 0.0
		for m in self._members:
			m.cache_fitness = m.organism.fitness()
			total_fitness += m.cache_fitness
		return total_fitness/self.size()
				
	@property
	def size(self):
		assert len(self._members) == self.population_size
		return self.population_size

	@property
	def generations(self):
		return self.generation_count

	@property
	def members(self):
		for m in self._members:
			yield m

	def count(self, entry):
		ids = [e.id for e in self._members]
		return ids.count(entry.id)
		
		
	def __str__(self):
		for m in self._members:
			print m
