import sys, os, math, string, random, unittest
import stats

class NotImplementedException(Exception):
	def __init__(self):
		myvar = 1
		
def probabilityOfFixation(Ne, s):
	return (1 - math.exp(2*s))/(1 - math.exp(2*Ne*s))

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
		self._sequence = sequence
		self._fitness = 1.0

	def copy(self):
		e = EvolvableSequence(self.sequence)
		return e

	@property
	def fitness(self):
		return self._fitness
	
	@fitness.setter
	def fitness(self,f):
		self._fitness = f
	
	@property
	def sequence(self):
		return self._sequence

	def spawn(self, mutator):
		res = SpawnResult()
		(mut_sequence, mutations) = mutator.mutate(self._sequence)
		offs = EvolvableSequence(mut_sequence)
		res.offspring = offs
		res.mutations = mutations
		return res

	@property
	def key(self):
		return self._sequence

	def __str__(self):
		return self._sequence

	def __eq__(self,x):
		return self._sequence == x._sequence
	
	def __getitem__(self,i):
		return self._sequence[i]
	
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
		self._refcount = 0
		self.coalescent = False
		self._id = None

	def addReference(self):
		self._refcount += 1
		return self._refcount

	def removeReference(self):
		self._refcount -= 1
		return self._refcount
	
	@property
	def count(self):
		return self._refcount

	@property
	def coalescent(self):
		return self.coalescent
	
	@coalescent.setter
	def coalescent(self, bool):
		self.coalescent = bool
	
	@property
	def id(self):
		return self._id

	@id.setter
	def id(self, the_id):
		self._id = the_id
	
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
				self.removeEntry(entry.parent)
	
	# Coalescent if 
	def isCoalescent(self, entry):
		res = entry.coalescent
		if not entry.coalescent:
			parent = entry.parent
			if not parent is None and self.isCoalescent(parent):
				if parent.count == 1:
					entry.coalescent = True
					res = True
		return res

class FixationResults:
	def __init__(self):
		self.fixed = None
		self.time_to_fixation = None
	

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
		if len(self._members) <= self.population_size:
			self._members.append(entry)
		else:
			self.inject(entry)

	def createOffspring(self, organism_entry, mutate=True):
		"""Creates offspring entry. Does not add to population."""
		# Begin assuming no mutation will occur
		new_entry = organism_entry
		if mutate:
			# Reproduce with mutation
			spawnres = organism_entry.organism.spawn(self.mutator)
			if spawnres.mutated:
				# Create an entry for the GeneBank
				new_entry = self.genebank.createEntry(spawnres.offspring, organism_entry, spawnres.offspring.fitness, self.generation_count)
		# Add to the GeneBank
		self.genebank.addEntry(new_entry)
		return new_entry

	def populate(self, organism):
		"""Add population_size copies of the specified organism"""
		# Make a parent for the whole population
		parent = organism
		parent_entry = self.genebank.createEntry(parent, None, parent.fitness, self.generation_count)
		parent_entry.coalescent = True
		# All new population
		for i in range(self.population_size):
			self.addMember(self.createOffspring(parent_entry, mutate=False))
		# Remove the parent -- it's no longer in the population
		self.genebank.removeEntry(parent_entry)

	def makeCumulativeProbabilities(self):
		total_fitness = 0.0
		# Cache fitness and compute total fitness
		for m in self._members:
			m.cache_fitness = m.organism.fitness
			total_fitness += m.cache_fitness
		# Create cumulative probabilities
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
			# Remove the previous generation
			for m in self._members:
				self.genebank.removeEntry(m)
			# Replace population with the buffer
			self._members = self._members_buffer
	
	def evolveUntilFixationOrLossOf(self, entry):
		"""Evolve until organism specified in entry is lost or fixed.
		Fixation is defined as when the entry becomes coalescent.
		"""
		start_generations = self.generations
		while entry.count > 0 and not self.genebank.isCoalescent(entry):
			#print entry.count, self.genebank.isCoalescent(entry), self.histogram()
			self.evolve(1)
		# Observe the results
		res = FixationResults()
		if entry.count == 0:
			res.fixed = False
			res.time_to_fixation = None
		else:
			res.fixed = True
			res.time_to_fixation = self.generations - start_generations
		return res

	def erase(self):
		"""Get rid of all information in the population."""
		self._members = []
		self.genebank.erase()
	
	def inject(self, organism):
		"""Put the supplied organism into a randomly chosen spot in the population as a spontaneous mutant."""
		slot = random.choice(range(len(self._members)))
		slot_entry = self._members[slot]
		new_entry = self.genebank.createEntry(organism, slot_entry, organism.fitness, self.generation_count)
		# Add injected organism as if it were a spontaneous mutant
		#new_entry.parent = slot_entry.parent
		#self.genebank.addEntry(new_entry.parent)
		self._members[slot] = new_entry
		self.genebank.addEntry(new_entry)
		#print self.genebank.isCoalescent(new_entry)
		# Kill existing organism
		self.genebank.removeEntry(slot_entry)
		return new_entry

	###########
	## Analytical methods
	#############
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
