import sys, os, math, string, random, collections, unittest
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


## Mutation

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

	# DAD: use FitnessEvaluator pattern?
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
		self._coalescent = False
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
		return self._coalescent
	
	@coalescent.setter
	def coalescent(self, bool):
		self._coalescent = bool
	
	@property
	def id(self):
		return self._id

	@id.setter
	def id(self, the_id):
		self._id = the_id
	
	def __str__(self):
		return "{} (f={}, par={}, n={}, birth={}, coal={})".format(self.id, self.fitness, self.parent.id, self._refcount, self.birthtime, self.coalescent)

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

	def __getitem__(self, id):
		return self.table[id]
		
	def addEntry(self, entry):
		entry.addReference()

	def removeEntry(self, entry):
		if not entry is None:
			if entry.removeReference()==0:
				# Refcount is zero: time to remove the organism from our world.
				del self.table[entry.id]
				self.removeEntry(entry.parent)
	
	def lineage(self, entry):
		lin = [entry]
		e = entry
		while not e.parent is None:
			lin.append(e.parent)
			e = e.parent
		return lin
	
	def __str__(self):
		return ','.join(["{}:{}".format(v.id,v.count) for v in self.table.values()])
	
class FixationResults:
	def __init__(self):
		self.fixed = None
		self.time_to_fixation = None

class Population:
	"""An evolving population."""
	def __init__(self, max_population_size, mutator, fitness_evaluator):
		raise NotImplementedException, "Must override."

	def evolve(self, num_generations):
		pass
	
	def populate(self, genotype):
		"""Populates the population with instances of the provided genotype."""
		pass
		
class SampleCounter(collections.Counter):
	def choice(self):
		"""Choose an element with probability equal to k/n, where k is the element count and n = sum(k)"""
		n = sum(self.values())
		ordered = self.most_common()
		r = random.randint(0,n-1)
		i = 0
		k = ordered[0][1]
		while r>k:
			i+=1
			k += ordered[i][1]
		return ordered[i][0]
		
	def elements(self):
		"""Iterate over elements, repeating each according to its frequency"""
		for (key, count) in self.items():
			for n in xrange(count):
				yield key
		

class WrightFisherPopulation(Population):
	"""An evolving population.
	A population consists of N (.population_size) individuals, which are themselves instances of M < N genotypes.
	Evolution proceeds by Wright-Fisher sampling. In generation t, N offspring individuals are generated. The probability that one
	individual offspring has, as its parent, an individual in generation t-1 is equal to the fitness of the parent.
	
	To accommodate large N, the implementation uses reference counting. If there are n_i
	instances of genotype i in the population, then the population will store a single instance of i with a count of n_i.
	
	To track genotypes
	"""
	def __init__(self, population_size, mutator):
		self.population_size = population_size
		# List of members of the population. These are GeneBankEntry objects, which wrap organisms to allow reference-counting.
		self._members = SampleCounter()
		# A buffer for building subsequent generations of the population.
		self._members_buffer = SampleCounter()
		# Class supporting the .mutate() method.
		self.mutator = mutator
		# The gene bank: the set of organisms along the line of descent.
		self.genebank = GeneBank()
		# The number of generations since the beginning of the simulation.
		self.generation_count = 0

	def addMember(self, entry):
		# Add to the population. Does not enforce population size.
		self._members[entry.id] += 1

	def removeMember(self, entry):
		# Remove from the population. Does not enforce population size.
		self._members[entry.id] -= 1

	def createOffspring(self, organism_entry, mutate=True):
		"""Creates offspring entry and add to GeneBank. Does not add to population."""
		# Begin assuming no mutation will occur
		new_entry = organism_entry
		if mutate:
			# Reproduce with mutation
			spawnres = organism_entry.organism.spawn(self.mutator)
			if spawnres.mutated:
				# Create an entry for the GeneBank
				new_entry = self.genebank.createEntry(spawnres.offspring, organism_entry, spawnres.offspring.fitness, self.generation_count)
		# DAD: this is wrong. Parent should be the parental genotype.
		#new_entry.parent = organism_entry
		# Add to the GeneBank
		self.genebank.addEntry(new_entry)
		return new_entry

	def populate(self, organism):
		"""Add population_size copies of the specified organism"""
		# Make a parent for the whole population
		parent = organism
		parent_entry = self.genebank.createEntry(parent, None, parent.fitness, self.generation_count)
		parent_entry.coalescent = True
		self.genebank.addEntry(parent_entry)
		# Now create the population individuals, who will have parent_entry as their parent
		offspring = self.genebank.createEntry(organism, parent_entry, organism.fitness, self.generation_count)
		# All new population
		for i in range(self.population_size):
			offs = self.createOffspring(offspring, mutate=False)
			self.addMember(offs)
		# Remove the parent -- it's no longer in the population
		self.genebank.removeEntry(parent_entry)
		return parent_entry

	def makeCumulativeProbabilities(self):
		"""Make sampleable list weighted by fitness and count (degeneracy) in population"""
		total_fitness = 0.0
		# Cache fitness and compute total fitness
		sorted_entries = []
		for (mk, m_count) in self._members.items():
			m = self.genebank[mk]
			# Cache the fitness...DAD we may need a FitnessEvaluator instance in here.
			m.cache_fitness = m.organism.fitness
			total_fitness += m.cache_fitness*m_count
			sorted_entries.append((m.cache_fitness*m_count, m))
		assert total_fitness > 0.0, "Total fitness = {} <= 0.0, aborting".format(total_fitness)
		# Create cumulative probabilities
		sorted_entries.sort(reverse=True)
		cum_probs = []
		cum_prob = 0.0
		for (fit, m) in sorted_entries:
			cum_prob += fit/total_fitness
			cum_probs.append(cum_prob)
		return cum_probs, sorted_entries

	def evolve(self, num_generations):
		for n in xrange(num_generations):
			# Clear out the buffer
			self._members_buffer = SampleCounter()
			cum_probs, sorted_entries = self.makeCumulativeProbabilities()
			# Increment generations
			self.generation_count += 1
			# Replicate into the new generation
			for nm in xrange(self.population_size):
				# Pick parent according to fitness: Wright-Fisher sampling
				parent_entry = sorted_entries[pickIndexByProbability(cum_probs, random.random())][1]
				# Reproduce with mutation; add to genebank
				spawn_entry = self.createOffspring(parent_entry)
				# Insert offspring into new generation
				self._members_buffer[spawn_entry.id] += 1
				#print nm, spawn_entry.id, parent_entry.id
			# Remove the previous generation
			# DAD: should shortcut and remove whole counts
			for m in self.members:
				#print m.count
				self.genebank.removeEntry(m)
			# Replace population with the buffer
			self._members = self._members_buffer
			self._members_buffer = None
	
	def evolveUntilFixationOrLossOf(self, entry):
		"""Evolve until organism specified in entry is lost or fixed.
		Fixation is defined as when the entry becomes coalescent.
		"""
		start_generations = self.generations
		while self.count(entry) > 0 and not self.isCoalescent(entry):
			#print self.count(entry), self.isCoalescent(entry), self.histogram()
			self.evolve(1)
		# Observe the results
		res = FixationResults()
		if self.count(entry) == 0:
			res.fixed = False
			res.time_to_fixation = None
		else:
			res.fixed = True
			res.time_to_fixation = self.generations - start_generations
		return res

	# Assess whether given organism, assumed to be in the current population, is
	# coalescent, i.e., shares a genotype with the last common ancestor of the population.
	# Organisms are born non-coalescent, can switch to become coalescent, and then
	# never lose that status.
	def isCoalescent(self, entry):
		# DAD: can cache...set entry.coalescent = True, only check LCA if not coalescent.
		# That is optimization; do not do prematurely... ;)
		lca = self.lastCommonAncestor()
		return entry == lca
	
	def old_stuff(self, entry):
		pid = -1
		if not entry.parent is None:
			pid = entry.parent.id
		print "check coal for id={} p={} c={} (n={})".format(entry.id, pid, entry.coalescent, self.count(entry))
		res = entry.coalescent
		if not entry.coalescent:
			# Check for coalescence. If parent is coalescent, and has a count of 1,
			# then we are coalescent as well.
			parent = entry.parent
			if not parent is None and parent.coalescent:
				if self.count(parent) == 1:
					entry.coalescent = True
					res = True
		return res
	
	
	def lastCommonAncestor(self):
		"""Get the last common ancestor of the population."""
		ids = set(self._members.keys())
		gb = self.genebank
		lineages = [[x.id for x in gb.lineage(gb[k])] for k in ids]
		s = [set(x) for x in lineages]
		intersect = reduce(lambda x,y: x.intersection(y), s)
		return gb[max(intersect)]
		while len(ids)>1:
			#print "lca:", ids
			parents = [gb[id].parent for id in ids]
			new_ids = [par.id for par in parents if not par is None]
			ids = set(new_ids)
		#print "lca:", ids
		return self.genebank[list(ids)[0]]

	def erase(self):
		"""Get rid of all information in the population."""
		self._members = []
		self.genebank.erase()
	
	def choice(self):
		return self.genebank[self._members.choice()]
	
	def inject(self, organism):
		"""Put the supplied organism into a randomly chosen spot in the population as a spontaneous mutant."""
		slot_entry = self.choice()
		new_entry = self.genebank.createEntry(organism, slot_entry.parent, organism.fitness, self.generation_count)
		self.genebank.addEntry(new_entry)
		# Add injected organism as if it were a spontaneous mutant
		#new_entry.parent = slot_entry.parent
		#self.genebank.addEntry(new_entry.parent)
		self.addMember(new_entry)
		self.removeMember(slot_entry)
		# Kill existing organism
		#self.genebank.removeEntry(slot_entry)
		return new_entry

	###########
	## Analytical methods
	#############
	def frequency(self, entry):
		"""Get the frequency of a particular organism."""
		return self.table[item.key()]

	def histogram(self):
		"""Return a sorted count of the number of each type of organism in the population"""
		return [x[1] for x in self._members.most_common()]

	def dominantOrganism(self):
		id = self._members.most_common(1)[0][0]
		return self.genebank[id]
		
	def averageFitness(self):
		total_fitness = 0.0
		for m in self._members:
			m.cache_fitness = m.organism.fitness()
			total_fitness += m.cache_fitness
		return total_fitness/self.size()
	
	## Properties
	@property
	def size(self):
		assert len(self._members) == self.population_size
		return self.population_size

	@property
	def generations(self):
		return self.generation_count

	@property
	def members(self):
		"""Return an iterator over the members of the population."""
		for m in self._members.elements():
			yield self.genebank[m]

	def count(self, entry):
		return self._members[entry.id]
		
	def __str__(self):
		for m in self._members:
			print m
