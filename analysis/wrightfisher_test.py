import sys, os, math, string, random, unittest, time
import stats, numpy
import wrightfisher as wf

class TraitOrganism(wf.Evolvable):
	def __init__(self):
		self.trait = 1.0
		self.key = randomSequence(12,string.letters)

	@property
	def fitness(self):
		return self.trait
	
	def spawn(self, mutation_rate):
		res = wf.SpawnResult()
		p = random.random()
		offs = TraitOrganism()
		res.offspring = offs
		if p < mutation_rate:
			offs.trait = self.trait + random.random()
			res.mutated = True
		else:
			offs.trait = self.trait
		return res

def randomSequence(n, alphabet):
	return ''.join(stats.sample_wr(alphabet, n))


## TEST CASES

class test001(unittest.TestCase):
	"""background"""
	def test_members(self):
		alphabet='ATGC'
		N = random.randint(10,100)
		pop = wf.WrightFisherPopulation(N, wf.SimpleMutator(0.01,alphabet), wf.SequenceFitnessEvaluator())
		pop.populate(TraitOrganism())
		for m in pop.members:
			self.assertTrue(pop.count(m)==N)
	
	def test_choice(self):
		sc = wf.SampleCounter()
		for al in string.letters:
			sc[al] = 0
		sc['a'] = 10
		c = sc.choice()
		self.assertTrue(c == 'a')
		
	def test_pick_index_by_prob(self):
		n = 100
		v = range(n)
		sv = sum(v)
		# Set up cumulative probabilities that we know for certain.
		cum_probs = []
		for x in v:
			cum_probs.append(x/float(n))
		# Now know that
		for x in range(n):
			#print x, wf.pickIndexByProbability(cum_probs, x/float(n))
			self.assertTrue(wf.pickIndexByProbability(cum_probs, x/float(n)) == x)

class test002(unittest.TestCase):
	"""basic run"""
	def test_basic_run(self):
		alphabet='ATGC'
		p = wf.WrightFisherPopulation(10, wf.SimpleMutator(0.01,alphabet), wf.SequenceFitnessEvaluator())
		p.populate(TraitOrganism())
		p.evolve(1000)
		#print len(p.members)
		#for o in p.members:
		#	print o.trait

class test003(unittest.TestCase):
	"""Sequence"""
	def test_run(self):
		alphabet = 'ATGC'
		p = wf.WrightFisherPopulation(10, wf.SimpleMutator(0.01,alphabet), wf.SequenceFitnessEvaluator())
		p.populate(wf.EvolvableSequence(randomSequence(20,alphabet)))
		p.evolve(1000)
		self.assertTrue(True)

class test004(unittest.TestCase):
	def test_zero_mutation(self):
		alphabet = 'ATGC'
		e = wf.EvolvableSequence(randomSequence(100,alphabet))
		# Zero mutation rate copy
		off = e.spawn(wf.SimpleMutator(0.0)).offspring
		# Should yield identical sequences
		self.assertTrue(e == off)

class test005(unittest.TestCase):
	"""Dominant organism"""
	def test_dominant_organism(self):
		alphabet = 'ATGC'
		pop = wf.WrightFisherPopulation(100,wf.SimpleMutator(0.0001,alphabet), wf.SequenceFitnessEvaluator())
		seq = wf.EvolvableSequence(randomSequence(100,alphabet))
		pop.populate(seq)
		pop.evolve(100)
		dom_seq_entry = pop.dominantOrganism()
		counts = []
		for m in pop.members:
			counts.append((pop.count(m),m))
		counts.sort(reverse=True)
		self.assertTrue(dom_seq_entry == counts[0][1])

class test006(unittest.TestCase):
	"""Histogram"""
	def test_histogram(self):
		alphabet = 'ATGC'
		n = 1000
		pop = wf.WrightFisherPopulation(n,wf.SimpleMutator(0.0001,alphabet), wf.SequenceFitnessEvaluator())
		seq = wf.EvolvableSequence(randomSequence(90,alphabet))
		pop.populate(seq)
		h = pop.histogram()
		self.assertTrue(h[0][1] == n)
		for i in range(10):
			pop.evolve(1)
			h = dict(pop.histogram())
			self.assertTrue(sum(h.values()) == n)

class test007(unittest.TestCase):
	"""Mutation info"""
	def test_mutation_info(self):
		alphabet = 'ATGC'
		mut = wf.SimpleMutator(0.1,alphabet)
		seq = wf.EvolvableSequence(randomSequence(50,alphabet))
		mutres = seq.spawn(mut)
		newseq = [x for x in mutres.offspring.sequence]
		for m in mutres.mutations:
			self.assertTrue(seq[m.location]==m.from_base)
			self.assertTrue(newseq[m.location]==m.to_base)
			newseq[m.location] = m.from_base
		self.assertTrue(''.join(newseq) == seq.sequence)
		
class test008(unittest.TestCase):
	"""Assured fixation"""
	def test_assured_fixation(self):
		alphabet = 'ATGC'
		dx = 0.1
		Ne = 20
		mu = 0.001
		pop = wf.WrightFisherPopulation(Ne,wf.SimpleMutator(mu,alphabet), wf.SequenceFitnessEvaluator())
		seq = wf.EvolvableSequence(randomSequence(100,alphabet))
		random.seed(3)
		numpy.random.seed(11)
		# No fitness
		seq.fitness = 0.0
		parent = pop.populate(seq)
		# Check to ensure coalescence
		self.assertTrue(pop.isCoalescent(parent))
		mutseq = wf.EvolvableSequence(randomSequence(100,alphabet))
		mutseq.fitness = 1.0
		mutentry = pop.inject(mutseq)
		cum_probs, sorted_entries = pop.makeCumulativeProbabilities()
		res = pop.evolveUntilFixationOrLossOf(mutentry)
		# Everyone else has zero fitness -- fixation is assured.
		self.assertTrue(res.fixed)
		# Fixation must have happened in a single generation.
		self.assertTrue(res.time == 1)
		# Parent must be the injected sequence
		fixed_entry = pop.choice()
		self.assertTrue(fixed_entry == mutentry)
		
class test009(unittest.TestCase):
	"""Predicting probability of fixation"""
	def test_prob_of_fixation(self):
		alphabet = 'ATGC'
		dx = 0.1
		Ne = 1000
		mu = 0.0001
		predicted_fixation_probability = wf.probabilityOfFixation(Ne, dx)
		n_fixations = 0
		n_total = 0
		n_trials = 3*int(1/dx)
		random.seed(111)
		#print "Aborting because this test takes a long time -- please do run occasionally!"
		#return
		for i in range(n_trials):
			pop = wf.WrightFisherPopulation(Ne,wf.SimpleMutator(mu,alphabet), wf.SequenceFitnessEvaluator())
			seq = wf.EvolvableSequence(randomSequence(100,alphabet))
			seq.fitness = 1.0
			parent = pop.populate(seq)
			# Check to ensure coalescence
			self.assertTrue(pop.isCoalescent(parent))
			mutseq = wf.EvolvableSequence(randomSequence(100,alphabet))
			mutseq.fitness = seq.fitness + dx
			mutentry = pop.inject(mutseq)
			res = pop.evolveUntilFixationOrLossOf(mutentry)
			n_total += 1
			if res.fixed:
				n_fixations += 1
		est_prob = n_fixations/float(n_total)
		p = predicted_fixation_probability
		exp_fixations = n_trials*p
		sd = math.sqrt(n_trials*p*(1.0-p))
		# Confirm that number of fixations is within 2 SD's of expectation
		self.assertTrue(n_fixations <= (exp_fixations+2*sd))
		self.assertTrue(n_fixations >= (exp_fixations-2*sd))
			
class test010(unittest.TestCase):
	"""LCA"""
	def test_lca(self):
		alphabet = 'ATGC'
		pop = wf.WrightFisherPopulation(100,wf.SimpleMutator(0.0001,alphabet), wf.SequenceFitnessEvaluator())
		pop.populate(wf.EvolvableSequence(randomSequence(100,alphabet)))
		for n in range(10):
			pop.evolve(1)
			#print pop.genebank
			e = pop.lastCommonAncestor()
			#print e
		
class test011(unittest.TestCase):
	"""pop size"""
	def test_pop_size(self):
		alphabet = 'ATGC'
		pop = wf.WrightFisherPopulation(100,wf.SimpleMutator(0.0001,alphabet), wf.SequenceFitnessEvaluator())
		pop.populate(wf.EvolvableSequence(randomSequence(100,alphabet)))
		for n in range(10):
			pop.evolve(1)
			self.assertTrue(sum(pop._members.values()) == pop.population_size)

class test012(unittest.TestCase):
	"""simple counting"""
	def test_simple_counting(self):
		alphabet = 'ATGC'
		dx = 0.1
		Ne = 20
		mu = 0.0001
		pop = wf.WrightFisherPopulation(Ne,wf.SimpleMutator(mu,alphabet), wf.SequenceFitnessEvaluator())
		seq = wf.EvolvableSequence(randomSequence(100,alphabet))
		random.seed(3)
		# No fitness
		seq.fitness = 0.0
		parent = pop.populate(seq)
		i = 0
		for m in pop.members:
			i += 1
		self.assertTrue(i==Ne)
		# Check to ensure coalescence
		self.assertTrue(pop.isCoalescent(parent))
		mutseq = wf.EvolvableSequence(randomSequence(100,alphabet))
		mutseq.fitness = 1.0
		mutentry = pop.inject(mutseq)
		i = 0
		for m in pop.members:
			i += 1
		self.assertTrue(i==Ne)

class test013(unittest.TestCase):
	"""large populations"""
	def test_large_pop(self):
		alphabet = 'ATGC'
		dx = 0.1
		mu = 0.00001
		n_gens = 100
		random.seed(3)
		seq = wf.EvolvableSequence(randomSequence(100,alphabet))
		for i in range(5):
			tstart = time.time()
			Ne = 10**i
			pop = wf.WrightFisherPopulation(Ne,wf.SimpleMutator(mu,alphabet), wf.SequenceFitnessEvaluator())
			pop.populate(seq)
			pop.evolve(n_gens)
			tend = time.time()
			print "# evolved {} generations at Ne={} (t={} sec)".format(n_gens, Ne, round(tend-tstart,3))
	
class test014(unittest.TestCase):
	"""SampleCounter"""
	def test_sample_average(self):
		sc = wf.SampleCounter()
		x = string.letters
		sumsq = 0
		ct = 0
		for i in range(len(x)):
			sc[x[i]] = i
			sumsq += i*i
			ct += i
		#print sc
		avg = sc.average(lambda y: x.index(y))
		self.assertTrue(avg==sumsq/float(ct))

class test015(unittest.TestCase):
	"""average fitness"""
	def test_average_fitness(self):
		alphabet = 'ATGC'
		dx = 0.1
		mu = 0.0001
		base_fitness = 0.01
		n_gens = 100
		random.seed(3)
		seq = wf.EvolvableSequence(randomSequence(100,alphabet), base_fitness)
		tstart = time.time()
		Ne = 100
		pop = wf.WrightFisherPopulation(Ne, wf.SimpleMutator(mu,alphabet), wf.SequenceFitnessEvaluator())
		pop.populate(seq)
		self.assertTrue(pop.averageFitness()==base_fitness)
		mutseq = wf.EvolvableSequence(randomSequence(100,alphabet), 1.0)
		pop.inject(mutseq)
		self.assertTrue(pop.averageFitness()==((Ne-1)*base_fitness + 1.0)/Ne)

class test016(unittest.TestCase):
	"""reference counting check"""
	def test_refcount(self):
		alphabet = 'ATGC'
		mu = 0.0001
		base_fitness = 1.0
		n_gens = 100
		random.seed(3)
		tstart = time.time()
		Ne = 100
		pop = wf.WrightFisherPopulation(Ne, wf.SimpleMutator(mu,alphabet), wf.SequenceFitnessEvaluator())
		seq = wf.EvolvableSequence(randomSequence(100,alphabet), base_fitness)
		# Zero-fitness mutant
		mutseq = wf.EvolvableSequence(randomSequence(100,alphabet), 0.0)
		pop.populate(seq)
		mutentry = pop.inject(mutseq)
		# We put one in: should be one.
		self.assertTrue(mutentry.count==1)
		pop.evolve(1)
		# Evolution should result in immediate loss of this mutant from the population and, thus, the genebank.
		self.assertTrue(mutentry.count==0)
		self.assertTrue(pop.genebank.getEntry(mutentry.id)==None)

if __name__=="__main__":
	unittest.main(verbosity=2)
	#t = test008("test_run")
	#t.test_run()
