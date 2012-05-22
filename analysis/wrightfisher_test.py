import sys, os, math, string, random, unittest, time
import stats
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
		pop = wf.WrightFisherPopulation(N, wf.SimpleMutator(0.01,alphabet))
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
		p = wf.WrightFisherPopulation(10, wf.SimpleMutator(0.01,alphabet))
		p.populate(TraitOrganism())
		p.evolve(1000)
		#print len(p.members)
		#for o in p.members:
		#	print o.trait

class test003(unittest.TestCase):
	"""Sequence"""
	def test_run(self):
		alphabet = 'ATGC'
		p = wf.WrightFisherPopulation(10, wf.SimpleMutator(0.01,alphabet))
		p.populate(wf.EvolvableSequence(randomSequence(20,alphabet)))
		p.evolve(1000)
		#for m in p.members:
		#	print m

class test004(unittest.TestCase):
	def test_run(self):
		alphabet = 'ATGC'
		e = wf.EvolvableSequence(randomSequence(100,alphabet))
		# Zero mutation rate copy
		off = e.spawn(wf.SimpleMutator(0.0)).offspring
		# Should yield identical sequences
		self.assertTrue(e == off)

class test005(unittest.TestCase):
	"""Tracking frequency"""
	def test_run(self):
		alphabet = 'ATGC'
		pop = wf.WrightFisherPopulation(100,wf.SimpleMutator(0.0001,alphabet))
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
	def test_run(self):
		alphabet = 'ATGC'
		n = 1000
		pop = wf.WrightFisherPopulation(n,wf.SimpleMutator(0.0001,alphabet))
		seq = wf.EvolvableSequence(randomSequence(90,alphabet))
		pop.populate(seq)
		h = pop.histogram()
		self.assertTrue(h[0] == n)
		for i in range(10):
			pop.evolve(1)
			h = pop.histogram()
			self.assertTrue(sum(h) == n)

class test007(unittest.TestCase):
	"""Mutation info"""
	def test_run(self):
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
	def test_run(self):
		alphabet = 'ATGC'
		dx = 0.1
		Ne = 20
		mu = 0.001
		pop = wf.WrightFisherPopulation(Ne,wf.SimpleMutator(mu,alphabet))
		seq = wf.EvolvableSequence(randomSequence(100,alphabet))
		random.seed(3)
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
	def test_run(self):
		alphabet = 'ATGC'
		dx = 0.1
		Ne = 1000
		mu = 0.0001
		predicted_fixation_probability = wf.probabilityOfFixation(Ne, dx)
		n_fixations = 0
		n_total = 0
		n_trials = 5*int(1/dx)
		random.seed(111)
		print "Aborting because this test takes too long"
		return
		for i in range(n_trials):
			pop = wf.WrightFisherPopulation(Ne,wf.SimpleMutator(mu,alphabet))
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
			#print n_total, n_trials, n_fixations, res.fixed, res.time, n_fixations/float(n_total), predicted_fixation_probability
		# DAD: mean is 
		est_prob = n_fixations/float(n_total)
		p = predicted_fixation_probability
		exp_fixations = n_trials*p
		sd = math.sqrt(n_trials*p*(1.0-p))
		# 
		self.assertTrue(n_fixations <= (exp_fixations+2*sd))
		self.assertTrue(n_fixations >= (exp_fixations-2*sd))
		#print n_fixations, n_total, n_fixations/float(n_total), predicted_fixation_probability
			
class test010(unittest.TestCase):
	"""LCA"""
	def test_lca(self):
		alphabet = 'ATGC'
		pop = wf.WrightFisherPopulation(100,wf.SimpleMutator(0.0001,alphabet))
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
		pop = wf.WrightFisherPopulation(100,wf.SimpleMutator(0.0001,alphabet))
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
		pop = wf.WrightFisherPopulation(Ne,wf.SimpleMutator(mu,alphabet))
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
		mu = 0.0001
		n_gens = 100
		random.seed(3)
		seq = wf.EvolvableSequence(randomSequence(100,alphabet))
		for i in range(4):
			tstart = time.time()
			Ne = 10**i
			pop = wf.WrightFisherPopulation(Ne,wf.SimpleMutator(mu,alphabet))
			pop.populate(seq)
			pop.evolve(n_gens)
			tend = time.time()
			print "# evolved {} generations at Ne={} (t={}ms)".format(n_gens, Ne, (tend-tstart)*1000)
	
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

if __name__=="__main__":
	unittest.main(verbosity=2)
	#t = test008("test_run")
	#t.test_run()
