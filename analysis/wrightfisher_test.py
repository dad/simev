import sys, os, math, string, random, unittest
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
	"""basic run"""
	def test_run(self):
		alphabet='ATGC'
		p = wf.Population(10, wf.SimpleMutator(0.01,alphabet))
		p.populate(TraitOrganism())
		p.evolve(1000)
		#print len(p.members)
		#for o in p.members:
		#	print o.trait

class test002(unittest.TestCase):
	"""pickByIndex"""
	def test_run(self):
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

class test003(unittest.TestCase):
	"""Sequence"""
	def test_run(self):
		alphabet = 'ATGC'
		p = wf.Population(10, wf.SimpleMutator(0.01,alphabet))
		p.populate(wf.EvolvableSequence(randomSequence(20,alphabet)))
		p.evolve(1000)
		#for m in p.members:
		#	print m

class test004(unittest.TestCase):
	def test_run(self):
		e = wf.EvolvableSequence(randomSequence(100,'ATGC'))
		# Zero mutation rate copy
		off = e.spawn(wf.SimpleMutator(0.0)).offspring
		# Should yield identical sequences
		self.assertTrue(e == off)

class test005(unittest.TestCase):
	"""Tracking frequency"""
	def test_run(self):
		alphabet = 'ATGC'
		pop = wf.Population(100,wf.SimpleMutator(0.0001,alphabet))
		seq = wf.EvolvableSequence(randomSequence(100,'ATGC'))
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
		pop = wf.Population(n,wf.SimpleMutator(0.0001,alphabet))
		seq = wf.EvolvableSequence(randomSequence(90,'ATGC'))
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
		seq = wf.EvolvableSequence(randomSequence(50,'ATGC'))
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
		mu = 0.0001
		predicted_fixation_probability = wf.probabilityOfFixation(Ne, dx)
		n_fixations = 0
		n_total = 0
		for i in range(5*int(1/dx)):
			pop = wf.Population(Ne,wf.SimpleMutator(0.0000001,alphabet))
			seq = wf.EvolvableSequence(randomSequence(100,'ATGC'))
			# No fitness
			seq.fitness = 0.0
			pop.populate(seq)
			mutseq = wf.EvolvableSequence(randomSequence(100,'ATGC'))
			mutseq.fitness = 1.0
			mutentry = pop.inject(mutseq)
			print ''
			for m in pop.members:
				print m.fitness
			res = pop.evolveUntilFixationOrLossOf(mutentry)
			# Everyone else has zero fitness -- fixation is assured.
			self.assertTrue(res.fixed)
		
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
		for i in range(5*int(1/dx)):
			pop = wf.Population(Ne,wf.SimpleMutator(0.0000001,alphabet))
			seq = wf.EvolvableSequence(randomSequence(100,'ATGC'))
			seq.fitness = 1.0
			pop.populate(seq)
			mutseq = wf.EvolvableSequence(randomSequence(100,'ATGC'))
			mutseq.fitness = seq.fitness + dx
			mutentry = pop.inject(mutseq)
			res = pop.evolveUntilFixationOrLossOf(mutentry)
			n_total += 1
			if res.fixed:
				n_fixations += 1
		print n_fixations, n_total, n_fixations/float(n_total), predicted_fixation_probability
			
		

if __name__=="__main__":
	unittest.main(verbosity=2)
