import sys, os, math, string, random, unittest
import stats
import wrightfisher as wf

class TraitOrganism(wf.Evolvable):
	def __init__(self):
		self.trait = 1.0
		self.key = randomSequence(12,string.letters)

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
			assert wf.pickIndexByProbability(cum_probs, x/float(n)) == x

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
		assert e == off

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
		assert dom_seq_entry == counts[0][1]

class test006(unittest.TestCase):
	"""Histogram"""
	def test_run(self):
		alphabet = 'ATGC'
		n = 1000
		pop = wf.Population(n,wf.SimpleMutator(0.0001,alphabet))
		seq = wf.EvolvableSequence(randomSequence(90,'ATGC'))
		pop.populate(seq)
		h = pop.histogram()
		assert h[0] == n
		for i in range(10):
			pop.evolve(1)
			h = pop.histogram()
			assert sum(h) == n

class test007(unittest.TestCase):
	"""Members iterator"""
	def test_run(self):
		alphabet = 'ATGC'
		n = 1000
		pop = wf.Population(n,wf.SimpleMutator(0.0001,alphabet))
		seq = wf.EvolvableSequence(randomSequence(90,'ATGC'))
		pop.populate(seq)
		

if __name__=="__main__":
	unittest.main(verbosity=2)
