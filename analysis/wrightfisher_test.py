import sys, os, math, string, random, unittest
import stats
import wrightfisher as wf

class Organism(wf.Evolvable):
	def __init__(self):
		self.trait = 1.0
	def fitness(self):
		return self.trait
	def spawn(self, mutation_rate):
		res = wf.SpawnResult()
		p = random.random()
		offs = Organism()
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
		p = wf.Population(0.01)
		for n in range(10):
			p.addMember(Organism())
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
		p = wf.Population(0.01)
		alphabet = 'ATGC'
		for n in range(10):
			p.addMember(wf.EvolvableSequence(randomSequence(20,alphabet),alphabet))
		p.evolve(1000)
		#for m in p.members:
		#	print m

class test004(unittest.TestCase):
	def test_run(self):
		e = wf.EvolvableSequence(randomSequence(100,'ATGC'),'ATGC')
		off = e.spawn(0.0).offspring
		assert e.eq(off)

if __name__=="__main__":
	unittest.main(verbosity=2)
