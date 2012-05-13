import sys, os, math, string, random, unittest
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
			wf.mutated = True
		else:
			offs.trait = self.trait
		return res

class test001(unittest.TestCase):
	"""basic run"""
	def test_run(self):
		print "here"
		p = wf.Population(10, Organism, 0.01)
		p.evolve(1000)
		for o in p.members:
			print o.trait

if __name__=="__main__":
	unittest.main(verbosity=2)
