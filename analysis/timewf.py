import random, stats
import wrightfisher as wf

def randomSequence(n, alphabet):
	return ''.join(stats.sample_wr(alphabet, n))

if __name__=="__main__":
	alphabet = 'ATGC'
	dx = 0.1
	mu = 0.000001
	n_gens = 100000
	Ne = 1000
	random.seed(3)
	seq = wf.EvolvableSequence(randomSequence(100,alphabet))
	pop = wf.WrightFisherPopulation(Ne,wf.SimpleMutator(mu,alphabet))
	pop.populate(seq)
	for i in range(n_gens):
		pop.evolve(1)
		print pop.histogram()
