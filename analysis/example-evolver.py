#! python

import sys, os, math, random, argparse, itertools, time
import scipy as sp
import util
import wrightfisher as wf

# A more detailed example of how the Wright-Fisher simulation can be used
# Run with:
# 	python3 example-evolver.py --help
# 	python3 example-evolver.py 100 10000 --output-frequency 10



class MySequenceFitnessEvaluator(wf.SequenceFitnessEvaluator):
	"""A simple example of how to create a custom fitness function for evolving organisms which have a sequence."""
	def fitness(self, organism, population):
		seq = organism.sequence # We can do this if we're using EvolvableSequences
		# A simple definition of sequence fitness based on simple sequence properties.
		# This fitness definition will lead to hill-climbing toward sequences consisting of all E and K.
		# 
		# Fitness is the fraction of E and K plus some constant small factor
		# so that sequences without E and K don't simply go extinct
		epsilon = 0.001
		fit = (seq.count('E') + seq.count('K') + epsilon)/float(len(seq))
		return fit

class MyCachedSequenceFitnessEvaluator(wf.SequenceFitnessEvaluator):

	def createCache(self):
		self._cache = {}

	"""A simple example of how to create a custom fitness function for evolving organisms which have a sequence."""
	def fitness(self, organism, population):
		seq = organism.sequence # We can do this if we're using EvolvableSequences
		# A simple definition of sequence fitness based on simple sequence properties.
		# This fitness definition will lead to hill-climbing toward sequences consisting of all E and K.
		# 
		# Fitness is the fraction of E and K plus some constant small factor
		# so that sequences without E and K don't simply go extinct

		# Cache the fitness once we've computed it, and only recompute if necessary.
		fit = None
		try:
			fit = self._cache[seq]
		except KeyError:
			epsilon = 0.001
			fit = (seq.count('E') + seq.count('K') + epsilon)/float(len(seq))
			self._cache[seq] = fit
		return fit

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Example use of Wright-Fisher evolution")
	# Required arguments
	parser.add_argument(dest="population_size", type=int, help="size of evolving population")
	parser.add_argument(dest="num_generations", type=int, help="number of generations to simulate")
	# Optional arguments
	parser.add_argument("--output-frequency", dest="output_frequency", type=int, default=1, help="frequency of output, in generations")
	parser.add_argument("--random-seed", dest="random_seed", type=int, default=None, help="seed for random number generator")
	parser.add_argument("-o", "--out", dest="out_fname", default=None, help="output filename")
	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()

	# Start up output
	if not options.out_fname is None:
		outf = open(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		# By default, write to stdout
		data_outs.addStream(sys.stdout)

	# Write out parameters
	data_outs.write("# Run started {}\n".format(util.timestamp()))
	data_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	def randomSequence(n, alphabet):
		""" A simple random sequence generator. """
		indices = range(len(alphabet))
		return ''.join([alphabet[i] for i in sp.random.choice(indices, size=n, replace=True)])

	# Set up parameters of the evolving population
	alphabet = 'ACDEFGHIKLMNPQRSTVWY'
	mutation_rate = 0.0001
	base_fitness = 1.0
	# Random seed
	if not options.random_seed is None:
		sp.random.seed(options.random_seed)

	seq = wf.EvolvableSequence(randomSequence(50,alphabet), base_fitness)
	# Start recording the time
	tstart = time.time()
	mutator = wf.SimpleMutator(mutation_rate,alphabet)
	fitness_evaluator = MyCachedSequenceFitnessEvaluator()
	fitness_evaluator.createCache()

	# Create the population
	pop = wf.WrightFisherPopulation(options.population_size, mutator, fitness_evaluator)
	pop.populate(seq)

	# Write output
	# This uses the util.DelimitedOutput() class, which produces self-documenting tab-delimited output
	dout = util.DelimitedOutput()
	dout.addHeader('generation','Generation (1-based)','d')
	dout.addHeader('average.fitness','Average fitness of population','f')
	dout.addHeader('lca','Last common ancestor of the population','s')
	# Write the header descriptions
	dout.describeHeader(data_outs)
	# Write the header fields
	dout.writeHeader(data_outs)
	n_written = 0

	generations_per_output = options.output_frequency
	# Total number of reporting iterations to run
	# The +1 is because the first iteration will write out the generation-zero starting values
	total_iterations = int(options.num_generations/float(generations_per_output))+1

	# Write out starting time
	data_outs.write("# Evolution finished {}\n".format(util.timestamp()))

	for n in range(total_iterations):
		# Create a dictionary of results
		result = dout.createResult(default=None)
		result['generation'] = pop.generations
		result['average.fitness'] = pop.averageFitness()
		lca_entry = pop.lastCommonAncestor()
		result['lca'] = lca_entry.organism.sequence
		# Parse the values, convert Nones to NA, etc.
		line = dout.formatLine(result)
		data_outs.write(line)
		n_written += 1
		# Evolve the population
		pop.evolve(generations_per_output)
	tstop = time.time()

	# Write out stopping time
	data_outs.write("# Evolution finished {stamp} ({elapsed:f} seconds elapsed)\n".format(stamp=util.timestamp(), elapsed=tstop-tstart))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

