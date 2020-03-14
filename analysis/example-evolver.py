#! python

import sys, os, math, random, argparse, itertools, time
import scipy as sp
import util
import wrightfisher as wf


class MySequenceFitnessEvaluator(wf.SequenceFitnessEvaluator):
	def __init__(self):
		pass

	"""Interface for fitness evaluation."""
	def fitness(self, organism, population):
		seq = organism.sequence # Relies on 
		# Fitness is the fraction of E and K plus some constant small factor
		# so that sequences without E and K don't simply go extinct
		epsilon = 0.001
		fit = (seq.count('E') + seq.count('K') + epsilon)/float(len(seq))
		return fit



if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Example use of Wright-Fisher evolution")
	# Required arguments
	parser.add_argument(dest="population_size", type=int, help="size of evolving population")
	parser.add_argument(dest="num_generations", type=int, help="number of generations to simulate")
	# Optional arguments
	parser.add_argument("--output-frequency", dest="output_frequency", type=int, default=1, help="frequency of output, in generations")
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

	'''
	# Read input
	if not os.path.isfile(options.in_fname):
		raise IOError("# Error: file {} does not exist".format(options.in_fname))
	with open(options.in_fname,'r') as inf:
		# Read a FASTA file?
		(headers, seqs) = biofile.readFASTA(inf)
	'''

	def randomSequence(n, alphabet):
		indices = range(len(alphabet))
		return ''.join([alphabet[i] for i in sp.random.choice(indices, size=n, replace=True)])


	# Set up parameters of the evolving population
	alphabet = 'ACDEFGHIKLMNPQRSTVWY'
	mutation_rate = 0.0001
	base_fitness = 1.0
	sp.random.seed(3)
	seq = wf.EvolvableSequence(randomSequence(50,alphabet), base_fitness)
	tstart = time.time()
	mutator = wf.SimpleMutator(mutation_rate,alphabet)
	fitness_evaluator = MySequenceFitnessEvaluator()

	# Create the population
	pop = wf.WrightFisherPopulation(options.population_size, mutator, fitness_evaluator)
	pop.populate(seq)

	# Write output
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
	total_iterations = int(options.num_generations/float(generations_per_output))

	# Write out starting time
	data_outs.write("# Evolution finished {}\n".format(util.timestamp()))

	for n in range(total_iterations):
		# Evolve the population by 
		pop.evolve(generations_per_output)

		# A dictionary of results, one result per addHeader call above
		result = dout.createResult(default=None)
		result['generation'] = pop.generations
		result['average.fitness'] = pop.averageFitness()
		lca_entry = pop.lastCommonAncestor()
		result['lca'] = lca_entry.organism.sequence
		# Parse the values, convert Nones to NA, etc.
		line = dout.formatLine(result)
		data_outs.write(line)
		n_written += 1
	tstop = time.time()

	# Write out stopping time
	data_outs.write("# Evolution finished {stamp} ({elapsed:f} seconds elapsed)\n".format(stamp=util.timestamp(), elapsed=tstop-tstart))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

