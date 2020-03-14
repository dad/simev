#! python

import sys, os, math, random, argparse, itertools
import scipy as sp
import util
import wrightfisher as wf

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Example use of Wright-Fisher evolution")
	# Required arguments
	parser.add_argument(dest="population_size", type=int, help="size of evolving population")
	parser.add_argument(dest="num_generations", type=int, help="number of generations to simulate")
	# Optional arguments
	#parser.add_argument("--randomize", dest="randomize", action="store_true", default=False, help="randomize each sequence?")
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

	# Set up parameters of the evolving population
	alphabet = 'ACDEFGHIKLMNPQRSTVWY'
	mutation_rate = 0.0001
	base_fitness = 1.0
	sp.random.seed(3)
	seq = wf.EvolvableSequence(randomSequence(100,alphabet), base_fitness)
	tstart = time.time()
	mutator = wf.SimpleMutator(mu,alphabet)
	fitness_evaluator = wf.SequenceFitnessEvaluator()

	# Create the population
	pop = wf.WrightFisherPopulation(options.population_size, mutator, fitness_evaluator)
	pop.populate(seq)

	# Write output
	dout = util.DelimitedOutput()
	dout.addHeader('from.amino.acid','Starting amino acid one-letter code','s')
	# Write the header descriptions
	dout.describeHeader(data_outs)
	# Write the header fields
	dout.writeHeader(data_outs)
	n_written = 0


	# Write out starting time
	data_outs.write("# Evolution finished {}\n".format(util.timestamp()))

	for n in range(options.num_generations):
		# A dictionary of results, one result per addHeader call above
		result = dout.createResult(default=None)
		pop.evolve(1)
		result['generation'] = n+1
		result['average.fitness'] = pop.averageFitness()

		# Parse the values, convert Nones to NA, etc.
		line = dout.formatLine(result)
		# A more manual approach:
		# line = format.format(column1="the answer is", column2=42)
		data_outs.write(line)
		n_written += 1


	# Write out stopping time
	data_outs.write("# Evolution finished {stamp} ({elapsed:f} elapsed time)\n".format(stamp=util.timestamp(), elapsed=tstop-tstart))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()

