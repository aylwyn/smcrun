#!/usr/bin/env python
# Aylwyn Scally 2015

import sys
import argparse
import os
import re
import os.path
import glob
#from numpy import median
import logging
from logging import error, warning, info, debug, critical

import aosutils

# TODO: start from bam files

def pairs(plist):
	for i, x in enumerate(plist[:-1]):
		for y in plist[i+1:]:
			yield((x, y))

def prep(args):
	debug(args)
	callmask = 'callmask.bed.gz'

	class Sample(object):
		def __init__(self, names, masks = []):
			self.names = names
			self.masks = masks
		
		def __add__(self, other):
			return(Sample(self.names + other.names, self.masks + other.masks))


	if args.s1name == 'vcf':
		vcf_input = [(os.path.abspath(x), '') for x in args.VCF_FILE]
		if args.vcf_list:
			for line in open(args.vcf_list):
				vcf_input.append(tuple(line.split()))

		inputsamps = []
		if args.samples:
			inputsamps = [Sample([x]) for x in args.samples.split(',')]
		if args.sample_list:
			for tok in (line.split() for line in open(args.sample_list)):
				inputsamps.append(Sample(tok[0:1], tok[1:]))
		if not inputsamps:
			inputsamps = [Sample(['*'])]

		runsamps = []
		if args.allpairs:
			runsamps = [pair[0] + pair[1] for pair in pairs(inputsamps)]
		else:
			runsamps = [inputsamps[0]]
			for samp in inputsamps[1:]:
				runsamps[0] += samp

		if args.pseudodip:
			pdip_arg = '--pseudodip'
		else:
			pdip_arg = ''

		if args.callmask:
			for samp in runsamps:
				samp.masks.append(os.path.abspath(args.callmask))

	elif args.s1name == 'ms':
		ms_file = os.path.abspath(args.MS_FILE)
		if not os.path.exists(ms_file):
			error('No such file: %s' % ms_file)
			return(2)
		runsamps = [Sample([ms_file])]

	debug('\n'.join([','.join(x.names) for x in runsamps]))
	debug('\n'.join([','.join(x.masks) for x in runsamps]))
	rdir = os.path.abspath('.')
	for samp in runsamps:
		os.chdir(rdir)
		if args.s1name == 'vcf':
			sname = '-'.join([os.path.splitext(os.path.basename(x))[0] for x in samp.names])
		elif args.s1name == 'ms':
			sname = os.path.basename(samp.names[0])
			if args.unphased:
				sname += '.unphased'

		maindir = sname + '.smcdir'
		if not os.path.exists(maindir):
			info('creating directory %s' % maindir)
			if not args.sim:
				os.mkdir(maindir)
		if os.path.exists(maindir):
			os.chdir(maindir)

		# make segsep/psmcfa dir
		if args.psmc:
			subdir = 'psmcfa'
		else:
			subdir = 'segsep'
		if not os.path.exists(subdir):
			info('creating directory %s/%s' % (maindir, subdir))
			if not args.sim:
				os.mkdir(subdir)
		elif not args.replace:
			warning('directory %s/%s exists; skipping; use --replace to override' % (maindir, subdir))
#			return(2)
			continue
		if os.path.exists(subdir):
			os.chdir(subdir)

		# mask argument (make merged mask if necessary)
		if len(samp.masks) > 1:
			mrgcmd = 'bedtools merge %s | gzip > %s' % (' '.join(['-i ' + x for x in samp.masks]), callmask)
			info('merging callmask files for %s' % (sname))
			aosutils.subcall(mrgcmd, args.sim, wait = True)
			cmask_arg = '--callmask=%s' % callmask
		elif len(samp.masks) == 1:
			cmask_arg = '--callmask=%s' % samp.masks[0]
		else:
			cmask_arg = ''

		# make segsep/psmcfa files
		if args.s1name == 'vcf':
			if not args.memory:
				args.memory = 3

			if args.psmc:
				outarg = '--psmcfa'
			else:
				outarg = '--segsep'

			if args.concat or args.merge:
				if args.concat:
					vcfname = 'all_chrs'
					bcfcmd = 'bcftools concat'
				elif args.merge:
					vcfname = 'all_samps'
					bcfcmd = 'bcftools merge'
				jobname = ':'.join((subdir, sname, vcfname))
				outname = '.'.join((vcfname, subdir))
				vcfs =  ' '.join([x[0] for x in vcf_input])
				jobname = ':'.join((subdir, sname))
				if samp.names[0] == '*':
					bcfview = '%s %s' % (bcfcmd, vcfs)
				else:
					bcfview = '%s %s | bcftools view -s %s - ' % (bcfcmd, vcfs, ','.join(samp.names))
				cmd = '%s "%s | vcf-proc.py %s %s %s" -o %s -M %d -j %s' % (args.submit, bcfview, outarg, pdip_arg, cmask_arg, outname, args.memory, jobname)
				if args.replace:
					cmd += ' --replace'
				if args.bsim:
					cmd += ' --sim'
				info('submitting \'%s\'' % (jobname))
				aosutils.subcall(cmd, args.sim, wait = True)
			else:
				if args.vcf_list and not (args.sim or os.path.exists('tmp')):
					os.mkdir('tmp')

				for vcf, repvcf in vcf_input:
					vcfname = os.path.basename(vcf).replace('.bgz', '').replace('.vcf', '')
					jobname = ':'.join((subdir, sname, vcfname))
					outname = '.'.join((vcfname, subdir))
					if samp.names[0] == '*':
						bcfview = 'bcftools view'
					else:
						bcfview = 'bcftools view -s %s' % (','.join(samp.names))
					if repvcf:
						tmprep = 'tmp/%s.rep.gz' % vcfname
						if not args.usetmp:
							cmd = '%s %s | gzip > %s' % (bcfview, repvcf, tmprep)
							info('extracting replacement calls from %s' % (repvcf))
							aosutils.subcall(cmd, args.sim, wait = True)
						reparg = '--replacecalls=%s' % tmprep
					else:
						reparg = ''
					cmd = '%s "%s %s | vcf-proc.py %s %s %s %s" -o %s -M %d -j %s' % (args.submit, bcfview, vcf, outarg, pdip_arg, cmask_arg, reparg, outname, args.memory, jobname)
					if args.replace:
						cmd += ' --replace'
					if args.bsim:
						cmd += ' --sim'
					info('submitting \'%s\'' % (jobname))
					aosutils.subcall(cmd, args.sim, wait = True)
		elif args.s1name == 'ms':
			jobname = ':'.join(('ms2smc', sname))
			if args.unphased:
				phasedarg = '--unphased'
			else:
				phasedarg = ''
			if args.psmc:
				cmd = 'ms2smc.py -l 1e7 %s --chrlen=%d --output=psmcfa %s' % (samp.names[0], args.chrlen, phasedarg)
			else:
				cmd = 'ms2smc.py -l 1e7 %s --chrlen=%d %s | awk \'{print > $1".segsep"}\'' % (samp.names[0], args.chrlen, phasedarg)
			info('running \'%s\'' % (jobname))
			aosutils.subcall(cmd, args.sim, wait = True)


def run(args): # run smc inference
	os.chdir(args.DIR)
	debug('In %s:' % args.DIR)
	if args.psmc:
		psmcfaf = glob.glob('psmcfa/*.psmcfa')[0]
		sname = 'psmc'
		jobname = ':'.join(('psmc', args.DIR))
		outf = '%s.out' % sname
		if not args.memory:
			args.memory = 3
		if args.decode:
			decarg = '-d'
		else:
			decarg = ''
		#TODO: cat args.nfiles psmcfa files into all.psmfca
		cmd = '%s "psmc -p \'%s\' %s %s" -o %s -M %d -j %s' % (args.submit, args.intervals, decarg, psmcfaf, outf, args.memory, jobname)
#		cmd = 'bsub.py "psmc -N25 -t15 -r5 -p \'%s\' %s" -o %s -M %d -j %s' % (args.intervals, psmcfaf, outf, args.memory, jobname)
		if args.replace:
			cmd += ' --replace'
		if args.bsim:
			cmd += ' --sim'
		if os.path.exists(outf) and not args.replace:
			warning('%s exists in %s; use --replace' % (outf, args.DIR))
		else:
			info('submitting \'%s\'' % (jobname))
			aosutils.subcall(cmd, args.sim, wait = True)
	else:
		infiles = glob.glob('segsep/*.segsep')
		infiles.sort(key=aosutils.natural_key)
		if args.nfiles:
			infiles = infiles[:args.nfiles]
		debug(infiles)

		if args.combined:
			infiles = [' '.join(infiles)]

		for jf, infile in enumerate(infiles):
			if args.combined:
				pref = 'combined_%d' % args.nfiles
			else:
				pref = os.path.basename(infile).replace('.segsep', '')
			sname = '%s.msmc' % pref
			outf = sname + '.out'
			jobname = ':'.join(('msmc', os.path.basename(os.path.normpath(args.DIR)), pref))
			if args.geneflow: # assume two samples
				if not args.memory:
					args.memory = 16
				cmd = '%s "msmc --fixedRecombination -P 0,0,1,1 -p \'%s\' -t %d -o %s %s" -o %s -M %d -t %d -q %s -j %s' % (args.submit, args.intervals, args.threads, sname, infile, outf, args.memory, args.threads, args.queue, jobname)
			else:
				if not args.memory:
					args.memory = 10
				cmd = '%s "msmc -t %d -o %s %s" -o %s -M %d -t %d -q %s -j %s' % (args.submit, args.threads, sname, outf, infile, args.memory, args.threads, args.queue, jobname)
			if args.replace:
				cmd += ' --replace'
			if args.bsim:
				cmd += ' --sim'
			if os.path.exists(outf) and not args.replace:
				warning('%s/%s exists; use --replace' % (args.DIR, outf))
			else:
				info('submitting \'%s\'' % (jobname))
				aosutils.subcall(cmd, args.sim, wait = True)

def plot(args): # make plots
	if not args.sim:
		os.chdir(args.DIR)
#	sname = os.path.splitext(os.path.basename(os.path.normpath(args.DIR)))[0]
	sname = os.path.basename(os.path.normpath(args.DIR))
	jobname = ':'.join(('smcplot', sname))

	if args.geneflow:
		cmd = 'smcplot.py -m run.final.txt -t %f -u %e --geneflow --hc -o %s-geneflow' % (args.tgen, args.mu, sname)
	else:
		cmd = 'smcplot.py -m run.final.txt -t %f -u %e --hc -o %s' % (args.tgen, args.mu, sname)
	info('running \'%s\'' % (jobname))
	aosutils.subcall(cmd, args.sim, wait = True)


pp = argparse.ArgumentParser(add_help=False)
pp.add_argument('--psmc', action='store_true', default = False, help = 'using psmc')
pp.add_argument('--bsub', action='store_true', default = False, help = 'submit jobs via bsub')
pp.add_argument('--replace', action='store_true', default = False, help = 'replace existing files')
pp.add_argument('--sim', action='store_true', default = False, help = 'dry run')
pp.add_argument('-v', '--verbose', action='store_true', default = False)#, help = 'dry run')
pp.add_argument('--debug', action='store_true', default = False, help=argparse.SUPPRESS)
pp.add_argument('--bsim', action='store_true', default = False, help=argparse.SUPPRESS)
pp.add_argument('-M', '--memory', type=int, default=0, help = 'GB of RAM to use')

p = argparse.ArgumentParser()
s = p.add_subparsers()#help='sub-command help')

p1 = s.add_parser('prep', help='prepare files for psmc/msmc analysis')#, add_help=False)
s1 = p1.add_subparsers(dest='s1name')#help='sub-command help')
p11 = s1.add_parser('ms', parents=[pp])#, help='prep help')
p11.add_argument('MS_FILE', help = 'ms simulation file')
p11.add_argument('--chrlen', default = 50e6, help='alleles are considered unphased in input') 
p11.add_argument('--unphased', action='store_true', default = False, help='alleles are considered unphased in input') 
p12 = s1.add_parser('vcf', parents=[pp])#, help='prep help')
p12.add_argument('VCF_FILE', nargs='*') 
p12.add_argument('-s', '--samples', help='comma-separated list of sample names in VCF_FILE') 
p12.add_argument('-S', '--sample_list', help='file containing list of sample names and (optionally) paths to mask files of uncallable regions') 
p12.add_argument('--allpairs', action='store_true', default = False, help='run on all pairs of sample names in list') 
p12.add_argument('--pseudodip', action='store_true', default = False, help='call vcf-proc.py with --pseudodip flag') 
p12.add_argument('--concat', action='store_true', default = False, help='concatenate input vcf files into single output (i.e. combine chromosomes; ignores replacement calls)') 
p12.add_argument('--merge', action='store_true', default = False, help='merge input vcf files into single output (i.e. combine samples; ignores replacement calls)') 
p12.add_argument('--callmask', default = '', help='call vcf-proc.py with --callmask=CALLMASK of uncallable regions applied to all samples (Note: this is additional to any masks specified via a SAMPLE_LIST') 
p12.add_argument('-f', '--vcf_list', help='file containing a list of input vcfs (one per line). For each one, a vcf of replacement calls may specified in a second column.') 
p12.add_argument('--usetmp', action='store_true', default=False, help='use existing replacement calls in tmp dir') 
p1.set_defaults(func=prep)

p2 = s.add_parser('run', parents=[pp], help='run psmc/msmc')
p2.add_argument('DIR')
p2.add_argument('--geneflow', action='store_true', default = False, help = 'infer gene flow with msmc')
p2.add_argument('--decode', action='store_true', default = False, help = 'decode Ne (psmc)')
p2.add_argument('-p', '--intervals', default = '6*1+7*2', help = 'time interval specification (psmc format); e.g. "4+25*2+4+6", "6*1+7*2"')
p2.add_argument('--combined', action='store_true', default = False, help = 'run on all segsep files combined')
p2.add_argument('-n', '--nfiles', type=int, default=0, help = 'number of segsep files to include')
p2.add_argument('-t', '--threads', type=int, default=4, help = 'number of threads to use')
p2.add_argument('-q', '--queue', default='normal', help = 'queue to use')
p2.set_defaults(func=run)

p3 = s.add_parser('plot', parents=[pp], help='plot psmc/msmc results')
p3.add_argument('DIR')
p3.add_argument('--geneflow', action='store_true', default = False, help = 'infer gene flow with msmc')
p3.add_argument('-t', '--tgen', type=int, default=30.0, help = 'generation time (y)')
p3.add_argument('-u', '--mu', type=int, default=1.25e-8, help = 'per-generation mutation rate (bp^-1)')
p3.set_defaults(func=plot)

args = p.parse_args()

loglevel = logging.WARNING
if args.verbose:
	loglevel = logging.INFO
if args.debug:
	loglevel = logging.DEBUG
logging.basicConfig(format = '%(module)s:%(lineno)d:%(levelname)s: %(message)s', level = loglevel)

if args.bsub:
	args.submit = 'submit.py bsub'
else:
	args.submit = 'submit.py nohup'

args.func(args)
