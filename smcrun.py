#!/software/bin/python
# Aylwyn Scally 2014

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

def prep(args):
	debug(args)
	# make dir
	if args.s1name == 'vcf':
		vcf_files = [os.path.abspath(x) for x in args.VCF_FILE]
		if args.samples:
			sname = '-'.join([os.path.splitext(os.path.basename(x))[0] for x in args.samples.split(',')])
		else:
			sname = os.path.splitext(args.VCF_FILE[0])[0]
	elif args.s1name == 'ms':
		ms_file = os.path.abspath(args.MS_FILE)
		sname = os.path.splitext(args.MS_FILE)[0]
	maindir = sname + '.smcdir'
	if not os.path.exists(maindir):
		info('creating directory %s' % maindir)
		if not args.sim:
			os.mkdir(maindir)
	if not args.sim and os.path.exists(maindir):
		os.chdir(maindir)

	# make segsep/psmcfa dir
	if args.psmc:
		subdir = 'psmcfa'
	else:
		subdir = 'segsep'
	if not os.path.exists(subdir):
		info('creating directory %s' % subdir)
		if not args.sim:
			os.mkdir(subdir)
	if not args.sim and os.path.exists(subdir):
		os.chdir(subdir)

	# make segsep/psmcfa files
	if args.s1name == 'vcf':
		for vcf in vcf_files:
			vcfname = os.path.basename(vcf).replace('.bgz', '').replace('.vcf', '')
			jobname = ':'.join((subdir, sname, vcfname))
			outname = '.'.join((vcfname, subdir))
			if args.samples:
				bcfview = 'bcftools view -s %s' % (args.samples)
			else:
				bcfview = 'bcftools view'
			if args.psmc:
				cmd = 'bsub.py "%s %s | vcfutils_noinfo.pl vcf2fq | fq2psmcfa -" -o psmcfa/$s -M 2 -j %s' % (bcfview, vcf, outname, jobname)
			else:
				cmd = 'bsub.py "%s %s | vcf-proc.py --segsep --alleles" -o %s -M 1 -t 2 -j %s' % (bcfview, vcf, outname, jobname)
			if args.bsim:
				cmd += ' --sim'
			info('submitting \'%s\'' % (jobname))
			aosutils.subcall(cmd, args.sim, wait = True)
	elif args.s1name == 'ms':
		jobname = ':'.join(('ms2psmc', sname))
		if args.phased:
			phasedarg = '--phased'
		else:
			phasedarg = ''
		if args.psmc:
			cmd = 'ms2psmc.py -l 1e7 %s --chrlen=1e8 --output=psmcfa %s' % (ms_file, phasedarg)
		else:
			cmd = 'ms2psmc.py -l 1e7 %s --chrlen=1e8 %s | awk \'{print > $1".segsep"}\'' % (ms_file, phasedarg)
		info('running \'%s\'' % (jobname))
		aosutils.subcall(cmd, args.sim, wait = True)

def run(args):
	# run smc inference
	if not args.sim:
		os.chdir(args.DIR)
	sname = os.path.splitext(os.path.basename(os.path.normpath(args.DIR)))[0]
	infiles = glob.glob('segsep/*.segsep')
	infiles.sort(key=aosutils.natural_key)
	if args.nfiles:
		infiles = infiles[:args.nfiles]
	debug(infiles)
	jobname = ':'.join(('smc', sname))
	sname = 'msmc'
	if args.psmc:
		if not args.memory:
			args.memory = 3
		#TODO: cat psmcfa files
		cmd = 'bsub.py "psmc -N25 -t15 -r5 -p \'4+25*2+4+6\' psmcfa/$s.psmcfa" -o %s.psmc -M %d -j %s' % (sname, sname, args.memory, jobname)
	else:
		if args.geneflow: # assume two samples
			if not args.memory:
				args.memory = 16
			cmd = 'bsub.py "msmc --fixedRecombination -P 0,0,1,1 -p 8*1+11*2 -t %d -o %s %s" -o %s.out -M %d -t %d -q long -j %s' % (args.threads, sname, ' '.join(infiles), sname, args.memory, args.threads, jobname)
		else:
			if not args.memory:
				args.memory = 10
			cmd = 'bsub.py "msmc -t %d -o %s segsep/*.segsep" -o %s.out -M %d -t %d -q long -j %s' % (args.threads, sname, sname, args.memory, args.threads, jobname)
	if args.bsim:
		cmd += ' --sim'
	info('submitting \'%s\'' % (jobname))
	aosutils.subcall(cmd, args.sim, wait = True)

def plot(args):
	# make plots
	if not args.sim:
		os.chdir(args.DIR)
	sname = os.path.splitext(os.path.basename(os.path.normpath(args.DIR)))[0]
	jobname = ':'.join(('smcplot', sname))

	if args.geneflow:
		cmd = 'smcplot.py -m run.final.txt -t %f -u %e --geneflow --hc -o %s-geneflow' % (args.tgen, args.mu, sname)
	else:
		cmd = 'smcplot.py -m run.final.txt -t %f -u %e --hc -o %s' % (args.tgen, args.mu, sname)
	info('running \'%s\'' % (jobname))
	aosutils.subcall(cmd, args.sim, wait = True)


pp = argparse.ArgumentParser(add_help=False)
pp.add_argument('--psmc', action='store_true', default = False, help = 'using psmc')
pp.add_argument('--sim', action='store_true', default = False, help = 'dry run')
pp.add_argument('-v', '--verbose', action='store_true', default = False)#, help = 'dry run')
pp.add_argument('--debug', action='store_true', default = False, help=argparse.SUPPRESS)
pp.add_argument('--bsim', action='store_true', default = False, help=argparse.SUPPRESS)

p = argparse.ArgumentParser()
s = p.add_subparsers()#help='sub-command help')

p1 = s.add_parser('prep', help='prepare files for psmc/msmc analysis')#, add_help=False)
s1 = p1.add_subparsers(dest='s1name')#help='sub-command help')
p11 = s1.add_parser('ms', parents=[pp])#, help='prep help')
p11.add_argument('MS_FILE', help = 'ms simulation file')
p11.add_argument('--phased', action='store_true', help='alleles are phased in input') 
p12 = s1.add_parser('vcf', parents=[pp])#, help='prep help')
p12.add_argument('VCF_FILE', nargs='+') 
p12.add_argument('-s', '--samples', help='comma-separated list of sample names in VCF_FILE') 
#p12.add_argument('--phased', action='store_true', help='alleles are phased in input') 
p1.set_defaults(func=prep)

p2 = s.add_parser('run', parents=[pp], help='run psmc/msmc')
p2.add_argument('DIR')
p2.add_argument('--geneflow', action='store_true', default = False, help = 'infer gene flow with msmc')
p2.add_argument('-n', '--nfiles', type=int, default=0, help = 'number of segsep files to include')
p2.add_argument('-t', '--threads', type=int, default=4, help = 'number of threads to use')
p2.add_argument('-M', '--memory', type=int, default=0, help = 'GB of RAM to use')
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

args.func(args)
