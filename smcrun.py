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
		if args.samples:
			sname = '-'.join([os.path.splitext(os.path.basename(x))[0] for x in args.samples.split(',')])
		else:
			sname = os.path.splitext(args.VCF_FILE[0])[0]
	elif args.s1name == 'ms':
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
		for vcf in args.VCF_FILE:
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
		if args.psmc:
			cmd = 'ms2psmc.py -l 1e7 %s --chrlen=1e8 --output=psmcfa' % (sname)
		else:
			cmd = 'ms2psmc.py -l 1e7 %s --chrlen=1e8 | awk \'{print > $1.segsep}\'' % (sname)
		info('submitting \'%s\'' % (jobname))
		aosutils.subcall(cmd, args.sim, wait = True)

def run(args):
	# run smc inference
	if not args.sim:
		os.chdir(args.DIR)
	sname = os.path.splitext(os.path.basename(args.DIR))[0]
	infiles = glob.glob('segsep/*.segsep')
	infiles.sort(key=aosutils.natural_key)
	if args.n:
		infiles = infiles[:args.n]
	debug(infiles)
	if args.psmc:
		jobname = ':'.join(('psmc', sname))
		#TODO: cat psmcfa files
		cmd = 'bsub.py "psmc -N25 -t15 -r5 -p \'4+25*2+4+6\' psmcfa/$s.psmcfa" -o %s.psmc -M 3 -j %s' % (sname, sname, jobname)
	else:
		jobname = ':'.join(('msmc', sname))
		if args.geneflow: # assume two samples
			cmd = 'bsub.py "msmc --fixedRecombination -P 0,0,1,1 -p 8*1+11*2 -t 8 -o %s %s" -o %s.msmc.out -M 20 -t 8 -q long -j %s' % (sname, ' '.join(infiles), sname, jobname)
		else:
			cmd = 'bsub.py "msmc -t 4 -o %s segsep/*.segsep" -o %s.msmc.out -M 10 -t 4 -q long -j %s' % (sname, sname, jobname)
	if args.bsim:
		cmd += ' --sim'
	info('submitting \'%s\'' % (jobname))
	aosutils.subcall(cmd, args.sim, wait = True)

	# make plots


pp = argparse.ArgumentParser(add_help=False)
pp.add_argument('--psmc', action='store_true', default = False, help = 'using psmc')
pp.add_argument('--sim', action='store_true', default = False, help = 'dry run')
pp.add_argument('-v', '--verbose', action='store_true', default = False)#, help = 'dry run')
pp.add_argument('--debug', action='store_true', default = False, help=argparse.SUPPRESS)
pp.add_argument('--bsim', action='store_true', default = False, help=argparse.SUPPRESS)

p = argparse.ArgumentParser()
s = p.add_subparsers()#help='sub-command help')

p1 = s.add_parser('prep', help='prepare files for smc analysis')
s1 = p1.add_subparsers(dest='s1name')#help='sub-command help')
p11 = s1.add_parser('ms', parents=[pp])#, help='prep help')
p11.add_argument('MS_FILE', help = 'ms simulation file')
p12 = s1.add_parser('vcf', parents=[pp])#, help='prep help')
p12.add_argument('VCF_FILE', nargs='+') 
p12.add_argument('-s', '--samples', help='comma-separated list of sample names in VCF_FILE') 
p12.add_argument('--phased', action='store_true', help='alleles are phased in VCF_FILE') 
p1.set_defaults(func=prep)

p2 = s.add_parser('run', parents=[pp])#, help='run help')
p2.add_argument('DIR')
p2.add_argument('--geneflow', action='store_true', default = False, help = 'infer gene flow with msmc')
p2.add_argument('-n', type=int, default=0, help = 'number of segsep files to include')
p2.set_defaults(func=run)

args = p.parse_args()

loglevel = logging.WARNING
if args.verbose:
	loglevel = logging.INFO
if args.debug:
	loglevel = logging.DEBUG
logging.basicConfig(format = '%(module)s:%(lineno)d:%(levelname)s: %(message)s', level = loglevel)

args.func(args)
