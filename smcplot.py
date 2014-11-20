#!/usr/bin/env python
# Aylwyn Scally 2014

import sys
import argparse
import os
import os.path
import glob
import pandas as pd
import matplotlib
import logging
from logging import error, warning, info, debug, critical

import aosutils

class SMCPlot(object):
	def __init__(self, ptype, path, pname, pcol='black', pls='-', plw='1', palph='1.0'):
		self.ptype = ptype
		self.path = path
		self.pname = pname
		if self.pname == '.':
			self.pname = os.path.basename(self.path).replace('.smcdir', '')
		self.pcol = pcol
		self.pls = pls
		self.plw = float(plw)/2.0
		self.palph = float(palph)

	def show(self):
		print("\t".join([self.ptype, self.path, self.pname, self.pcol]))

p = argparse.ArgumentParser()
p.add_argument('--hc', action='store_true', default = False)
p.add_argument('--nodisplay', action='store_true', default = False)
p.add_argument('--title', default = '')
#p.add_argument('--geneflow', action='store_true', default = False)
#p.add_argument('--coalrates', action='store_true', default = False)
#p.add_argument('--nescale', action='store_true', default = False)
p.add_argument('--showvals', action='store_true', default = False)
p.add_argument('--noplot', action='store_true', default = False)
#p.add_argument('--psmc', action='store_true', default = False)
#p.add_argument('--msmc', action='store_true', default = False)
p.add_argument('--run_psmcplot', action='store_true', default = False)
#p.add_argument('-m', '--msmcfile', default = [], action='append')
#p.add_argument('-M', '--msmcdir', default = [], action='append')
#p.add_argument('-p', '--psmcfile', default = [], action='append')
#p.add_argument('-P', '--psmcdir', default = [], action='append')
p.add_argument('-u', '--mugen', type=float, default = 1.25e-8)
p.add_argument('-t', '--tgen', type=float, default = 25.0)
p.add_argument('--maxy', type=float, default = 0.0)
p.add_argument('--maxx', type=float, default = 0.0)
p.add_argument('--legend', type=int, default = 0, help='<0: no legend; else location specified as in matplotlib Legend()')
p.add_argument('--label', default='')
p.add_argument('orient', nargs='?', choices=['vertical', 'horizontal'], default='vertical')
#p.add_argument('-s', '--simfile', default='')
p.add_argument('-f', '--plotlist', nargs = '*', default='', help = 'Each line:  ptype, path, name, colour, linestyle, linewidth')
p.add_argument('-o', '--outname', default='smc')
p.add_argument('-v', '--verbose', action='store_true', default = False)

args = p.parse_args()

if args.nodisplay:
	matplotlib.use('Agg')
import matplotlib.pyplot as plt

loglevel = logging.WARNING
if args.verbose:
	loglevel = logging.INFO
logging.basicConfig(format = '%(module)s:%(lineno)d:%(levelname)s: %(message)s', level = loglevel)

if args.run_psmcplot:# or args.msmcdir or args.msmcfile:
	info('Scaling with mu = %e per gen and tgen = %.1f' % (args.mugen, args.tgen))

nplots = len(args.plotlist)

if args.hc:
	from matplotlib import rcParams
	rcParams['axes.labelsize'] = 6
	rcParams['xtick.labelsize'] = 6
	rcParams['ytick.labelsize'] = 6
	rcParams['legend.fontsize'] = 6
	fsize = [2.25, 2.25*nplots]

if args.orient == 'vertical':
	fig, ax = plt.subplots(nplots, 1, squeeze = False)
elif args.orient ==  'horizontal':
	fig, ax = plt.subplots(1, nplots, squeeze = False)
axv = ax.reshape(-1)

sdir = os.path.abspath('.')
for ip, ppath in enumerate(args.plotlist):
	os.chdir(sdir)
	pdir, pfile = os.path.split(ppath)
	if pdir:
		info('switching to directory %s' % pdir)
		os.chdir(pdir)

	inputlist = []
	plotfiles = []
	for tok in (line.split() for line in open(pfile)):
		if tok[0].startswith('maxx'):
			args.maxx = float(tok[1])
		elif tok[0].startswith('maxy'):
			args.maxy = float(tok[1])
		elif tok[0].startswith('legend'):
			args.legend = int(tok[1])
		elif tok[0].startswith('label'):
			args.label = tok[1]
		else:
			inputlist.append(SMCPlot(*tok))

	#msmcdirs = args.msmcdir
	#psmcdirs = args.psmcdir
	#if args.msmc:
	#	msmcdirs += args
	#if args.psmc:
	#	psmcdirs += args
	#for x in msmcdirs:
	#	inputlist.append(SMCPlot('msmcdir', x, os.path.basename(x)))
	#for x in psmcdirs:
	#	inputlist.append(SMCPlot('psmcdir', x, os.path.basename(x)))
	#
	#for x in args.msmcfile:
	#	plotfiles.append(SMCPlot('msmcfile', x, x))
	#for x in args.psmcfile:
	#	plotfiles.append(SMCPlot('psmcfile', x, x))
	#for x in args.simfile:
	#	plotfiles.append(SMCPlot('simfile', x, x))

	for splot in inputlist:
		resfiles = []
		if splot.ptype in ('msmcfile', 'psmcfile', 'simfile'):
			plotfiles.append(splot)
			continue

		if splot.ptype == 'msmcdir':
			resfiles = glob.glob(os.path.join(splot.path, '*.final.txt'))

		if splot.ptype == 'psmcdir':
			if args.run_psmcplot:
				cwd = os.path.abspath('.')
				os.chdir(splot.path)
				for psmcfile in glob.glob('*psmc.out'):
					pref = psmcfile.replace('psmc.out', '') 
					if not pref:
						pref = '1'
					cmd = 'psmc_plot.pl -u %e -g %f -R %s %s' % (args.mugen, args.tgen, pref, psmcfile)
					info('running psmc_plot.pl in %s' % splot.path)
					aosutils.subcall(cmd, sim=False, wait = True)
				os.chdir(cwd)
			resfiles = glob.glob(os.path.join(splot.path, '*.0.txt'))

		if not resfiles:
			error('no results files found in %s' % splot.path)
		else:
			if len(resfiles) > 1:
				error('multiple results files found in %s; using first' % splot.path)
			splot.ptype = splot.ptype.replace('dir', 'file')
			splot.path = resfiles[0]
			plotfiles.append(splot)

	plt.sca(axv[ip])

	maxy = 0
	maxx = 0
	for ir, splot in enumerate(plotfiles):
		debug(splot.pname)

		if splot.ptype == 'msmcfile':
			tb=pd.read_table(splot.path, header = 0)
			rx = args.tgen*(tb.ix[:,1])/args.mugen
			rx[0] = max(rx[0], 1) # to account for log x axis
			li = len(rx) - 1
			if args.coalrates:
				for ncol, lab in [(3, 'l00'), (4, 'l01'), (5, 'l11')]:#, (3, 'l00')] 
					if args.nescale:
						ry = (1/tb.ix[:,ncol])/(2*args.mugen)
					else:
						ry = tb.ix[:,ncol]
					plt.step(rx, ry, label = lab, color=coalrates_col[lab])
		#			plt.text(rx[li], ry[li], sname, fontsize=6)
					maxy = max(maxy, max(ry))
			else:
				if args.geneflow:
					ry = 2 * tb.ix[:,4] / (tb.ix[:,3] + tb.ix[:,5])
				else:
					ry = (1/tb.ix[:,3])/(2*args.mugen)
				plt.step(rx, ry, label = splot.pname, color=splot.pcol, ls=splot.pls, lw=splot.plw, where='post', alpha=splot.palph)
		#	li = len(rx) - 1
		#	plt.text(args.tgen*(tb.ix[li,1])/args.mugen, (1/tb.ix[li,3])/(2*args.mugen), sname, color=icol, fontsize=6)
		#	plt.text(rx[li], ry[li], sname, color=icol, fontsize=6)

		if splot.ptype == 'psmcfile':
			tb=pd.read_table(splot.path, header = None)

			rx = tb.ix[:,0]
			ry = 1e4*tb.ix[:,1]/1e3
			if args.showvals:
				for ix, iy in zip(rx, ry):
					print('%f\t%f' % (ix, iy))
			plt.step(rx, ry, label = splot.pname, color=splot.pcol, ls=splot.pls, lw=splot.plw, where='post', alpha=splot.palph)

		if splot.ptype == 'simfile':
			tb=pd.read_table(splot.path, header = None)

			rx = tb.ix[:,0]
			ry = tb.ix[:,1]
			plt.step(rx, ry, label = splot.pname, color=splot.pcol, ls='--', where='post', alpha=splot.palph)

		maxx = max(maxx, max(rx))
		maxy = max(maxy, max(ry))

	plt.xscale('log')
	plt.xlim(1e3, maxx**1.1)
	plt.xlabel('Time (years)')
#	if args.geneflow:
#		plt.ylim(0, 1.1)
#		plt.ylabel(r'Relative coalescent rate')
#	else:
	if args.maxy <= 0.0:
		plt.ylim(0, 1.1*maxy)
	else:
		plt.ylim(0, 1.1*args.maxy)
	if args.maxx <= 0.0:
		plt.xlim(0, 1.1*maxx)
	else:
		plt.xlim(0, 1.1*args.maxx)
#	if not args.coalrates:
	plt.ylabel(r'$N_e/1000$')

	if args.title:
		plt.title(args.title)#, fontsize = 'small')

#	if args.coalrates:
#		plt.legend(('l00', 'l01', 'l11'), loc = 2, prop={'size':8})#, ncol = 2), fontsize = 'xx-small')
	if args.legend >= 0:
		# remove duplicate labels
		handles, labels = plt.gca().get_legend_handles_labels()
		by_label = dict(zip(labels, handles))
		lg = plt.legend(by_label.values(), by_label.keys(), loc=args.legend)
		lg.get_frame().set_linewidth(0.5) 
	#	plt.legend(loc = args.legend)#, prop={'size':8})#, ncol = 2), fontsize = 'xx-small')

	if args.label:
		plt.text(0.03, 0.97, args.label, transform=axv[ip].transAxes, fontsize=9, fontweight='bold', va='top')

os.chdir(sdir)

plt.subplots_adjust(left=0.35/fsize[0],bottom=0.375/fsize[1], right=1 - 0.075/fsize[0], top= 1 - 0.075/fsize[1], hspace = 0.25)

if not args.noplot:
	if args.hc:
		fig.set_size_inches(fsize[0], fsize[1])
		fname = args.outname + '.pdf'

		sys.stderr.write('hardcopy in %s\n' % fname)
		plt.savefig(fname)
	else:
		plt.show()
