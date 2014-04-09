#!/usr/bin/env python
# Aylwyn Scally 2014

import sys
import optparse
import os
import os.path
import glob
import pandas as pd
import matplotlib.pyplot as plt
import logging
from logging import error, warning, info, debug, critical

import aosutils

class SMCPlot(object):
	def __init__(self, ptype, path, pname, pcol='black'):#, pline)
		self.ptype = ptype
		self.path = path
		self.pname = pname
		self.pcol = pcol
#		self.pline = pline

	def show(self):
		print("\t".join([self.ptype, self.path, self.pname, self.pcol]))

p = optparse.OptionParser()
p.add_option('--hc', action='store_true', default = False)
#p.add_option('--colours', default = '')
p.add_option('--title', default = '')
p.add_option('--recycle_colours', action='store_true', default = False)
p.add_option('--geneflow', action='store_true', default = False)
p.add_option('--coalrates', action='store_true', default = False)
p.add_option('--nescale', action='store_true', default = False)
p.add_option('--showvals', action='store_true', default = False)
p.add_option('--noplot', action='store_true', default = False)
p.add_option('--psmc', action='store_true', default = False)
p.add_option('--msmc', action='store_true', default = False)
p.add_option('--run_psmcplot', action='store_true', default = False)
p.add_option('-m', '--msmcfile', default = [], action='append')
p.add_option('-M', '--msmcdir', default = [], action='append')
p.add_option('-p', '--psmcfile', default = [], action='append')
p.add_option('-P', '--psmcdir', default = [], action='append')
p.add_option('-u', '--mugen', type='float', default = 1.25e-8)
p.add_option('-t', '--tgen', type='float', default = 25.0)
p.add_option('--maxy', type='float', default = 0.0)
p.add_option('-s', '--simfile', default='')
p.add_option('-f', '--plotlist', default='')
p.add_option('-o', '--outname', default='')
p.add_option('-v', '--verbose', action='store_true', default = False)

opt, args = p.parse_args()

#if len(args) > 0:
#	rfiles = args
loglevel = logging.WARNING
if opt.verbose:
    loglevel = logging.INFO
#if opt.debug:
#    loglevel = logging.DEBUG
logging.basicConfig(format = '%(module)s:%(lineno)d:%(levelname)s: %(message)s', level = loglevel)

info('Scaling with mu = %e per gen and tgen = %.1f' % (opt.mugen, opt.tgen))

inputlist = []
plotfiles = []

if opt.plotlist:
	for tok in (line.split() for line in open(opt.plotlist)):
		inputlist.append(SMCPlot(*tok))

msmcdirs = opt.msmcdir
psmcdirs = opt.psmcdir
if opt.msmc:
	msmcdirs += args
if opt.psmc:
	psmcdirs += args
for x in msmcdirs:
	inputlist.append(SMCPlot('msmcdir', x, os.path.basename(x)))
for x in psmcdirs:
	inputlist.append(SMCPlot('psmcdir', x, os.path.basename(x)))

for x in opt.msmcfile:
	plotfiles.append(SMCPlot('msmcfile', x, x))
for x in opt.psmcfile:
	plotfiles.append(SMCPlot('psmcfile', x, x))
for x in opt.simfile:
	plotfiles.append(SMCPlot('simfile', x, x))

for splot in inputlist:
	resfiles = []
	if splot.ptype in ('msmcfile', 'psmcfile', 'simfile'):
		plotfiles.append(splot)
		continue

	if splot.ptype == 'msmcdir':
		resfiles = glob.glob(os.path.join(splot.path, '*.final.txt'))

	if splot.ptype == 'psmcdir':
		if opt.run_psmcplot:
			cwd = os.path.abspath('.')
			os.chdir(splot.path)
			for psmcfile in glob.glob('*psmc.out'):
				pref = psmcfile.replace('psmc.out', '') 
				if not pref:
					pref = '1'
				cmd = 'psmc_plot.pl -u %e -g %f -R %s %s' % (opt.mugen, opt.tgen, pref, psmcfile)
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
#TODO: add filename to splot.pname
#		splot.show()
		plotfiles.append(splot)

maxy = 0
maxx = 0
for ir, splot in enumerate(plotfiles):
#	print(splot.ptype)
	print(splot.pname)

#		if coldict:
#			icol = coldict[sname]
#		ic += 1
#		print(icol)

	if splot.ptype == 'msmcfile':
		tb=pd.read_table(splot.path, header = 0)
		rx = opt.tgen*(tb.ix[:,1])/opt.mugen
		rx[0] = max(rx[0], 1) # to account for log x axis
		li = len(rx) - 1
		if opt.coalrates:
			for ncol, lab in [(3, 'l00'), (4, 'l01'), (5, 'l11')]:#, (3, 'l00')] 
				if opt.nescale:
					ry = (1/tb.ix[:,ncol])/(2*opt.mugen)
				else:
					ry = tb.ix[:,ncol]
				plt.step(rx, ry, label = lab, color=coalrates_col[lab])
	#			plt.text(rx[li], ry[li], sname, fontsize=6)
				maxy = max(maxy, max(ry))
		else:
			if opt.geneflow:
				ry = 2 * tb.ix[:,4] / (tb.ix[:,3] + tb.ix[:,5])
			else:
				ry = (1/tb.ix[:,3])/(2*opt.mugen)
			plt.step(rx, ry, label = splot.pname, color=splot.pcol)
	#	li = len(rx) - 1
	#	plt.text(opt.tgen*(tb.ix[li,1])/opt.mugen, (1/tb.ix[li,3])/(2*opt.mugen), sname, color=icol, fontsize=6)
	#	plt.text(rx[li], ry[li], sname, color=icol, fontsize=6)

	if splot.ptype == 'psmcfile':
		tb=pd.read_table(splot.path, header = None)

		rx = tb.ix[:,0]
		ry = 1e4*tb.ix[:,1]
		if opt.showvals:
			for ix, iy in zip(rx, ry):
				print('%f\t%f' % (ix, iy))
		plt.step(rx, ry, label = splot.pname, color=splot.pcol)
#		li = len(rx) - 1
#		plt.text(rx[li], ry[li], splot.pname, color=splot.pcol, fontsize=6)

	if splot.ptype == 'simfile':
		tb=pd.read_table(splot.path, header = None)

		rx = tb.ix[:,0]
		ry = tb.ix[:,1]
		plt.step(rx, ry, label = splot.pname, color=splot.pcol, ls='--', where='post')
#		li = len(rx) - 1
#		plt.text(rx[li], ry[li], sname, color=icol, fontsize=6)

	maxx = max(maxx, max(rx))
	maxy = max(maxy, max(ry))

#	plt.title(sname)
#	plt.xscale('log')
#	plt.xlim(1e3, 2e6)
#	plt.ylim(0, 1.1)
#	plt.show()

plt.xscale('log')
plt.xlim(1e3, maxx**1.1)
#plt.xlim(1e3, 1e7)
#plt.ylim(0, 40e3)
plt.xlabel('Time (y)')
if opt.geneflow:
	plt.ylim(0, 1.1)
	plt.ylabel(r'Relative coalescent rate')
else:
	if opt.maxy <= 0.0:
		plt.ylim(0, 1.1*maxy)
	else:
		plt.ylim(0, 1.1*opt.maxy)
	if not opt.coalrates:
		plt.ylabel(r'$N_e$')

#if opt.coalrates:
#	plt.legend(('l00', 'l01', 'l11'), loc = 2, prop={'size':8})#, ncol = 2), fontsize = 'xx-small')
#else:
#	plt.legend(loc = 2, prop={'size':8})#, ncol = 2), fontsize = 'xx-small')

if opt.title:
	plt.title(opt.title, fontsize = 'small')

if not opt.noplot:
	if opt.hc:
		if opt.outname:
			fname = opt.outname + '.pdf'
		else:
			if opt.coalrates:
				fname = 'msmc-coalrates' + '.pdf'
			else:
				fname = 'smc' + '.pdf'

		sys.stderr.write('hardcopy in %s\n' % fname)
		plt.savefig(fname)
	else:
		plt.show()
