#!/usr/bin/env python
# Aylwyn Scally 2014

import sys
import optparse
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

p = optparse.OptionParser()
p.add_option('--hc', action='store_true', default = False)
p.add_option('--nodisplay', action='store_true', default = False)
#p.add_option('--colours', default = '')
p.add_option('--title', default = '')
#p.add_option('--recycle_colours', action='store_true', default = False)
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
p.add_option('--maxx', type='float', default = 0.0)
p.add_option('--legend', type='int', default = 0, help='0: no legend; else location specified as in matplotlib Legend()')
p.add_option('-s', '--simfile', default='')
p.add_option('-f', '--plotlist', default='', help = 'Each line:  ptype, path, name, colour, linestyle, linewidth')
p.add_option('--plotlist_bold', default='', help = 'plot file for bold lines (format as PLOTLIST)')
p.add_option('-o', '--outname', default='')
p.add_option('-v', '--verbose', action='store_true', default = False)

opt, args = p.parse_args()

if opt.nodisplay:
	matplotlib.use('Agg')
import matplotlib.pyplot as plt

loglevel = logging.WARNING
if opt.verbose:
	loglevel = logging.INFO
logging.basicConfig(format = '%(module)s:%(lineno)d:%(levelname)s: %(message)s', level = loglevel)

if opt.run_psmcplot or opt.msmcdir or opt.msmcfile:
	info('Scaling with mu = %e per gen and tgen = %.1f' % (opt.mugen, opt.tgen))

if opt.hc:
	from matplotlib import rcParams
	rcParams['axes.labelsize'] = 6
	rcParams['xtick.labelsize'] = 6
	rcParams['ytick.labelsize'] = 6
	rcParams['legend.fontsize'] = 6

inputlist = []
plotfiles = []

if opt.plotlist_bold:
	lalph = '0.3'
else:
	lalph = '1.0'

if opt.plotlist:
	for tok in (line.split() for line in open(opt.plotlist)):
		inputlist.append(SMCPlot(*tok + [lalph]))

if opt.plotlist_bold:
	for tok in (line.split() for line in open(opt.plotlist_bold)):
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
		plotfiles.append(splot)

fig = plt.figure()
if opt.hc:
	fsize = [2.25, 2.25]
	mar = [0.15, 0.15, 0.037, 0.036]
	margins = [mar[0], mar[1], 1-mar[2]-mar[0], 1-mar[3]-mar[1]]
	ax = plt.Axes(fig, margins)
	fig.add_axes(ax)

maxy = 0
maxx = 0
for ir, splot in enumerate(plotfiles):
	info(splot.pname)

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
			plt.step(rx, ry, label = splot.pname, color=splot.pcol, ls=splot.pls, lw=splot.plw, where='post', alpha=splot.palph)
	#	li = len(rx) - 1
	#	plt.text(opt.tgen*(tb.ix[li,1])/opt.mugen, (1/tb.ix[li,3])/(2*opt.mugen), sname, color=icol, fontsize=6)
	#	plt.text(rx[li], ry[li], sname, color=icol, fontsize=6)

	if splot.ptype == 'psmcfile':
		tb=pd.read_table(splot.path, header = None)

		rx = tb.ix[:,0]
		ry = 1e4*tb.ix[:,1]/1e3
		if opt.showvals:
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
plt.xlabel('Time (y)')
if opt.geneflow:
	plt.ylim(0, 1.1)
	plt.ylabel(r'Relative coalescent rate')
else:
	if opt.maxy <= 0.0:
		plt.ylim(0, 1.1*maxy)
	else:
		plt.ylim(0, 1.1*opt.maxy)
	if opt.maxx <= 0.0:
		plt.xlim(0, 1.1*maxx)
	else:
		plt.xlim(0, 1.1*opt.maxx)
	if not opt.coalrates:
		plt.ylabel(r'$N_e/1000$')

if opt.title:
	plt.title(opt.title)#, fontsize = 'small')

if opt.coalrates:
	plt.legend(('l00', 'l01', 'l11'), loc = 2, prop={'size':8})#, ncol = 2), fontsize = 'xx-small')
elif opt.legend:
	# remove duplicate labels
	handles, labels = plt.gca().get_legend_handles_labels()
	by_label = dict(zip(labels, handles))
	plt.legend(by_label.values(), by_label.keys(), loc=opt.legend)
#	plt.legend(loc = opt.legend)#, prop={'size':8})#, ncol = 2), fontsize = 'xx-small')

#plt.tight_layout()

if not opt.noplot:
	if opt.hc:
		fig.set_size_inches(fsize[0], fsize[1])
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
