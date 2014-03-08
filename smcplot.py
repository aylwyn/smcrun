#!/usr/bin/env python
# Aylwyn Scally 2014

import sys
import optparse
import os.path
import glob
import pandas as pd
import matplotlib.pyplot as plt

#colist = ['0xFFFFB300', '0xFF803E75', '0xFFFF6800', '0xFFA6BDD7', '0xFFC10020', '0xFFCEA262', '0xFF817066', '0xFF007D34', '0xFFF6768E', '0xFF00538A', '0xFFFF7A5C', '0xFF53377A', '0xFFFF8E00', '0xFFB32851', '0xFFF4C800', '0xFF7F180D', '0xFF93AA00', '0xFF593315', '0xFFF13A13', '0xFF232C16']
colist = [(0, 0, 1.0), (1.0, 0, 0), (0, 1.0, 0), (1.0, 1.0, 0), (1.0, 0, 1.0), (1.0, 0.5, 0.5), (0.5, 0.5, 0.5), (0.5, 0, 0), (1.0, 0.5, 0)]

p = optparse.OptionParser()
p.add_option('--hc', action='store_true', default = False)
p.add_option('--colours', default = '')
p.add_option('--recycle_colours', action='store_true', default = False)
p.add_option('--geneflow', action='store_true', default = False)
p.add_option('--coalrates', action='store_true', default = False)
p.add_option('--showvals', action='store_true', default = False)
p.add_option('--noplot', action='store_true', default = False)
p.add_option('-m', '--msmcfile', default = [], action='append')
p.add_option('-M', '--msmcdir', default = '')
p.add_option('-p', '--psmcfile', default = [], action='append')
p.add_option('-P', '--psmcdir', default = '')
p.add_option('-u', '--mugen', type='float', default = 1.25e-8)
p.add_option('-t', '--tgen', type='float', default = 25.0)
p.add_option('-s', '--simfile', default='')
p.add_option('-o', '--outname', default='psmc')

opt, args = p.parse_args()

#if len(args) > 0:
#	rfiles = args

if opt.msmcdir:
	opt.msmcfile += glob.glob(os.path.join(opt.msmcdir, '*.msmc.final.txt'))
if opt.psmcdir:
	opt.psmcfile += glob.glob(os.path.join(opt.msmcdir, '*.psmc.0.txt'))

if opt.colours:
	# read coldict file TODO
	coldict = {'gorberber': 'red', 'gorbergra': 'green', 'gorgorgor': 'blue'}
else:
	coldict = {}

#for rfile in glob.glob(os.path.join(rdir, dpref+'*.[1-7].msmc.final.txt')):
maxy = 0
maxx = 0
ic = 0
for ir, rfile in enumerate(opt.msmcfile):
	tb=pd.read_table(rfile, header = 0)
#	sname = os.path.splitext(os.path.basename(rfile))[0]
	sname = '.'.join(os.path.basename(rfile).split('.')[0:2])
	sname = rfile.split('.')[0]
	spname = sname.split('_')[0]
	print(sname)

	if coldict:
		icol = coldict[spname]
	else:
		icol = colist[ic]
	ic += 1

	rx = opt.tgen*(tb.ix[:,1])/opt.mugen
	rx[0] = max(rx[0], 1)
	li = len(rx) - 1
	if opt.coalrates:
		for ncol, lab in [(3, 'l00'), (4, 'l01'), (5, 'l11')]:#, (3, 'l00')] 
			ry = tb.ix[:,ncol]
			plt.step(rx, ry, label = lab) #TODO: fix initial point
#			plt.text(rx[li], ry[li], sname, fontsize=6)
	else:
		if opt.geneflow:
			ry = 2 * tb.ix[:,4] / (tb.ix[:,3] + tb.ix[:,5])
		else:
			ry = (1/tb.ix[:,3])/(2*opt.mugen)
		plt.step(rx, ry, label = sname, color=icol) #TODO: fix initial point
#	plt.step(opt.tgen*(tb.ix[:,1])/opt.mugen, (1/tb.ix[:,3])/(2*opt.mugen), label = sname, color=icol)
#	li = len(tb.ix[:,1]) - 1
	li = len(rx) - 1
#	plt.text(opt.tgen*(tb.ix[li,1])/opt.mugen, (1/tb.ix[li,3])/(2*opt.mugen), sname, color=icol, fontsize=6)
#	plt.text(rx[li], ry[li], sname, color=icol, fontsize=6)
	maxx = max(maxx, max(rx))
	maxy = max(maxy, max(ry))

if opt.recycle_colours:
	ic = 0

for ir, rfile in enumerate(opt.psmcfile):
	tb=pd.read_table(rfile, header = None)
	sname = os.path.splitext(os.path.basename(rfile))[0]
	sname = '.'.join(os.path.basename(rfile).split('.')[0:2])
	spname = sname.split('_')[0]
	print(sname)

	if coldict:
		icol = coldict[spname]
	else:
		icol = colist[ic]
	ic += 1

	rx = tb.ix[:,0]
	ry = 1e4*tb.ix[:,1]
	if opt.showvals:
		for ix, iy in zip(rx, ry):
			print('%f\t%f' % (ix, iy))
	plt.step(rx, ry, label = sname, color=icol)
	li = len(rx) - 1
	plt.text(rx[li], ry[li], sname, color=icol, fontsize=6)
	maxx = max(maxx, max(rx))
	maxy = max(maxy, max(ry))

if opt.simfile:
	tb=pd.read_table(opt.simfile, header = None)
	sname = 'sim'

	icol = 'black'
	rx = tb.ix[:,0]
	ry = tb.ix[:,1]
	plt.step(rx, ry, label = sname, color=icol, ls='--', where='post')
	li = len(rx) - 1
	plt.text(rx[li], ry[li], sname, color=icol, fontsize=6)
	maxx = max(maxx, max(rx))
	maxy = max(maxy, max(ry))


plt.xscale('log')
plt.xlim(1e3, maxx**1.1)
#plt.xlim(1e3, 1e7)
#plt.ylim(0, 40e3)
plt.xlabel('Time (y)')
if opt.geneflow:
	plt.ylim(0, 1.1)
	plt.ylabel(r'Relative coalescent rate')
else:
	plt.ylim(0, 1.1*maxy)
	plt.ylabel(r'$N_e$')

plt.legend(loc = 4, prop={'size':8})#, ncol = 2), fontsize = 'xx-small')

if not opt.noplot:
	if opt.hc:
		fname = opt.outname + '.pdf'
		sys.stderr.write('hardcopy in %s\n' % fname)
		plt.savefig(fname)
	else:
		plt.show()
