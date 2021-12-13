#!/usr/bin/env python

"""This script converts between different tool output formats, the tools include, RFMIX, LOTER, LAMPLD, PCADMIX, SUPPORTMIX, CHROMOPAINTER, phased ELAI to RFMIX, phased MULTIMIX to LOTER, RFMIX, WINPOP and LAIT. unphased ELAI to WINPOP, and LAIT"""

import os
import sys
import fileinput
import gzip
import subprocess
import random
import numpy as np
import itertools
import csv
from functools import reduce

def convert_12(parameters):
	""" RFMIX (SNP/row & Hap/col) to LOTER (Hap/row & SNP/Cols)."""
	mapp = {}
	pop_index=[parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
	for i in pop_index:
		mapp[str(i+1)] = str(i)															#Convert interger dict keys to strings
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+"/rfmix/"+str(ch)+'/'
		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.rfmixloter.txt','w')		#Loter output format
		with open(cpath+parameters['outfile']+'.'+str(ch)+'.0.Viterbi.txt') as f:		#Read rfmix output format
			lis = [x.split() for x in f]
			for x in zip(*lis):										#transpose sequence of iterables:zip(*original_list)
				fw.write(' '.join([mapp[i] for i in x]))
				fw.write('\n')
		fw.close()


def convert_16(parameters):
	"""RFMIX to WINPOP"""
	nanc = len(parameters['anc_pop'])					# Number of ancestral populations
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/rfmix/'+str(ch)+'/'
		fp = open(cpath+parameters['outfile']+'.'+str(ch)+'.0.Viterbi.txt')
		nsnp = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath+ parameters['outfile']+'.'+str(ch)+".0.Viterbi.txt | sed '/^\s*#/d;/^\s*$/d' | wc -l ; exit 0"])
		nsnp = int(nsnp) # Number of SNPs

		rw = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath + parameters['outfile']+'.'+str(ch)+".0.Viterbi.txt | sed '/^\s*#/d;/^\s*$/d' | head -1 ; exit 0"]) 
		rw = str(rw)
		n = int(len(rw.strip().split())/2) # number of individuals

		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.rfmixwinpop.txt','w+')
		for i in range(n):
			for j in range(nanc):
				fw.write(("%d: "%(i,))+'\t'.join(nsnp*['0.0'])+'\n')

		fw.seek(0)
		llen = len(fw.readline())
		llenc = reduce(lambda x, y: x + y, [nanc*[llen + len(str(i))-1] for i in range(n)])
		fw.seek(0)
		fp.seek(0)

		snp = 0
		for line in fp:
			if not line.strip(): continue
			tline = [s.strip() for s in line.strip().split()]
			snp += 1
			for i in range(n):
				fw.seek(0)
				if tline[2*i]==tline[2*i+1]: 
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step)
					fw.write('1.0')
				else:
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
					fw.seek(0); j = int(tline[2*i+1])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
		fw.close(); fp.close()


def convert_210(parameters):
	"""LOTER to LAIT. ONLY works for len(anc_pop) <= 3. 2 stage conversion: LOTER to LAMPLD haps (remove space in LOTER format) & LAMPLD haps to LAIT by perl"""
	if len(parameters['anc_pop']) == 2 or len(parameters['anc_pop']) == 3:
		for ch in parameters['CHR']:
			cpath = parameters['infolder']+'/'+parameters['admix']+"/loter/"+str(ch)+'/'

			fo=open(cpath+parameters['outfile']+'.'+str(ch)+'.ancestry.txt')
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt','w')
			for line in fo:
				data=line.strip().split()
				fw.write(''.join(data)+'\n')
			fw.close()
			fo.close()
			##LAMPLD hap to LAIT (haplotypes in rows and SNPs in columns)
			os.system('perl ' +parameters['wrkdir']+'/franc_util/soft/lampld/standardizeOutput.pl lampld  '+str(len(parameters['anc_pop'])) +' '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt '+ cpath+parameters['outfile']+'.'+str(ch)+'.loterlait.txt')
			os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt')
	else:
		print("\nTo have LAIT output, requires a <= 3 way admixtures  ...\n")
		sys.exit(2)


def convert_110(parameters):
	"""RFMIX to LAIT. ONLY works for len(anc_pop) <= 3"""
	if len(parameters['anc_pop']) == 2 or len(parameters['anc_pop']) == 3:
		for ch in parameters['CHR']:
			cpath = parameters['infolder']+'/'+parameters['admix']+'/rfmix/'+str(ch)+'/'
			mapp = {}
			pop_index=[parameters['anc_pop'].index(i) for i in parameters['anc_pop']]	
			##RFMIX to LOTER
			for i in pop_index:
				mapp[str(i+1)] = str(i)													#Converts interger dict keys to strings
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.rfmixlo.txt','w')		#Loter output format
			with open(cpath+parameters['outfile']+'.'+str(ch)+'.0.Viterbi.txt') as f:	#Read rfmix output format
				lis = [x.split() for x in f]
				for x in zip(*lis):												#transpose seq of iterables:zip(*original_list)
					fw.write(' '.join([mapp[i] for i in x]))
					fw.write('\n')
			fw.close()
			##LOTER to LAIT
			fo = open(cpath+parameters['outfile']+'.'+str(ch)+'.rfmixlo.txt')
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt','w')
			for line in fo:
				data=line.strip().split()
				fw.write(''.join(data)+'\n')
			fw.close()
			fo.close()
			os.system('perl ' +parameters['wrkdir']+'/franc_util/soft/lampld/standardizeOutput.pl lampld  '+str(len(parameters['anc_pop'])) +' '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt '+ cpath+parameters['outfile']+'.'+str(ch)+'.rfmixlait.txt')
			os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.rfmixlo.txt')		#LAMPLD hap to LAIT (hap/row & SNP/col)
	else:
		print("\nTo have LAIT output, requires a <= 3 way admixtures  ...\n")
		sys.exit(2)


"""LOTER to RFMIX, WINPOP and LAIT (for 2 and 3 way admixtures)"""
def convert_21(parameters):
	"""LOTER (haps in rows, SNPs in columns) to RFMIX (SNPs in rows & haplotypes in columns)"""
	mapp = {}
	pop_index=[parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
	for i in pop_index:
		mapp[str(i)] = str(i + 1)															#Convert interger dict keys to strings
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+"/loter/"+str(ch)+'/'
		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.loterrfmix.txt','w')			#RFMIX file format
		with open(cpath+parameters['outfile']+'.'+str(ch)+'.ancestry.txt') as f:			#Read LOTER output
			lis = [x.split() for x in f] 													#split file contents into a list
			for x in zip(*lis):													#transpose seq of iterables:zip(*original_list)
				fw.write(' '.join([mapp[i] for i in x]))
				fw.write('\n')
		fw.close()


def convert_26(parameters):
	"""LOTER to WINPOP."""
	mapp = {}
	pop_index=[parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
	for i in pop_index:
		mapp[str(i)] = str(i + 1)														#Convert interger dict keys to strings
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+"/loter/"+str(ch)+'/'

		##LOTER to RFMIX
		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.loterrf.txt','w')			#RFMIX file format
		with open(cpath+parameters['outfile']+'.'+str(ch)+'.ancestry.txt') as f:		#Read LOTER output
			lis = [x.split() for x in f] 												#split file contents into a list
			for x in zip(*lis):														#transpose seq of iterables:zip(*original_list)
				fw.write(' '.join([mapp[i] for i in x]))
				fw.write('\n')
		fw.close()

		##RFMIX to WINPOP
		nanc = len(parameters['anc_pop'])					# Number of ancestral populations
		fp = open(cpath+parameters['outfile']+'.'+str(ch)+'.loterrf.txt')
		nsnp = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath+ parameters['outfile']+'.'+str(ch)+".loterrf.txt | sed '/^\s*#/d;/^\s*$/d' | wc -l ; exit 0"])
		nsnp = int(nsnp) # Number of SNPs

		rw = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath + parameters['outfile']+'.'+str(ch)+".loterrf.txt | sed '/^\s*#/d;/^\s*$/d' | head -1 ; exit 0"]) 
		rw = str(rw)
		n = int(len(rw.strip().split())/2) # number of individuals

		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.loterwinpop.txt','w+')
		for i in range(n):
			for j in range(nanc):
				fw.write(("%d: "%(i,))+'\t'.join(nsnp*['0.0'])+'\n')

		fw.seek(0)
		llen = len(fw.readline())
		llenc = reduce(lambda x, y: x + y, [nanc*[llen + len(str(i))-1] for i in range(n)])
		fw.seek(0)
		fp.seek(0)

		snp = 0
		for line in fp:
			if not line.strip(): continue
			tline = [s.strip() for s in line.strip().split()]
			snp += 1
			for i in range(n):
				fw.seek(0)
				if tline[2*i]==tline[2*i+1]: 
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step)
					fw.write('1.0')
				else:
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
					fw.seek(0); j = int(tline[2*i+1])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
		fw.close(); fp.close()
		os.system('rm -rf '+ cpath+parameters['outfile']+'.'+str(ch)+'.loterrf.txt')


def convert_31(parameters):
	"""PCADMIX to RFMIX. A 4 stage process: (1) exclude 1st column (2) Transpose to windows in rows & haps in columns, (3) Decode windows to SNPs in rows & (4) map 0-(K-1), PCADMIX to 1-K, RFMIX."""
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/pcadmix/'+str(ch)+'/'
		os.system("awk '{for (i=1; i<=NF-1; i++) $i = $(i+1); NF-=1; print}' " + cpath+parameters['outfile']+'.'+str(ch)+'.vit.txt > '+cpath+parameters['outfile']+'.'+str(ch)+'.vitRFMIX')						#Get rid of first column 
		##Transpose haplotype/row & windows/column to window/row and haplotype/column
		fo = cpath+parameters['outfile']+'.'+str(ch)+'.vitRFMIX'
		fw2 = open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIX','w')
		with open(fo) as f:
			lis = [x.split() for x in f]
			for x in zip(*lis):													#transpose seq of iterables:zip(*original_list)
				fw2.write(' '.join([i for i in x]))
				fw2.write('\n')
		fw2.close()
		#Decode from win to SNPs
		window = {}
		for line in fileinput.input(cpath+parameters['outfile']+'.'+str(ch)+'.markers.txt'):#loop thru analyzed markers for win len
			data = line.split()
			window[fileinput.lineno()] = data[1:]

		dat = open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIX').readlines()

		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIX.txt','w')
		for i in range(len(window)):
			for j in range(len(window[i+1])):
				fw.write(''.join(dat[i]))
		fw.close()

		##Map the 0 - K-1 to 1 - K
		mapp = {}
		pop_index = [parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
		for i in pop_index:
			mapp[str(i)] = str(i+1)												#Convert interger dict keys to strings
		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.pcadmixrfmix.txt','w')
		fo = open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIX.txt')
		for line in fo:
			tline=line.strip().split()
			fw.write(' '.join([mapp[i] for i in tline]))
			fw.write('\n')
		fw.close()
		fo.close()

		os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.vitRFMIX '+ cpath+parameters['outfile']+'.'+str(ch) + '.RFMIX '+ cpath+parameters['outfile']+'.'+str(ch)+'.RFMIX.txt')				#Delete unnecessary files


def convert_32(parameters):
	"""PCADMIX to LOTER. A 4 stages: (1) exclude first col (2) Transpose to window/row & hap/col, 
	(3) Convert window/row to SNP/row & (4) Transpose to  hap/row & SNP/column."""
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/pcadmix/'+str(ch)+'/'
		os.system("awk '{for (i=1; i<=NF-1; i++) $i = $(i+1); NF-=1; print}' " + cpath+parameters['outfile']+'.'+str(ch)+'.vit.txt > '+cpath+parameters['outfile']+'.'+str(ch)+'.vitRFMIX')						#Get rid of first col
		##Transpose from hap/row & win/col to win/row & hap/col
		fo = cpath+parameters['outfile']+'.'+str(ch)+'.vitRFMIX'
		fw2 = open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIXlike','w')
		with open(fo) as f:
			lis = [x.split() for x in f]
			for x in zip(*lis):													#transpose seq of iterables:zip(*original_list)
				fw2.write(' '.join([i for i in x]))
				fw2.write('\n')
		fw2.close()
		##Decode from windows to SNPs
		window = {}
		for line in fileinput.input(cpath+parameters['outfile']+'.'+str(ch)+'.markers.txt'):#loop thru analyzed markers 4 win len
			data = line.split()
			window[fileinput.lineno()] = data[1:]

		dat = open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIXlike').readlines()
		fw  =open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIXlike.txt','w')
		for i in range(len(window)):
			for j in range(len(window[i+1])):
				fw.write(''.join(dat[i]))
		fw.close()
		fo = cpath+parameters['outfile']+'.'+str(ch)+'.RFMIXlike.txt'
		fw = open( cpath+parameters['outfile']+'.'+str(ch)+'.pcadmixloter.txt','w')
		with open(fo) as f:
			lis = [x.split() for x in f]
			for x in zip(*lis):													#transpose seq of iterables:zip(*original_list)
				fw.write(' '.join([i for i in x]))
				fw.write('\n')
		fw.close()
		os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.vitRFMIX '+ cpath+parameters['outfile']+'.'+ str(ch)+'.RFMIXlike '+ cpath+parameters['outfile']+'.'+str(ch)+'.RFMIXlike.txt')		#Delete unnecessary files


def convert_36(parameters):
	"""PCADMIX to WINPOP"""
	for ch in parameters['CHR']:
		##PCADMIX to RFMIX
		cpath = parameters['infolder']+'/'+parameters['admix']+'/pcadmix/'+str(ch)+'/'
		os.system("awk '{for (i=1; i<=NF-1; i++) $i = $(i+1); NF-=1; print}' " + cpath+parameters['outfile']+'.'+str(ch)+'.vit.txt > '+cpath+parameters['outfile']+'.'+str(ch)+'.vitRFMIX')						#remove first col 
		##Transpose hap/row & wins/col to win/row and hap/col
		fo = cpath+parameters['outfile']+'.'+str(ch)+'.vitRFMIX'
		fw2 = open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIX','w')
		with open(fo) as f:
			lis = [x.split() for x in f]
			for x in zip(*lis):													#transpose seq of iterables:zip(*original_list)
				fw2.write(' '.join([i for i in x]))
				fw2.write('\n')
		fw2.close()
		##Decode from win to SNPs
		window = {}
		for line in fileinput.input(cpath+parameters['outfile']+'.'+str(ch)+'.markers.txt'):#loop thru analyzed markers 4 win len
			data = line.split()
			window[fileinput.lineno()] = data[1:]

		dat = open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIX').readlines()
		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIX.txt','w')
		for i in range(len(window)):
			for j in range(len(window[i+1])):
				fw.write(''.join(dat[i]))
		fw.close()
		##Map the 0 - K-1 to 1 - K
		mapp = {}
		pop_index = [parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
		for i in pop_index:
			mapp[str(i)] = str(i+1)												#Convert interger dict keys to strings
		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.pcadmixrfmix.txt','w')
		fo = open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIX.txt')
		for line in fo:
			tline=line.strip().split()
			fw.write(' '.join([mapp[i] for i in tline]))
			fw.write('\n')
		fw.close()
		fo.close()

		##RFMIX to WINPOP
		nanc = len(parameters['anc_pop'])					# Number of anc pop
		fp = open(cpath+parameters['outfile']+'.'+str(ch)+'.pcadmixrfmix.txt')
		nsnp = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath+ parameters['outfile']+'.'+str(ch)+".pcadmixrfmix.txt | sed '/^\s*#/d;/^\s*$/d' | wc -l ; exit 0"])
		nsnp = int(nsnp) # Number of SNPs

		rw = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath + parameters['outfile']+'.'+str(ch)+".pcadmixrfmix.txt | sed '/^\s*#/d;/^\s*$/d' | head -1 ; exit 0"])
		rw = str(rw)
		n = int(len(rw.strip().split())/2) # number of individuals

		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.pcadmixwinpop.txt','w+')
		for i in range(n):
			for j in range(nanc):
				fw.write(("%d: "%(i,))+'\t'.join(nsnp*['0.0'])+'\n')

		fw.seek(0)
		llen = len(fw.readline())
		llenc = reduce(lambda x, y: x + y, [nanc*[llen + len(str(i))-1] for i in range(n)])
		fw.seek(0)
		fp.seek(0)

		snp = 0
		for line in fp:
			if not line.strip(): continue
			tline = [s.strip() for s in line.strip().split()]
			snp += 1
			for i in range(n):
				fw.seek(0)
				if tline[2*i]==tline[2*i+1]: 
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step)
					fw.write('1.0')
				else:
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
					fw.seek(0); j = int(tline[2*i+1])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
		fw.close(); fp.close()
		os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.vitRFMIX '+ cpath+parameters['outfile']+'.'+str(ch) + '.RFMIX '+ cpath+parameters['outfile']+'.'+str(ch)+'.RFMIX.txt '+ cpath+parameters['outfile']+'.'+str(ch) + '.pcadmixrfmix.txt')#Delete unnecessary files

def convert_310(parameters):
	"""PCADMIX to LAIT. ONLY works if the len(anc_pop) <= 3. (1) PCADMIX to LOTER, (2) LOTER to LAIT """
	if len(parameters['anc_pop']) == 2 or len(parameters['anc_pop']) == 3:
		for ch in parameters['CHR']:
			cpath = parameters['infolder']+'/'+parameters['admix']+'/pcadmix/'+str(ch)+'/'

			##PCADMIX to LOTER
			os.system("awk '{for (i=1; i<=NF-1; i++) $i = $(i+1); NF-=1; print}' " + cpath+parameters['outfile']+'.'+str(ch)+'.vit.txt > '+cpath+parameters['outfile']+'.'+str(ch)+'.vitRFMIX')						#Get rid of first col
			##Transpose from hap/row & win/col to win/row & hap/col
			fo = cpath+parameters['outfile']+'.'+str(ch)+'.vitRFMIX'
			fw2 = open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIXlike','w')
			with open(fo) as f:
				lis = [x.split() for x in f]
				for x in zip(*lis):												#transpose seq of iterables:zip(*original_list)
					fw2.write(' '.join([i for i in x]))
					fw2.write('\n')
			fw2.close()
			#Decode from windows to SNPs
			window = {}
			for line in fileinput.input(cpath+parameters['outfile']+'.'+str(ch)+'.markers.txt'):#determin win len, analyzed markers
				data = line.split()
				window[fileinput.lineno()] = data[1:]

			dat = open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIXlike').readlines()
			fw  =open(cpath+parameters['outfile']+'.'+str(ch)+'.RFMIXlike.txt','w')
			for i in range(len(window)):
				for j in range(len(window[i+1])):
					fw.write(''.join(dat[i]))
			fw.close()
			fo = cpath+parameters['outfile']+'.'+str(ch)+'.RFMIXlike.txt'
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.pcadmixlo.txt','w')
			with open(fo) as f:
				lis = [x.split() for x in f]
				for x in zip(*lis):												#transpose seq of iterables:zip(*original_list)
					fw.write(' '.join([i for i in x]))
					fw.write('\n')
			fw.close()
			os.system('rm -r '+cpath+parameters['outfile']+'.'+str(ch)+'.vitRFMIX '+ cpath+parameters['outfile']+'.'+ str(ch)+'.RFMIXlike '+ cpath+parameters['outfile']+'.'+str(ch)+'.RFMIXlike.txt')		#Delete unnecessary files

			##convert LOTER to LAIT
			fo=open(cpath+parameters['outfile']+'.'+str(ch)+'.pcadmixlo.txt')
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt','w')
			for line in fo:
				data=line.strip().split()
				fw.write(''.join(data)+'\n')
			fw.close()
			fo.close()
			##LAMPLD hap to LAIT (haplotypes in rows and SNPs in columns)
			os.system('perl ' +parameters['wrkdir']+'/franc_util/soft/lampld/standardizeOutput.pl lampld '+str(len(parameters['anc_pop'])) +' '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt '+ cpath+parameters['outfile']+'.'+str(ch)+'.pcadmixlait.txt')
			os.system('rm -r '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.pcadmixlo.txt')

	else:
		print("\nLAIT standardizes ancestry in <= 3 way admixtures  ...\n")
		sys.exit(2)


def convert_41(parameters):
	"""SUPPORTMIX to RFMIX. SUPPORTMIX. LA is in .tped file, has win/row & hap/col, 1st 4 cols are head. 4 stages conversion. 
	(1) number of SNPs in each window, (2) Exclude 1st 4 cols, (3) anc/win to anc/SNP & (4) map 0-(K-1) to 1-K"""
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/supportmix/'+str(ch)+'/'
		os.system('cp '+parameters['infolder']+'/'+parameters['admix']+'.'+str(ch)+'.bim '+cpath+'output/')
		##Determining number of SNPs in each window
		SNP = []; MARKERS = {}
		for line in fileinput.input(cpath+'output/'+parameters['admix']+'.'+str(ch)+'.bim'):#read SNPs file with all SNPs in order
			data = line.split()
			SNP.append(data[1])												#keep track of SNP ID's in SNPs file be4 processing
		for line in fileinput.input(cpath+'output/'+parameters['outfile']+'.'+str(ch)+'.positions.txt'):#loop thru SVM position file
			data = line.split()
			if fileinput.lineno () > 1:										#exclude header
				if data[0] in SNP:											#Check if the start rsID is in the SNP list
					start = SNP.index(data[0])								#assigning start to the index of data[0] 
				if data[1] in SNP:											#Check if data[1] i.e, "rsId" in pos is in SNPs file
					end = SNP.index(data[1])								#assign rsID end to index of data[1]
				MARKERS[fileinput.lineno()-1] = SNP[start:end+1]			#considering line #s as keys assign values as start:end

		##Remove 1st 4 cols of tped & save in temporary file 
		os.system("awk '{for (i=1; i<=NF-4; i++) $i = $(i+4); NF-=4; print}' " + cpath+'output/'+parameters['outfile']+'.'+str(ch)+'.tped' + ' > '+cpath+parameters['admix']+'.'+str(ch)+'.tpedlo')
		##Decoding from windows to SNPs/rows i.e writting each window the number of SNP size in that window
		fw = open(cpath+parameters['admix']+'.'+str(ch)+'.RFMIXlike', 'w')

		data = open(cpath+parameters['admix']+'.'+str(ch)+'.tpedlo').readlines()
		for i in range(len(MARKERS)):
			for j in range(len(MARKERS[i+1])):
				fw.write(''.join(data[i]))
		fw.close()
		##Map 0-K-1:SUPPORTMIX to 1-K:RFMIX
		mapp = {}
		pop_index = [parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
		for i in pop_index:
			mapp[str(i)] = str(i+1)										#Convert interger dictionary keys to strings
		fw = open(cpath+ parameters['outfile']+'.'+str(ch)+'.supportmixrfmix.txt','w')
		fo = open(cpath+parameters['admix']+'.'+str(ch)+'.RFMIXlike')
		for line in fo:
			tline=line.strip().split()
			fw.write(' '.join([mapp[i] for i in tline]))				#fw.write(' ' + ' '.join(map(lambda x: mapp[x], g)))
			fw.write('\n')
		fw.close()
		fo.close()
		os.system('rm -rf '+cpath+parameters['admix']+'.'+str(ch)+'.RFMIXlike '+ cpath+parameters['admix']+'.'+str(ch)+'.tpedlo')


def convert_42(parameters):
	"""SUPPORTMIX to LOTER. A 4 stage process (1) find # of SNPs/win, (2) Remove 1st 4 cols, (3) Convert win to SNP/row & 
	(4) transpose to hap/row & SNP/column."""
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/supportmix/'+str(ch)+'/'
		SNP = []; MARKERS = {}											#SNPs for each win & all as lists & markers as dictionary
		os.system('cp '+parameters['infolder']+parameters['admix']+'.'+str(ch)+'.bim '+cpath+'output/')
		##Determining number of SNPs /window
		for line in fileinput.input(cpath+'output/'+parameters['admix']+'.'+str(ch)+'.bim'):	#read SNPs file under analysis
			data = line.split()
			SNP.append(data[0])											#keep track of SNP ID's in SNPs file be4 processing
		for line in fileinput.input(cpath+'output/'+parameters['outfile']+'.'+str(ch)+'.positions.txt'):#loop thru SVM pos file
			data = line.split()
			if fileinput.lineno () > 1:									#exclude header
				if data[0] in SNP:										#Check if the start rsID is in the SNP list
					start = SNP.index(data[0])							#assigning start to the index of data[0] 
				if data[1] in SNP:										#Check if data[1]; end "rsId" in .pos is in SNPs file
					end = SNP.index(data[1])							#assign end to the index of data[1]
				MARKERS[fileinput.lineno()-1] = SNP[start:end+1]		#consider line #s as keys assign values as start:end

		##Remove 1st 4 cols of supportmix tped file & save in temporary file
		os.system("awk '{for (i=1; i<=NF-4; i++) $i = $(i+4); NF-=4; print}' " + cpath+'output/'+parameters['outfile']+'.'+str(ch)+'.tped' + ' > '+cpath+parameters['admix']+'.'+str(ch)+'.tpedlo')
		##Decode win/row to SNP/row i.e write the # SNPs in each window
		fw = open(cpath+parameters['admix']+'.'+str(ch)+'.RFMIXlike', 'w')

		data = open(cpath+parameters['admix']+'.'+str(ch)+'.tpedlo').readlines()
		for i in range(len(MARKERS)):
			for j in range(len(MARKERS[i+1])):
				fw.write(''.join(data[i]))
		fw.close()
		##Transposing
		fw2=open(cpath+ parameters['outfile']+'.'+str(ch)+'.supportmixloter.txt', 'w')
		with open(cpath+parameters['admix']+'.'+str(ch)+'.RFMIXlike') as f:
			lis = [x.split() for x in f]
			for x in zip(*lis):
				fw2.write(' '.join([i for i in x]))
				fw2.write('\n')
		fw2.close()
		os.system('rm -rf '+cpath+parameters['admix']+'.'+str(ch)+'.RFMIXlike '+cpath+parameters['admix']+'.'+str(ch)+'.tpedlo '+cpath+'output/'+parameters['outfile']+'.'+str(ch)+'.positions.txt '+cpath+'output/'+parameters['admix']+'.'+str(ch)+'.bim')


def convert_46(parameters):
	"""SUPPORTMIX to WINPOP. """
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/supportmix/'+str(ch)+'/'
		os.system('cp '+parameters['infolder']+parameters['admix']+'.'+str(ch)+'.bim '+cpath+'output/')
		#SUPPORTMIX to RFMIX
		SNP = []; MARKERS = {}
		for line in fileinput.input(cpath+'output/'+parameters['admix']+'.'+str(ch)+'.bim'):#read SNPs file with all SNPs
			data = line.split()
			SNP.append(data[1])												#keep track of SNP ID's in SNPs file be4 processing
		for line in fileinput.input(cpath+'output/'+parameters['outfile']+'.'+str(ch)+'.positions.txt'):#loop thru SVM position file
			data = line.split()
			if fileinput.lineno () > 1:										#exclude header
				if data[0] in SNP:											#Check if the start rsID is in the SNP list
					start = SNP.index(data[0])								#assign start to the index of data[0] 
				if data[1] in SNP:											#Check if "rsId" in pos is in SNPs file
					end = SNP.index(data[1])								#assign rsID end to index of data[1]
				MARKERS[fileinput.lineno()-1] = SNP[start:end+1]			#line #s as keys store start:end as values

		os.system("awk '{for (i=1; i<=NF-4; i++) $i = $(i+4); NF-=4; print}' " + cpath+'output/'+parameters['outfile']+'.'+str(ch)+'.tped' + ' > '+cpath+parameters['admix']+'.'+str(ch)+'.tpedlo')		##Remove 1st 4 cols of tped & save in temporary file
		##Decoding from wins to SNPs/rows i.e writting each window the number of SNP size in that window
		fw = open(cpath+parameters['admix']+'.'+str(ch)+'.RFMIXlike', 'w')

		data = open(cpath+parameters['admix']+'.'+str(ch)+'.tpedlo').readlines()
		for i in range(len(MARKERS)):
			for j in range(len(MARKERS[i+1])):
				fw.write(''.join(data[i]))
		fw.close()
		##Map 0-K-1:SUPPORTMIX to 1-K:RFMIX
		mapp = {}
		pop_index = [parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
		for i in pop_index:
			mapp[str(i)] = str(i+1)										#Convert interger dictionary keys to strings
		fw = open(cpath+ parameters['outfile']+'.'+str(ch)+'.supportmixrf.txt','w')
		fo = open(cpath+parameters['admix']+'.'+str(ch)+'.RFMIXlike')
		for line in fo:
			tline=line.strip().split()
			fw.write(' '.join([mapp[i] for i in tline]))				#fw.write(' ' + ' '.join(map(lambda x: mapp[x], g)))
			fw.write('\n')
		fw.close()
		fo.close()

		##RFMIX to WINPOP
		nanc = len(parameters['anc_pop'])					# Number of ancestral populations
		fp = open(cpath+parameters['outfile']+'.'+str(ch)+'.supportmixrf.txt')
		nsnp = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath+ parameters['outfile']+'.'+str(ch)+".supportmixrf.txt | sed '/^\s*#/d;/^\s*$/d' | wc -l ; exit 0"])
		nsnp = int(nsnp) # Number of SNPs

		rw = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath + parameters['outfile']+'.'+str(ch)+".supportmixrf.txt | sed '/^\s*#/d;/^\s*$/d' | head -1 ; exit 0"]) #fp.readline()
		rw = str(rw)
		n = int(len(rw.strip().split())/2) # number of individuals

		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.supportmixwinpop.txt','w+')
		for i in range(n):
			for j in range(nanc):
				fw.write(("%d: "%(i,))+'\t'.join(nsnp*['0.0'])+'\n')

		fw.seek(0)
		llen = len(fw.readline())
		llenc = reduce(lambda x, y: x + y, [nanc*[llen + len(str(i))-1] for i in range(n)])
		fw.seek(0)
		fp.seek(0)

		snp = 0
		for line in fp:
			if not line.strip(): continue
			tline = [s.strip() for s in line.strip().split()]
			snp += 1
			for i in range(n):
				fw.seek(0)
				if tline[2*i]==tline[2*i+1]: 
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step)
					fw.write('1.0')
				else:
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
					fw.seek(0); j = int(tline[2*i+1])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
		fw.close(); fp.close()
		os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.Vit.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.supportmixrf.txt '+cpath+'output/'+parameters['admix']+'.'+str(ch)+'.bim '+cpath+'output/'+parameters['outfile']+'.'+str(ch)+'.positions.txt '+cpath+parameters['admix']+'.'+str(ch)+'.RFMIXlike '+ cpath+parameters['admix']+'.'+str(ch)+'.tpedlo')


def convert_410(parameters):
	"""SUPPORTMIX to LAIT. ONLY works if the len(+parameters['anc_pop']) <= 3. (1) SUPPORTMIX to LOTER, (2) LOTER to LAIT """
	if len(parameters['anc_pop']) == 2 or len(parameters['anc_pop']) == 3:
		for ch in parameters['CHR']:
			#SUPPORTMIX to RFMIX
			cpath = parameters['infolder']+'/'+parameters['admix']+'/supportmix/'+str(ch)+'/'
			os.system('cp '+parameters['infolder']+parameters['admix']+'.'+str(ch)+'.bim '+cpath+'output/')
			##Determining number of SNPs in each window
			SNP = []; MARKERS = {}
			for line in fileinput.input(cpath+'output/'+parameters['admix']+'.'+str(ch)+'.bim'):#read SNPs file
				data = line.split()
				SNP.append(data[1])												#track of SNP ID's in SNPs file be4 processing
			for line in fileinput.input(cpath+'output/'+parameters['outfile']+'.'+str(ch)+'.positions.txt'):#position file
				data = line.split()
				if fileinput.lineno () > 1:										#exclude header
					if data[0] in SNP:											#Check if the start rsID is in the SNP list
						start = SNP.index(data[0])								#assigning start to the index of data[0] 
					if data[1] in SNP:											#Check if data[1] i.e, "rsId" in pos in SNPs file
						end = SNP.index(data[1])								#assign rsID end to index of data[1]
					MARKERS[fileinput.lineno()-1] = SNP[start:end+1]			#line #s as keys assign values start:end

			##Remove 1st 4 cols of tped & save in temporary file 
			os.system("awk '{for (i=1; i<=NF-4; i++) $i = $(i+4); NF-=4; print}' " + cpath+'output/'+parameters['outfile']+'.'+str(ch)+'.tped' + ' > '+cpath+parameters['admix']+'.'+str(ch)+'.tpedlo')
			##Decoding from win to SNPs/rows i.e writting number of SNPs in each window
			fw = open(cpath+parameters['admix']+'.'+str(ch)+'.RFMIXlike', 'w')

			data = open(cpath+parameters['admix']+'.'+str(ch)+'.tpedlo').readlines()

			for i in range(len(MARKERS)):
				for j in range(len(MARKERS[i+1])):
					fw.write(''.join(data[i]))
			fw.close()
			##Map 0-K-1:SUPPORTMIX to 1-K:RFMIX
			mapp = {}
			pop_index = [parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
			for i in pop_index:
				mapp[str(i)] = str(i+1)										#Convert interger dictionary keys to strings
			fw = open(cpath+ parameters['outfile']+'.'+str(ch)+'.supportmixrf.txt','w')
			fo = open(cpath+parameters['admix']+'.'+str(ch)+'.RFMIXlike')
			for line in fo:
				tline=line.strip().split()
				fw.write(' '.join([mapp[i] for i in tline]))
				fw.write('\n')
			fw.close()
			fo.close()
			os.system('rm -rf '+cpath+parameters['admix']+'.'+str(ch)+'.RFMIXlike '+ cpath+parameters['admix']+'.'+str(ch)+'.tpedlo')

			##RFMIX to LAIT
			mapp = {}
			pop_index=[parameters['anc_pop'].index(i) for i in parameters['anc_pop']]	
			##RFMIX to LOTER
			for i in pop_index:
				mapp[str(i+1)] = str(i)														#Convert interger dict keys to strings
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.rfmixlo.txt','w')			#Loter output format
			with open(cpath+parameters['outfile']+'.'+str(ch)+'.supportmixrf.txt') as f:	#Read rfmix output format
				lis = [x.split() for x in f]
				for x in zip(*lis):
					fw.write(' '.join([mapp[i] for i in x]))
					fw.write('\n')
			fw.close()
			##LOTER to lampldhap
			fo = open(cpath+parameters['outfile']+'.'+str(ch)+'.rfmixlo.txt')
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt','w')
			for line in fo:
				data=line.strip().split()
				fw.write(''.join(data)+'\n')
			fw.close()
			fo.close()
			##lampldhap to LAIT
			os.system('perl ' +parameters['wrkdir']+'/franc_util/soft/lampld/standardizeOutput.pl lampld  '+str(len(parameters['anc_pop'])) +' '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt '+ cpath+parameters['outfile']+'.'+str(ch)+'.supportmixlait.txt')
			os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.rfmixlo.txt '+ cpath+parameters['outfile']+'.'+str(ch)+'.supportmixlo.txt')

	else:
		print("\nLAIT standardizes ancestry in <= 3 way admixtures  ...\n")
		sys.exit(2)


def convert_610(parameters):
	"""WINPOP to LAIT"""
	if len(parameters['anc_pop']) == 2 or len(parameters['anc_pop']) == 3:
		for ch in parameters['CHR']:
			cpath = parameters['infolder']+'/'+parameters['admix']+'/winpop/'+str(ch)+'/'
			os.system('perl ' +parameters['wrkdir']+'/franc_util/soft/lampld/standardizeOutput.pl lamp '+str(len(parameters['anc_pop'])) +' '+ cpath+'ancestry_lamp4.out ' + cpath+parameters['outfile']+'.'+str(ch)+'.winpoplait.txt')
	else:
		print("\nLAIT standardizes ancestry in <= 3 way admixtures  ...\n")
		sys.exit(2)

def convert_72(parameters):
	"""LAMPLD to LOTER"""
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/lampld/'+str(ch)+'/'
		os.system(parameters['wrkdir']+'/franc_util/soft/lampld/convertLAMPLDout.pl ' + cpath +parameters['outfile']+'.'+str(ch)+'.out '+ cpath +parameters['outfile']+'.'+str(ch)+'.hap')
		fo = open(cpath+parameters['outfile']+'.'+str(ch)+'.hap')
		fw=open(cpath +parameters['outfile']+'.'+str(ch)+'.lampldloter.txt','w')
		for line in fo:
			fw.write(' '.join(line))
		fw.close()
		fo.close()
		os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.hap')


def convert_71(parameters):
	"""LAMPLD to RFMIX"""
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/lampld/'+str(ch)+'/'
		##LAMPLD to LOTER
		os.system(parameters['wrkdir']+'/franc_util/soft/lampld/convertLAMPLDout.pl ' + cpath +parameters['outfile']+'.'+str(ch)+'.out '+ cpath +parameters['outfile']+'.'+str(ch)+'.hap')
		fo = open(cpath+parameters['outfile']+'.'+str(ch)+'.hap')
		fw=open(cpath +parameters['outfile']+'.'+str(ch)+'.lampldlo.txt','w')
		for line in fo:
			fw.write(' '.join(line))
		fw.close()
		fo.close()

		##LOTER to RFMIX
		mapp = {}
		pop_index=[parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
		for i in pop_index:
			mapp[str(i)] = str(i + 1)
		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.lampldrfmix.txt','w')			#RFMIX file format
		with open(cpath+parameters['outfile']+'.'+str(ch)+'.lampldlo.txt') as f:			#Read LOTER output
			lis = [x.split() for x in f]
			for x in zip(*lis):													#transpose seq of iterables:zip(*original_list)
				fw.write(' '.join([mapp[i] for i in x]))
				fw.write('\n')
		fw.close()
		os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldlo.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.hap')


def convert_76(parameters):
	"""LAMPLD to WINPOP"""
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/lampld/'+str(ch)+'/'
		##LAMPLD to RFMIX
		##1. LAMPLD to LOTER
		os.system(parameters['wrkdir']+'/franc_util/soft/lampld/convertLAMPLDout.pl ' + cpath +parameters['outfile']+'.'+str(ch)+'.out '+ cpath +parameters['outfile']+'.'+str(ch)+'.hap')
		fo = open(cpath+parameters['outfile']+'.'+str(ch)+'.hap')
		fw=open(cpath +parameters['outfile']+'.'+str(ch)+'.lampldlo.txt','w')
		for line in fo:
			fw.write(' '.join(line))
		fw.close()
		fo.close()
		##2. LOTER to RFMIX
		mapp = {}
		pop_index=[parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
		for i in pop_index:
			mapp[str(i)] = str(i + 1)
		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.lampldrf.txt','w')			#RFMIX file format
		with open(cpath+parameters['outfile']+'.'+str(ch)+'.lampldlo.txt') as f:		#Read LOTER output
			lis = [x.split() for x in f] 												#split file contents into a list
			for x in zip(*lis):													#transpose seq of iterables:zip(*original_list)
				fw.write(' '.join([mapp[i] for i in x]))
				fw.write('\n')
		fw.close()

		##RFMIX to WINPOP
		nanc = len(parameters['anc_pop'])					# Number of ancestral populations
		fp = open(cpath+parameters['outfile']+'.'+str(ch)+'.lampldrf.txt')
		nsnp = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath+ parameters['outfile']+'.'+str(ch)+".lampldrf.txt | sed '/^\s*#/d;/^\s*$/d' | wc -l ; exit 0"])
		nsnp = int(nsnp) # Number of SNPs

		rw = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath + parameters['outfile']+'.'+str(ch)+".lampldrf.txt | sed '/^\s*#/d;/^\s*$/d' | head -1 ; exit 0"]) #fp.readline()
		rw = str(rw)
		n = int(len(rw.strip().split())/2) # number of individuals

		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.lampldwinpop.txt','w+')
		for i in range(n):
			for j in range(nanc):
				fw.write(("%d: "%(i,))+'\t'.join(nsnp*['0.0'])+'\n')

		fw.seek(0)
		llen = len(fw.readline())
		llenc = reduce(lambda x, y: x + y, [nanc*[llen + len(str(i))-1] for i in range(n)])
		fw.seek(0)
		fp.seek(0)

		snp = 0
		for line in fp:
			if not line.strip(): continue
			tline = [s.strip() for s in line.strip().split()]
			snp += 1
			for i in range(n):
				fw.seek(0)
				if tline[2*i]==tline[2*i+1]: 
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step)
					fw.write('1.0')
				else:
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
					fw.seek(0); j = int(tline[2*i+1])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
		fw.close(); fp.close()
		os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldlo.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.hap '+cpath+parameters['outfile']+'.'+str(ch)+'.Vit.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldrf.txt')


def convert_710(parameters):
	"""Converts LAMPLD to LOTER. (1) Convert from LAMPLD to LAMPLD haps then (2) LAMPLD haps to LAIT"""
	if len(parameters['anc_pop']) == 2 or len(parameters['anc_pop']) == 3:
		for ch in parameters['CHR']:
			cpath = parameters['infolder']+'/'+parameters['admix']+'/lampld/'+str(ch)+'/'
			##LAMPLD haps to LAIT
			os.system(parameters['wrkdir']+'/franc_util/soft/lampld/convertLAMPLDout.pl ' + cpath +parameters['outfile']+'.'+str(ch)+'.out '+ cpath +parameters['outfile']+'.'+str(ch)+'.hap')
			os.system('perl ' +parameters['wrkdir']+'/franc_util/soft/lampld/standardizeOutput.pl lampld  '+str(len(parameters['anc_pop'])) +' '+ cpath +parameters['outfile']+'.'+str(ch)+'.hap '+ cpath+parameters['outfile']+'.'+str(ch)+'.lampldlait.txt')
			#Delete unnecessary files
			os.system('rm -rf '+ cpath+parameters['outfile']+'.'+str(ch)+'.hap')
	else:
		print("\nLAIT standardizes ancestry in <= 3 way admixtures  ...\n")
		sys.exit(2)


def convert_82(parameters):
	"""ELAI phased to LOTER"""
	nanc = len(parameters['anc_pop'])				# Number of ancestral populations
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/elai/'+str(ch)+'/'
		if parameters['phased'] == 'YES':
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.elairloter.txt','w')
			fp = open(cpath+parameters['outfile']+'.'+str(ch)+'.ps21.txt')
			for line in fp:
				tstr = ''
				tmp = line.strip().split()
				for i in range(int(len(tmp)/nanc)):		#iterating thru number of SNPs
					tancprop = [float(s) for s in tmp[(nanc*i):(nanc*(i+1))]]			#slicing SNPs to have individual SNPs e.g for 2 pop we have [0:K] as SNP1, entries [K:4] SNP2 etc
					maxanc = max(tancprop)					# capture the haplotype with a maximum
					ancindex = [s[0] for s in enumerate(tancprop) if s[1] == maxanc] #enumerate gives the sliced values indices, 
					tstr += str(*random.sample(ancindex, 1)) #why are we considering a random sample
				fw.write(' '.join(tstr)+'\n')
			fw.close(); fp.close()
		else:
		##Check if its a 2 or 3 way admixture to convert to LAIT
			if len(parameters['anc_pop']) == 2 or len(parameters['anc_pop']) == 3:
				##ELAI unphased to LAIT
				os.system('perl ' +parameters['wrkdir']+'/franc_util/soft/lampld/standardizeOutput.pl elai unphased  '+cpath+parameters['outfile']+'.'+str(ch)+'.ps21.txt '+ cpath+parameters['outfile']+'.'+str(ch)+'.elailait.txt')


def convert_81(parameters):
	"""ELAI phased to RFMIX"""
	nanc = len(parameters['anc_pop'])				# Number of ancestral populations
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/elai/'+str(ch)+'/'
		if parameters['phased'] == 'YES':
			##ELAI to LOTER
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.elailoter.txt','w')
			fp = open(cpath+parameters['outfile']+'.'+str(ch)+'.ps21.txt')
			for line in fp:
				tstr = ''
				tmp = line.strip().split()
				for i in range(int(len(tmp)/nanc)):		#iterating thru number of SNPs
					tancprop = [float(s) for s in tmp[(nanc*i):(nanc*(i+1))]]			#slicing SNPs to have individual SNPs e.g for 2 pop we have [0:K] as SNP1, entries [K:4] SNP2 etc
					maxanc = max(tancprop)					# capture the haplotype with a maximum
					ancindex = [s[0] for s in enumerate(tancprop) if s[1] == maxanc] #enumerate gives the sliced values indices, 
					tstr += str(*random.sample(ancindex, 1)) #why are we considering a random sample
				fw.write(' '.join(tstr)+'\n')
			fw.close(); fp.close()

			##LOTER to RFMIX
			mapp = {}
			pop_index=[parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
			for i in pop_index:
				mapp[str(i)] = str(i + 1)
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.elairfmix.txt','w')			#RFMIX file format
			with open(cpath+parameters['outfile']+'.'+str(ch)+'.elailoter.txt') as f:		#Read LOTER output
				lis = [x.split() for x in f] 												#split file contents into a list
				for x in zip(*lis):													#transpose seq of iterables:zip(*original_list)
					fw.write(' '.join([mapp[i] for i in x]))
					fw.write('\n')
			fw.close()
			os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.elailoter.txt')
		else:
			##Check if its a 2 or 3 way admixture to convert to LAIT
			if len(parameters['anc_pop']) == 2 or len(parameters['anc_pop']) == 3:
				##ELAI unphased to LAIT
				os.system('perl ' +parameters['wrkdir']+'/franc_util/soft/lampld/standardizeOutput.pl elai unphased  '+cpath+parameters['outfile']+'.'+str(ch)+'.ps21.txt '+ cpath+parameters['outfile']+'.'+str(ch)+'.elailait.txt')


def convert_86(parameters):
	"""ELAI phased to RFMIX"""
	nanc = len(parameters['anc_pop'])				# Number of ancestral populations
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/elai/'+str(ch)+'/'
		if parameters['phased'] == 'YES':
			##ELAI to LOTER
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.elailoter.txt','w')
			fp = open(cpath+parameters['outfile']+'.'+str(ch)+'.ps21.txt')
			for line in fp:
				tstr = ''
				tmp = line.strip().split()
				for i in range(int(len(tmp)/nanc)):		#iterating thru number of SNPs
					tancprop = [float(s) for s in tmp[(nanc*i):(nanc*(i+1))]]			#slicing SNPs to have indiv SNPs e.g for 2 pop we have [0:2] as SNP1, entries [2:4] SNP2 etc
					maxanc = max(tancprop)					# capture the haplotype with a maximum
					ancindex = [s[0] for s in enumerate(tancprop) if s[1] == maxanc] #enumerate gives the sliced values indices, 
					tstr += str(*random.sample(ancindex, 1)) #why are we considering a random sample
				fw.write(' '.join(tstr)+'\n')
			fw.close(); fp.close()

			##LOTER to RFMIX
			mapp = {}
			pop_index=[parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
			for i in pop_index:
				mapp[str(i)] = str(i + 1)
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.elairfmix.txt','w')			#RFMIX file format
			with open(cpath+parameters['outfile']+'.'+str(ch)+'.elailoter.txt') as f:		#Read LOTER output
				lis = [x.split() for x in f] 												#split file contents into a list
				for x in zip(*lis):													#transpose seq of iterables:zip(*original_list)
					fw.write(' '.join([mapp[i] for i in x]))
					fw.write('\n')
			fw.close()

			##RFMIX to WINPOP
			fp = open(cpath+parameters['outfile']+'.'+str(ch)+'.elairfmix.txt')
			nsnp = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath+ parameters['outfile']+'.'+str(ch)+".elairfmix.txt | sed '/^\s*#/d;/^\s*$/d' | wc -l ; exit 0"])
			nsnp = int(nsnp) # Number of SNPs

			rw = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath + parameters['outfile']+'.'+str(ch)+".elairfmix.txt | sed '/^\s*#/d;/^\s*$/d' | head -1 ; exit 0"]) #fp.readline()
			rw = str(rw)
			n = int(len(rw.strip().split())/2) # number of individuals

			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.elaiwinpop.txt','w+')
			for i in range(n):
				for j in range(nanc):
					fw.write(("%d: "%(i,))+'\t'.join(nsnp*['0.0'])+'\n')

			fw.seek(0)
			llen = len(fw.readline())
			llenc = reduce(lambda x, y: x + y, [nanc*[llen + len(str(i))-1] for i in range(n)])
			fw.seek(0)
			fp.seek(0)

			snp = 0
			for line in fp:
				if not line.strip(): continue
				tline = [s.strip() for s in line.strip().split()]
				snp += 1
				for i in range(n):
					fw.seek(0)
					if tline[2*i]==tline[2*i+1]: 
						j = int(tline[2*i])
						step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
						fw.seek(step)
						fw.write('1.0')
					else:
						j = int(tline[2*i])
						step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
						fw.seek(step); fw.write('0.5')
						fw.seek(0); j = int(tline[2*i+1])
						step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
						fw.seek(step); fw.write('0.5')
			fw.close(); fp.close()
			os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.elailoter.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.elairfmix.txt')

		else:
			##Check if its a 2 or 3 way admixture to convert to LAIT
			if len(parameters['anc_pop']) == 2 or len(parameters['anc_pop']) == 3:
				##ELAI unphased to LAIT
				os.system('perl ' +parameters['wrkdir']+'/franc_util/soft/lampld/standardizeOutput.pl elai unphased  '+cpath+parameters['outfile']+'.'+str(ch)+'.ps21.txt '+ cpath+parameters['outfile']+'.'+str(ch)+'.elailait.txt')


def convert_810(parameters):
	"""Converts ELAI to LAIT. a is the ELAI output format and b is the LAIT format"""
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/elai/'+str(ch)+'/'
		if len(parameters['anc_pop']) == 2 or len(parameters['anc_pop']) == 3:
			if parameters['phased'] == 'YES':
				##ELAI (unphased) to LAIT 
				os.system('perl ' +parameters['wrkdir']+'/franc_util/soft/lampld/standardizeOutput.pl elai phased  '+cpath+parameters['outfile']+'.'+str(ch)+'.ps21.txt '+ cpath+parameters['outfile']+'.'+str(ch)+'.elailait.txt')
			else:
				os.system('perl ' +parameters['wrkdir']+'/franc_util/soft/lampld/standardizeOutput.pl elai unphased  '+cpath+parameters['outfile']+'.'+str(ch)+'.ps21.txt '+ cpath+parameters['outfile']+'.'+str(ch)+'.elailait.txt')
		else:
			print("\nLAIT standardizes ancestry in <= 3 way admixtures  ...\n")
			sys.exit(2)


import re
def sorted_nicely(l):
	"""This function sort the list in ascending order"""
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	return sorted(l, key = alphanum_key)


def convert_91(parameters):
	"""MULTIMIX phased and resolved to LOTER. Arguments is the multimixrfmix file to write"""
	import glob
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/multimix/'+str(ch)+'/'
		##MULTIMIX to RFMIX
		file_list = glob.glob(cpath+'out_'+parameters['model']+'/resolved_calls_indiv*')
		z=sorted_nicely(file_list)														#sort to ind1_haplo1,ind1_haplo2,...
		cols = [4]																		# extracting column 4 add more columns here
		data = []
		for f in z:
			data.append(np.loadtxt(f, usecols=cols))
		arr = np.vstack(data)
		np.savetxt(cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt',arr,fmt="%i") 	#saving array in a file
		fo = cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt'
		fw2 = open(cpath+parameters['outfile']+'.'+str(ch) +'.multimixrfmix.txt','w')
		with open(fo) as f:
			lis = [x.split() for x in f]
			for x in zip(*lis):												#transpose sequence of iterables:zip(*original_list)
				fw2.write(' '.join([i for i in x]))
				fw2.write('\n')
		fw2.close()
		os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt')


def convert_92(parameters):
	"""MULTIMIX phased and resolved to LOTER"""
	import glob
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/multimix/'+str(ch)+'/'
		#MULTIMIX to RFMIX		
		file_list = glob.glob(cpath+'out_'+parameters['model']+'/resolved_calls_indiv*')
		z=sorted_nicely(file_list)														#sort to ind1_haplo1,ind1_haplo2,...
		cols = [4]																		# extracting column 4 add more columns here
		data = []
		for f in z:
			data.append(np.loadtxt(f, usecols=cols))
		arr = np.vstack(data)
		np.savetxt(cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt',arr,fmt="%i") 	#saving array in a file
		fo = cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt'
		fw2 = open(cpath+parameters['outfile']+'.'+str(ch) +'.multimixrfm.txt','w')
		with open(fo) as f:
			lis = [x.split() for x in f]
			for x in zip(*lis):
				fw2.write(' '.join([i for i in x]))
				fw2.write('\n')
		fw2.close()
		os.system('rm -r '+cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt')
		##RFMIX to LOTER	
		mapp = {}
		pop_index=[parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
		for i in pop_index:
			mapp[str(i+1)] = str(i)											#Convert interger dict keys to strings
		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.elailoter.txt','w')					#Loter output format
		with open(cpath+parameters['outfile']+'.'+str(ch)+'.multimixrfm.txt') as f:				#Read rfmix output format
			lis = [x.split() for x in f]
			for x in zip(*lis):
				fw.write(' '.join([mapp[i] for i in x]))
				fw.write('\n')
		fw.close()
		os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt '+ cpath+parameters['outfile']+'.'+str(ch) +'.multimixrfm.txt')


def convert_96(parameters):
	"MULTIMIX phased and resolved to WINPOP"
	import glob
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/multimix/'+str(ch)+'/'
		#MULTIMIX to RFMIX
		file_list = glob.glob(cpath+'out_'+parameters['model']+'/resolved_calls_indiv*')
		z=sorted_nicely(file_list)														#sort to ind1_haplo1,ind1_haplo2,...
		cols = [4]																		# extracting column 4 add more columns here
		data = []
		for f in z:
			data.append(np.loadtxt(f, usecols=cols))
		arr = np.vstack(data)
		np.savetxt(cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt',arr,fmt="%i") 	#saving array in a file
		fo = cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt'
		fw2 = open(cpath+parameters['outfile']+'.'+str(ch) +'.multimixrfmix.txt','w')
		with open(fo) as f:
			lis = [x.split() for x in f]
			for x in zip(*lis):
				fw2.write(' '.join([i for i in x]))
				fw2.write('\n')
		fw2.close()

		##RFMIX to WINPOP
		nanc = len(parameters['anc_pop'])					# Number of ancestral populations
		fp = open(cpath+parameters['outfile']+'.'+str(ch)+'.multimixrfmix.txt')
		nsnp = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath+ parameters['outfile']+'.'+str(ch)+".multimixrfmix.txt | sed '/^\s*#/d;/^\s*$/d' | wc -l ; exit 0"])
		nsnp = int(nsnp) # Number of SNPs

		rw = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath + parameters['outfile']+'.'+str(ch)+".multimixrfmix.txt | sed '/^\s*#/d;/^\s*$/d' | head -1 ; exit 0"]) #fp.readline()
		rw = str(rw)
		n = int(len(rw.strip().split())/2) # number of individuals

		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.multimixwinpop.txt','w+')
		for i in range(n):
			for j in range(nanc):
				fw.write(("%d: "%(i,))+'\t'.join(nsnp*['0.0'])+'\n')

		fw.seek(0)
		llen = len(fw.readline())
		llenc = reduce(lambda x, y: x + y, [nanc*[llen + len(str(i))-1] for i in range(n)])
		fw.seek(0)
		fp.seek(0)

		snp = 0
		for line in fp:
			if not line.strip(): continue
			tline = [s.strip() for s in line.strip().split()]
			snp += 1
			for i in range(n):
				fw.seek(0)
				if tline[2*i]==tline[2*i+1]: 
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step)
					fw.write('1.0')
				else:
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
					fw.seek(0); j = int(tline[2*i+1])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
		fw.close(); fp.close()
		os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt'+cpath+parameters['outfile']+'.'+str(ch)+'.multimixrfmix.txt')


def convert_910(parameters):
	""" MULTIMIX phased and resolved to LAIT"""
	import glob
	if len(parameters['anc_pop']) == 2 or len(parameters['anc_pop']) == 3:
		for ch in parameters['CHR']:
			cpath = parameters['infolder']+'/'+parameters['admix']+'/multimix/'+str(ch)+'/'
			#MULTIMIX to RFMIX
			file_list = glob.glob(cpath+'out_'+parameters['model']+'/resolved_calls_indiv*')
			z=sorted_nicely(file_list)							#sort to ind1_haplo1,ind1_haplo2,...
			cols = [4]											# extracting column 4 add more columns here
			data = []
			for f in z:
				data.append(np.loadtxt(f, usecols=cols))
			arr = np.vstack(data)
			np.savetxt(cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt',arr,fmt="%i") 	#saving array in a file
			fo = cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt'
			fw2 = open(cpath+parameters['outfile']+'.'+str(ch) +'.multimixrfmix.txt','w')
			with open(fo) as f:
				lis = [x.split() for x in f]
				for x in zip(*lis):								#transpose seq of iterables:zip(*original_list)
					fw2.write(' '.join([i for i in x]))
					fw2.write('\n')
			fw2.close()
			os.system('rm -r '+cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt')
			##RFMIX to LAIT, 1.RFMIX to LOTER & 2.LOTER to LAIT
			mapp = {}
			pop_index=[parameters['anc_pop'].index(i) for i in parameters['anc_pop']]	
			for i in pop_index:
				mapp[str(i+1)] = str(i)													#Converts interger dict keys to strings
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.rfmixlo.txt','w')		#Loter output format
			with open(cpath+parameters['outfile']+'.'+str(ch)+'.multimixrfmix.txt') as f:	#Read rfmix output format
				lis = [x.split() for x in f]
				for x in zip(*lis):												#transpose seq of iterables:zip(*original_list)
					fw.write(' '.join([mapp[i] for i in x]))
					fw.write('\n')
			fw.close()
			fo = open(cpath+parameters['outfile']+'.'+str(ch)+'.rfmixlo.txt')
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt','w')
			for line in fo:
				data=line.strip().split()
				fw.write(''.join(data)+'\n')
			fw.close()
			fo.close()
			os.system('perl ' +parameters['wrkdir']+'/franc_util/soft/lampld/standardizeOutput.pl lampld  '+str(len(parameters['anc_pop'])) +' '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt '+ cpath+parameters['outfile']+'.'+str(ch)+'.multimixlait.txt')

			os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.rfmixlo.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.multimixrf.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.multimixrfmix.txt')
	else:
		print("\nLAIT standardizes ancestry in <= 3 way admixtures  ...\n")
		sys.exit(2)


def convert_52(parameters):
	"""CHROMOPAINTER to LOTER"""
	##CHROMOPAINTER to LOTER
	nanc = len(parameters['anc_pop'])
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/chromopainter/'+str(ch)+'/'
		try:
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt','w')
			with gzip.open(cpath+parameters['outfile']+'.'+str(ch)+'.copyprobsperlocus.out.gz', 'rb') as f:
				header = f.readline().strip().split()[1:]
				tstr = ''
				for line in f:
					tline = line.strip()
					if tline.startswith('HAP'): # This means that we have a new hap
						if tstr: 
							fw.write(' '.join(tstr)+'\n')
							tstr = ''
						continue
					tline = line.strip().split()
					tline = [float(s.strip()) for s in tline[1:(nanc+1)]]
					maxanc = max(tline)
					ancindex = [s[0] for s in enumerate(tline) if s[1] == maxanc]
					tstr += str(*random.sample(ancindex, 1))
				fw.write(' '.join(tstr)+'\n')
			fw.close()
			nlines = subprocess.check_output(["/bin/sh", "-c", "wc -l %s; exit 0"%(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt',)])
			nlines = int(nlines.split()[0].strip())
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt','w')
			for i in range(nlines//2,0,-1):
				hap1 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i-1,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
				hap2 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
				fw.write('%s\n%s\n'%(hap1.strip(), hap2.strip()))
			fw.close()
		except:				
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt','w')
			with gzip.open(cpath+parameters['outfile']+'.'+str(ch)+'.copyprobsperlocus.out.gz', 'rb') as f:
				header = f.readline().strip().split()[1:]
				tstr = ''
				for line in f:
					tline = line.strip()
					if tline.startswith(b'HAP'): # This means that we have a new hap
						if tstr: 
							fw.write(' '.join(tstr)+'\n')
							tstr = ''
						continue
					tline = line.strip().split()
					tline = [float(s.strip()) for s in tline[1:(nanc+1)]]
					maxanc = max(tline)
					ancindex = [s[0] for s in enumerate(tline) if s[1] == maxanc]
					tstr += str(*random.sample(ancindex, 1))
				fw.write(' '.join(tstr)+'\n')
			fw.close()

			nlines = subprocess.check_output(["/bin/sh", "-c", "wc -l %s; exit 0"%(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt',)])
			nlines = int(nlines.split()[0].strip())
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterbytes.txt','w')
			for i in range(nlines//2,0,-1):
				hap1 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i-1,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
				hap2 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
				fw.write('%s\n%s\n'%(hap1.strip(), hap2.strip()))
			fw.close()
			fo=open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterbytes.txt')
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt','w')
			for line in fo:
				fw.write(line[2:-4] +'\n')
			fw.close();fo.close()
			os.system('rm -r '+cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterbytes.txt')
		os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt ')


def convert_51(parameters):
	"""CHROMOPAINTER TO RFMIX"""
	##CHROMOPAINTER to LOTER
	nanc = len(parameters['anc_pop'])
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/chromopainter/'+str(ch)+'/'
		try:
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt','w')
			with gzip.open(cpath+parameters['outfile']+'.'+str(ch)+'.copyprobsperlocus.out.gz', 'rb') as f:
				header = f.readline().strip().split()[1:]
				tstr = ''
				for line in f:
					tline = line.strip()
					if tline.startswith('HAP'): # This means that we have a new hap
						if tstr: 
							fw.write(' '.join(tstr)+'\n')
							tstr = ''
						continue
					tline = line.strip().split()
					tline = [float(s.strip()) for s in tline[1:(nanc+1)]]
					maxanc = max(tline)
					ancindex = [s[0] for s in enumerate(tline) if s[1] == maxanc]
					tstr += str(*random.sample(ancindex, 1))
				fw.write(' '.join(tstr)+'\n')
			fw.close()
			nlines = subprocess.check_output(["/bin/sh", "-c", "wc -l %s; exit 0"%(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt',)])
			nlines = int(nlines.split()[0].strip())
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt','w')
			for i in range(nlines//2,0,-1):
				hap1 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i-1,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
				hap2 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
				fw.write('%s\n%s\n'%(hap1.strip(), hap2.strip()))
			fw.close()
		except:				
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt','w')
			with gzip.open(cpath+parameters['outfile']+'.'+str(ch)+'.copyprobsperlocus.out.gz', 'rb') as f:
				header = f.readline().strip().split()[1:]
				tstr = ''
				for line in f:
					tline = line.strip()
					if tline.startswith(b'HAP'): # This means that we have a new hap
						if tstr: 
							fw.write(' '.join(tstr)+'\n')
							tstr = ''
						continue
					tline = line.strip().split()
					tline = [float(s.strip()) for s in tline[1:(nanc+1)]]
					maxanc = max(tline)
					ancindex = [s[0] for s in enumerate(tline) if s[1] == maxanc]
					tstr += str(*random.sample(ancindex, 1))
				fw.write(' '.join(tstr)+'\n')
			fw.close()

			nlines = subprocess.check_output(["/bin/sh", "-c", "wc -l %s; exit 0"%(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt',)])
			nlines = int(nlines.split()[0].strip())
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterbytes.txt','w')
			for i in range(nlines//2,0,-1):
				hap1 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i-1,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
				hap2 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
				fw.write('%s\n%s\n'%(hap1.strip(), hap2.strip()))
			fw.close()
			fo=open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterbytes.txt')
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt','w')
			for line in fo:
				fw.write(line[2:-4] +'\n')
			fw.close();fo.close()
			os.system('rm -r '+cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterbytes.txt')
		##LOTER to RFMIX 
		mapp = {}
		pop_index=[parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
		for i in pop_index:
			##Converting the interger dictionary keys to strings
			mapp[str(i)] = str(i + 1)
		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterrfmix.txt','w')				#RFMIX file format to write
		with open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt') as f:		#Reading LOTER output
			lis = [x.split() for x in f] 										#spliting contents of a file into a list
			for x in zip(*lis):													#transpose sequence of iterables:zip(*original_list)
				fw.write(' '.join([mapp[i] for i in x]))
				fw.write('\n')
		fw.close()
		os.system('rm -r '+cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt ')


def convert_56(parameters):
	"""CHROMOPAINTER to WINPOP"""
	##CHROMOPAINTER to LOTER
	nanc = len(parameters['anc_pop'])
	for ch in parameters['CHR']:
		cpath = parameters['infolder']+'/'+parameters['admix']+'/chromopainter/'+str(ch)+'/'
		try:
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt','w')
			with gzip.open(cpath+parameters['outfile']+'.'+str(ch)+'.copyprobsperlocus.out.gz', 'rb') as f:
				header = f.readline().strip().split()[1:]
				tstr = ''
				for line in f:
					tline = line.strip()
					if tline.startswith('HAP'): # This means that we have a new hap
						if tstr: 
							fw.write(' '.join(tstr)+'\n')
							tstr = ''
						continue
					tline = line.strip().split()
					tline = [float(s.strip()) for s in tline[1:(nanc+1)]]
					maxanc = max(tline)
					ancindex = [s[0] for s in enumerate(tline) if s[1] == maxanc]
					tstr += str(*random.sample(ancindex, 1))
				fw.write(' '.join(tstr)+'\n')
			fw.close()
			nlines = subprocess.check_output(["/bin/sh", "-c", "wc -l %s; exit 0"%(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt',)])
			nlines = int(nlines.split()[0].strip())
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt','w')
			for i in range(nlines//2,0,-1):
				hap1 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i-1,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
				hap2 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
				fw.write('%s\n%s\n'%(hap1.strip(), hap2.strip()))
			fw.close()
		except:				
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt','w')
			with gzip.open(cpath+parameters['outfile']+'.'+str(ch)+'.copyprobsperlocus.out.gz', 'rb') as f:
				header = f.readline().strip().split()[1:]
				tstr = ''
				for line in f:
					tline = line.strip()
					if tline.startswith(b'HAP'): # This means that we have a new hap
						if tstr: 
							fw.write(' '.join(tstr)+'\n')
							tstr = ''
						continue
					tline = line.strip().split()
					tline = [float(s.strip()) for s in tline[1:(nanc+1)]]
					maxanc = max(tline)
					ancindex = [s[0] for s in enumerate(tline) if s[1] == maxanc]
					tstr += str(*random.sample(ancindex, 1))
				fw.write(' '.join(tstr)+'\n')
			fw.close()

			nlines = subprocess.check_output(["/bin/sh", "-c", "wc -l %s; exit 0"%(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt',)])
			nlines = int(nlines.split()[0].strip())
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterbytes.txt','w')
			for i in range(nlines//2,0,-1):
				hap1 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i-1,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
				hap2 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
				fw.write('%s\n%s\n'%(hap1.strip(), hap2.strip()))
			fw.close()
			fo=open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterbytes.txt')
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt','w')
			for line in fo:
				fw.write(line[2:-4] +'\n')
			fw.close();fo.close()
			os.system('rm -r '+cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterbytes.txt')

		##LOTER to RFMIX 
		mapp = {}
		pop_index=[parameters['anc_pop'].index(i) for i in parameters['anc_pop']]
		for i in pop_index:
			##Converting the interger dictionary keys to strings
			mapp[str(i)] = str(i + 1)
		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterrfmix.txt','w')				#RFMIX file format to write
		with open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt') as f:		#Reading LOTER output
			lis = [x.split() for x in f] 										#spliting contents of a file into a list
			for x in zip(*lis):													#transpose sequence of iterables:zip(*original_list)
				fw.write(' '.join([mapp[i] for i in x]))
				fw.write('\n')
		fw.close()
		os.system('rm -r '+cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt ')	##Delete uneccessary files

		#RFMIX to WINPOP
		nanc = len(parameters['anc_pop'])					# Number of ancestral populations
		fp = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterrfmix.txt')
		nsnp = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath+ parameters['outfile']+'.'+str(ch)+".chromopainterrfmix.txt | sed '/^\s*#/d;/^\s*$/d' | wc -l ; exit 0"])
		nsnp = int(nsnp) # Number of SNPs

		rw = subprocess.check_output(["/bin/sh", "-c", "cat "+cpath + parameters['outfile']+'.'+str(ch)+".chromopainterrfmix.txt | sed '/^\s*#/d;/^\s*$/d' | head -1 ; exit 0"]) #fp.readline()
		rw = str(rw)
		n = int(len(rw.strip().split())/2) # number of individuals

		fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterwinpop.txt','w+')
		for i in range(n):
			for j in range(nanc):
				fw.write(("%d: "%(i,))+' '.join(nsnp*['0.0'])+'\n')

		fw.seek(0)
		llen = len(fw.readline())
		llenc = reduce(lambda x, y: x + y, [nanc*[llen + len(str(i))-1] for i in range(n)])
		fw.seek(0)
		fp.seek(0)

		snp = 0
		for line in fp:
			if not line.strip(): continue
			tline = [s.strip() for s in line.strip().split()]
			snp += 1
			for i in range(n):
				fw.seek(0)
				if tline[2*i]==tline[2*i+1]: 
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step)
					fw.write('1.0')
				else:
					j = int(tline[2*i])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
					fw.seek(0); j = int(tline[2*i+1])
					step = sum(llenc[:nanc*i]) + sum(llenc[(nanc*i):(nanc*i+j-1)]) + 4*snp-2+len(str(i))
					fw.seek(step); fw.write('0.5')
		fw.close(); fp.close()
		os.system('rm -r '+cpath+ parameters['outfile']+'.'+str(ch) +'.chromopainterrfmix.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')


"""This function only works if the len(anc_pop) <= 3"""
def convert_510(parameters):
	""" CHROMOPAINTER to LAIT"""
	if len(parameters['anc_pop']) == 2 or len(parameters['anc_pop']) == 3:
		for ch in parameters['CHR']:
			cpath = parameters['infolder']+'/'+parameters['admix']+'/chromopainter/'+str(ch)+'/'
			## CHROMOPAINTER to LOTER
			nanc = len(parameters['anc_pop'])
			try:
				fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt','w')
				with gzip.open(cpath+parameters['outfile']+'.'+str(ch)+'.copyprobsperlocus.out.gz', 'rb') as f:
					header = f.readline().strip().split()[1:]
					tstr = ''
					for line in f:
						tline = line.strip()
						if tline.startswith('HAP'): # This means that we have a new hap
							if tstr: 
								fw.write(' '.join(tstr)+'\n')
								tstr = ''
							continue
						tline = line.strip().split()
						tline = [float(s.strip()) for s in tline[1:(nanc+1)]]
						maxanc = max(tline)
						ancindex = [s[0] for s in enumerate(tline) if s[1] == maxanc]
						tstr += str(*random.sample(ancindex, 1))
					fw.write(' '.join(tstr)+'\n')
				fw.close()
				nlines = subprocess.check_output(["/bin/sh", "-c", "wc -l %s; exit 0"%(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt',)])
				nlines = int(nlines.split()[0].strip())
				fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt','w')
				for i in range(nlines//2,0,-1):
					hap1 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i-1,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
					hap2 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
					fw.write('%s\n%s\n'%(hap1.strip(), hap2.strip()))
				fw.close()
			except:				
				fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt','w')
				with gzip.open(cpath+parameters['outfile']+'.'+str(ch)+'.copyprobsperlocus.out.gz', 'rb') as f:
					header = f.readline().strip().split()[1:]
					tstr = ''
					for line in f:
						tline = line.strip()
						if tline.startswith(b'HAP'): # This means that we have a new hap
							if tstr: 
								fw.write(' '.join(tstr)+'\n')
								tstr = ''
							continue
						tline = line.strip().split()
						tline = [float(s.strip()) for s in tline[1:(nanc+1)]]
						maxanc = max(tline)
						ancindex = [s[0] for s in enumerate(tline) if s[1] == maxanc]
						tstr += str(*random.sample(ancindex, 1))
					fw.write(' '.join(tstr)+'\n')
				fw.close()

				nlines = subprocess.check_output(["/bin/sh", "-c", "wc -l %s; exit 0"%(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt',)])
				nlines = int(nlines.split()[0].strip())
				fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterbytes.txt','w')
				for i in range(nlines//2,0,-1):
					hap1 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i-1,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
					hap2 = str(subprocess.check_output(["/bin/sh", "-c", "sed '%dq;d' %s; exit 0"%(2*i,cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')]))
					fw.write('%s\n%s\n'%(hap1.strip(), hap2.strip()))
				fw.close()
				fo=open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterbytes.txt')
				fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt','w')
				for line in fo:
					fw.write(line[2:-4] +'\n')
				fw.close();fo.close()
				os.system('rm -r '+cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterbytes.txt')

			## LOTER to LAIT
			fo=open(cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt')
			fw = open(cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt','w')
			for line in fo:
				data=line.strip().split()
				fw.write(''.join(data)+'\n')
			fw.close()
			fo.close()
			##LAMPLD hap to LAIT (haplotypes in rows and SNPs in columns)
			os.system('perl ' +parameters['wrkdir']+'/franc_util/soft/lampld/standardizeOutput.pl lampld  '+str(len(parameters['anc_pop'])) +' '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt '+ cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterlait.txt')
			os.system('rm -rf '+cpath+parameters['outfile']+'.'+str(ch)+'.lampldhap.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloter.txt '+cpath+parameters['outfile']+'.'+str(ch)+'.chromopainterloterreverse.txt')
	else:
		print("\nLAIT standardizes ancestry in <= 3 way admixtures  ...\n")
		sys.exit(2)


if __name__=='__main__':
	pass

