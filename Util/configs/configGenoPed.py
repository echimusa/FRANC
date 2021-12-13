#!/usr/bin/env python

"""Converting OXFORD/EIGENSTRAT data format to PED. It creates the parameter file, run convertf and then convert the six column map file the PLINK map file with four column. The user can then use PLINK --make-bed and --chr  to create the desired PLINK bim/bed format and split to different chromosomes. The user should provide the working directory path, pop = list of populations to convert and the working directory contaninig softwares"""

from __future__ import print_function
import os,sys

wrkdir = "/home/ephie/"								#Directory path containing the OXFORD data format
pop = ["CEU","YRI","SIM2"]							#provide the ancestral population name prefix
soft = wrkdir +'soft'								#this "soft" is a directory inside the ``inputs" 

for p in pop:
	DataStruct = (wrkdir,p, wrkdir,p, wrkdir,p,wrkdir,p, wrkdir,p, wrkdir,p)
	StrStruct = open("GenoPed.txt").read()%DataStruct
	param = open(wrkdir+p+".parOXFORD_PED.GENO", 'w')
	print(StrStruct, file = param)
	param.close()

	os.system(soft+'convertf -p '+wrkdir+p+".parOXFORD_PED.GENO")	#Run CONVERTF

	fo = open(wrkdir+p+'.pedsnp')
	fw = open(wrkdir+p+".map", 'w')				#opening file to write and giving it the name
	for line in fo:
		tline=line.strip()
		if tline.startswith("#") or not tline: continue
		tline = [x.strip() for x in tline.split()]
		fw.write(tline[0] +'\t' + tline[1] + '\t0\t' + tline[3] + '\n')
	fo.close()
	fw.close()

