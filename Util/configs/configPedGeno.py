#!/usr/bin/env python

"""Converting PLINK PED data format to OXFORD/EIGENSTRAT. It creates the parameter file, and run convertf. The user should provide the working directory path, pop = list of populations to convert and the working directory contaninig softwares"""

from __future__ import print_function
import os,sys

infolder = "/home/ephie/"
pop = ["CEU","YRI","CHB","KHS","GIH"]
soft = infolder +'soft/'

for p in pop:
	##PLINK ped to OXFORD 
	DataStruct = (infolder, admix, admix,ch, infolder, admix, admix, ch, infolder, admix, admix, ch, infolder,admix,tool,admix,ch, infolder, admix, tool, admix, ch, infolder,admix, tool,admix)
	StrStruct = open("configPedGeno.txt").read()%DataStruct
	param = open(mpath+str(ch)+'/parPED_OXFORD.'+str(ch)+'.PED','w')
	print(StrStruct, file = param)
	param.close()

	os.system(soft+'convertf -p '+mpath+str(ch)+'/parPED_OXFORD.'+str(ch)+'.PED')
