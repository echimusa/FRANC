#!/usr/bin/python
# -*- coding: utf8 -*-

from __future__ import print_function
import os, sys, subprocess, fileinput, gzip, time, datetime
from convertall import *

__version__ = "v1.0"
__author__ = "E. Geza, N. Mulder, E.R Chimusa and G.K Mazandu"
__email__ = "kuzamunu@aims.ac.za and ephie@aims.ac.za"
__status__ = "Production"

try:
	from itertools import izip
except:
	izip = zip

def readuserinput(name_of_form):
	"""This function gets all the information provided in the parameter file (FrancPar.txt) to run franc. The form consists of all the information the user provides to run any of the nine local ancestry deconvolution tools. FrancPar.txt is a file accompanying franc with information in lines separated by = sign"""
	parameters = {}; Tools = ['lampld','winpop',"elai","supportmix", "pcadmix","chromopainter","rfmix","loter"]#tools analyzed
	cline = 0
	try:
		parfile = open(name_of_form)								#reading parameter file
		for lines in parfile:
			tline = lines.strip()
			cline += 1
			if not tline or tline.startswith('#'): continue			#ignore blank and comments
			elif tline.startswith('anc_pop'):
				anc_pop = tline.split('=')[1].strip()
				anc_pop = anc_pop.split('#')[0].strip().upper()		#ancestral populations file name prefix
				anc_pop = [a.strip() for a in anc_pop.split(':')]	#ancestral populations splitted by ":"
				parameters['anc_pop'] = anc_pop
		
			elif tline.startswith('anc_prop'):
				anc_prop = tline.split('=')[1].strip()
				anc_prop = anc_prop.split('#')[0].strip()
				anc_prop = [float(a.strip()) for a in anc_prop.split(',')]
				parameters['anc_prop'] = anc_prop
	
			elif tline.startswith('CHR'):
				CHR = tline.split('=')[1].strip()
				CHR = CHR.split('#')[0].strip()
				parameters['CHR'] = set([int(s.strip()) for s in CHR.split(',')])

			#Generations since admixture: SUPPORTMIX, WINPOP, RFMIX, ELAI
			elif tline.startswith('G'):										
				G = tline.split('=')[1].strip()
				G = G .split('#')[0].strip()
				parameters['G'] = G

			elif tline.startswith('recombrate'):
				recombrate = tline.split('=')[1].strip()
				recombrate = recombrate.split('#')[0].strip()
				parameters['recombrate'] = recombrate
			
			elif tline.startswith('offset'):
				offset= tline.split('=')[1].strip()
				offset = offset.split('#')[0].strip()
				offset = float(offset)
				parameters['offset'] = offset

			elif tline.startswith('ldcutoff'):
				ldcutoff = tline.split('=')[1].strip()
				ldcutoff = ldcutoff.split('#')[0].strip()
				ldcutoff = float(ldcutoff)
				parameters['ldcutoff'] = ldcutoff
		
			elif tline.startswith('phased'):
				phased = tline.split('=')[1].strip()
				phased = phased.split('#')[0].strip()
				parameters['phased'] = phased

			elif tline.startswith('PopPhased'):
				PopPhased = tline.split('=')[1].strip()
				PopPhased = PopPhased.split('#')[0].strip()
				parameters['PopPhased'] = PopPhased

			elif tline.startswith('parallel'):
				parallel = tline.split('=')[1].strip()
				parallel = parallel.split('#')[0].strip()
				parameters['parallel'] = parallel

			elif tline.startswith('w'):
				w = tline.split('=')[1].strip()
				w = w .split('#')[0].strip()
				w = int(w)
				parameters['w'] = w

			elif tline.startswith('resolve'):
				resolve = tline.split('=')[1].strip()
				resolve = resolve.split('#')[0].strip()
				parameters['resolve'] = resolve

			elif tline.startswith('model'):
				model = tline.split('=')[1].strip()
				model = model.split('#')[0].strip()
				parameters['model'] = model

	except:
		print ("\nPlease recheck your parameter file, and try again ...\n")
		sys.exit(2)
	
	return parameters 
	

class franc(object):
	def __init__(self, parameters):
		self.infolder = parameters['infolder']
		self.admix = parameters['admix']
		self.tool = parameters['tool']
		self.mode = parameters['mode']
		self.outfile = parameters['outfile']
		self.anc_pop = parameters['anc_pop']
		self.anc_prop = parameters['anc_prop']
		self.window = parameters['w']
		self.CHR = parameters['CHR']
		self.G = int(parameters['G'])
		self.recombrate = parameters['recombrate']
		self.offset = parameters['offset']
		self.ldcutoff = parameters['ldcutoff']
		self.phased = parameters['phased']
		self.PopPhased = parameters['PopPhased']
		self.parallel = parameters['parallel']
		self.resolve = parameters['resolve']
		self.model = parameters['model']
		self.outform = parameters['outformat']
		self.wrkdir = parameters['wrkdir']
		self.soft = self.wrkdir+'/franc_util/soft/'
		self.all = parameters


	"""#####********************************************************************************************************************###
	###															SUPPORTMIX														
	####************************************************************************************************************************"""
	def supportmix(self):
		mpath = self.infolder + '/' + self.admix+'/supportmix/'
		os.chdir(mpath)
		for ch in self.CHR:
			if not os.path.exists(mpath+str(ch)+'/data'+str(ch)+'/'):
				os.system('mkdir '+mpath+str(ch)+'/data'+str(ch))
			if not os.path.exists(mpath+str(ch)+'/data'+str(ch)+'/tpedData/'):
				os.system('mkdir '+mpath+str(ch)+'/data'+str(ch)+'/tpedData')
			if not os.path.exists(mpath+str(ch)+'/data'+str(ch)+'/hapmap2Subset/'):
				os.system('mkdir '+mpath+str(ch)+'/data'+str(ch)+'/hapmap2Subset')
			if not os.path.exists(mpath+str(ch)+'/output/'):
				os.system('mkdir '+mpath+str(ch)+'/output')

			##SUPPORTMIX Configuration file manipulation##If the window size in the parameter file =20, then we default the window size to 100SNPs for SUPPORTMIX
			if self.window == 20:
				w = 100
			else:
				w = self.window
			t = "supportmix"
			
			"""This part is a code to create the configuration file for SUPPORTMIX"""
			DataStruct = (self.infolder, self.admix, t, ch, ch, self.infolder, self.admix, t, ch, ch, self.admix,ch) 
			dyndata = ""
			n = len(self.anc_pop)
			for i in range(n):
				dyndata += "sample%d = %s/%s/%s/%d/data%d/tpedData/%s.%d.tped.gz\n"%(i+1, self.infolder, self.admix, t, ch, ch, self.anc_pop[i], ch)
			DataStruct += (dyndata, str(self.anc_pop)[1:-1], self.infolder, self.admix, t, ch, self.admix, ch, w, ch, ch, self.infolder, self.admix, t, ch, ch, self.G)
			StrStruct = open(self.wrkdir+"/franc_util/configs/configsvm.txt").read()%DataStruct
			param = open(mpath+str(ch)+'/svmConf_'+self.admix+'.'+str(ch)+'.cfg','wt')
			print(StrStruct, file = param)
			param.close()

			pop = self.anc_pop + [self.admix]
			for p in pop:
				##PLINK bed to VCF
				os.system(self.soft+'plink --bfile '+ self.infolder+'/'+p+'.'+str(ch)+' --chr '+ str(ch) +' --recode vcf --out '+mpath+p+'.'+str(ch))
				##Phase ancestral VCF with EAGLE
				os.system(self.soft+'eagle --vcf='+mpath+p+'.'+str(ch)+'.vcf --geneticMapFile='+self.infolder+'/genetic_map_phase.txt.gz --chrom='+str(ch)+' --outPrefix='+mpath+p+'.'+str(ch)+'.eagle --numThreads=4 --genoErrProb 0.003 --pbwtOnly 2>&1 | tee '+mpath+t+'_'+self.admix+'.'+str(ch)+'.log')
				#Unzipping VCF files
				os.system('gzip -d '+mpath+p+'.'+str(ch)+'.eagle.vcf.gz')
				##Phased VCF to PLINK bed
				os.system(self.soft+'plink --vcf '+mpath+p+'.'+str(ch)+'.eagle.vcf --double-id --vcf-require-gt --recode --make-bed --out '+mpath+p+'.'+str(ch)+'.eagle')
				##PLINK bed to PLINK tped
				os.system(self.soft+'plink --bfile '+ mpath+p+'.'+str(ch)+'.eagle --recode transpose --out '+mpath+p+'.'+str(ch)+'.eagle')
				##Move PLINK tped to right Dir
				os.system('mv '+mpath+p+'.'+str(ch)+'.eagle.tped '+mpath+str(ch)+'/data'+str(ch)+'/tpedData/'+p+'.'+str(ch)+'.tped') 
				##Convert to right PLINK tfam files for all populations (admixed and ancestral)
				fo = open(mpath+p+'.'+str(ch)+'.eagle.tfam')
				fw = open(mpath+str(ch)+'/data'+str(ch)+'/tpedData/'+p+'.'+str(ch)+'.tfam','w')
				for line in fo:
					data=line.split()
					fw.write(p + " " + data[1] + " " + data[2] + " " + data[3] + " " + data[4] + "\n")
				fw.close()
				fo.close()
				##zip SUPPORTMIX tped and tfam files
				os.system('gzip '+mpath+str(ch)+'/data'+str(ch)+'/tpedData/'+p+'.'+str(ch)+'.*')

			##Copy the genetic maps to respective folders
			os.system('cp '+self.infolder+'/genetic_map_chr'+str(ch)+'.txt '+mpath+str(ch)+'/data'+str(ch)+'/hapmap2Subset/genetic_map_chr'+str(ch)+'_txt')

			###**********************************************
			##RUNNING SUPPORTMIX
			##***********************************************
			os.chdir(self.soft)
			os.system('unzip -qo supportmix.zip')
			os.system('mv '+self.soft+'SupportMixDistribution '+self.soft+'supportmix')
			os.system(self.soft+'supportmix/Application/SupportMix -C '+mpath+str(ch)+'/svmConf_'+self.admix+'.'+str(ch)+'.cfg')

			##deleting unnecessary files
			for p in pop:
				os.system('rm -rf '+mpath+p+'.'+str(ch)+'.vcf '+mpath+p+'.'+str(ch)+'.eagle.vcf.gz '+mpath+p+'.'+str(ch)+'.eagle.vcf ' + mpath+p+'.'+str(ch)+'.eagle.bed '+ mpath+p+'.'+str(ch)+'.eagle.bim '+ mpath+p+'.'+str(ch)+'.eagle.fam '+ mpath+p+'.'+str(ch)+'.eagle.ped '+ mpath+p+'.'+str(ch)+'.eagle.map '+ mpath+p+'.'+str(ch)+'.eagle.log '+ mpath+p+'.'+str(ch)+'.eagle.nosex ')
		os.system('rm -rf '+self.soft+'supportmix')

	"""#####********************************************************************************************************************###
	###															WINPOP														
	#####******************************************************************************************************************"""

	def winpop(self):
		t = "winpop"
		os.chdir(self.infolder+'/'+self.admix+'/winpop')
		mpath = self.infolder+'/'+self.admix+'/'+t+'/'
		for ch in self.CHR:
			cpath = mpath + str(ch)+'/'

			"""ADMIXED Genotypes"""
			##PLINK bed to PLINK ped
			os.system(self.soft+'plink --bfile '+self.infolder+'/'+self.admix+'.'+str(ch)+ ' --recode --out '+ self.infolder+'/'+self.admix+'/'+self.admix+'.'+str(ch))

			##PLINK ped to OXFORD 
			DataStruct = (self.infolder,self.admix,self.admix,ch, self.infolder,self.admix,self.admix,ch, self.infolder,self.admix,self.admix,ch, self.infolder,self.admix,t,self.admix,ch, self.infolder,self.admix,t,self.admix,ch, self.infolder,self.admix, t,self.admix)
			StrStruct = open(self.wrkdir+"/franc_util/configs/configPedGeno.txt").read()%DataStruct
			param = open(cpath+'parPED_OXFORD.'+str(ch)+'.PED','w')
			print(StrStruct, file = param)
			param.close()

			os.system(self.soft+'convertf -p '+cpath+'parPED_OXFORD.'+str(ch)+'.PED')
			##Inverting admixed genotypes
			fo=mpath+self.admix+'.'+str(ch)+'.geno'
			fw=open(cpath+'geno.'+str(ch)+'.txt','w')
			rows={}
			for line in fileinput.input(fo):
				data=line.split()
				rows[fileinput.lineno()]=data[0]
				sample=len(data[0])
			for i in range(sample):
				for des in rows:
					if rows[des][i] in ["2","1","0"]:					#if genotypes aren't missing
						fw.writelines(rows[des][i] +" ")				#tab seperated genotypes for WINPOP
					else:
						fw.writelines("-1" +" ")						#missing values in WINPOP
				fw.write("\n")
			fw.close()
			##Ancestral populations
			for anc in self.anc_pop:
				##MAF
				os.system(self.soft+'plink --bfile '+self.infolder+'/'+anc+'.'+str(ch) +' --freq --out '+mpath+anc+'.'+str(ch))
				##Calculate the reference allele & writing output of each CHR
				fo = mpath+anc+'.'+str(ch)+'.frq'
				fw = open(mpath+str(ch)+'/'+anc+'.'+str(ch)+'.prob','w')
				for line in fileinput.input(fo):
					data = line.split()
					if fileinput.lineno()>1:
						fw.writelines(str(1.0-float(data[4]))+"\n")
				fw.close()

			##MANIPULATING POS FILE
			os.system("awk '{print $4}' "+ self.infolder+'/'+self.anc_pop[0]+'.'+str(ch)+'.bim > '+mpath+str(ch)+'/chr'+str(ch)+'.pos')
			##configuration file for WINPOP
			pfil = []
			for p in self.anc_pop:
				pfil.append(mpath+str(ch)+'/'+p+'.'+str(ch)+'.prob')
			F = ",".join(pfil)
			anc_prop = str(self.anc_prop)[1:-1]
			anc_prop = anc_prop.replace(' ','')
			DataStruct = (len(self.anc_pop),F,self.infolder,self.admix,t,ch,ch,self.infolder,self.admix,t,ch,ch,self.infolder, self.admix,t,ch,self.admix,ch,self.G,anc_prop)
			StrStruct = open(self.wrkdir+"/franc_util/configs/par_winpop.txt").read()%DataStruct
			param=open(mpath+str(ch)+"/config.txt","w")
			print(StrStruct, file = param)
			param.close()

			####************************************************************************
			####RUNNING WINPOP 
			####************************************************************************
			"""Checking inputs for WINPOP, required"""
			if len(self.anc_pop) != len(self.anc_prop):					#compare # of anc_pop to # of anc_prop
				print("\nPlease number of ancestral populations should tally the ancestral proportions to run in %s tool. .. number of ancestral population is %d, while proportions are provided for %d "%(t,len(self.anc_pop),len(self.anc_prop)))
				sys.exit(2)

			os.system('cp -r '+self.soft+'lamp/bin '+mpath+str(ch)+'/')
			os.chdir(mpath+str(ch)+'/')
			os.system(mpath+str(ch)+'/bin/lamp '+mpath+str(ch)+'/config.txt')

			##Removing unneccesary files
			for anc in self.anc_pop:
				os.system('rm -rf '+mpath+anc+'.'+str(ch)+'.frq '+cpath+anc+'.'+str(ch)+'.prob '+mpath+anc+'.'+str(ch)+'.log')
			os.system('rm -rf '+mpath+str(ch)+'/chr'+str(ch)+'.pos '+cpath+'geno.'+str(ch)+'.txt '+mpath+self.admix+'.'+str(ch)+'.geno '+cpath+'parPED_OXFORD.'+str(ch)+'.PED '+self.infolder +'/'+self.admix+'/'+self.admix+'.* '+cpath+"bin "+cpath+'frequency-estimates.*')


	"""#####********************************************************************************************************************###
	###															LAMPLD														
	####************************************************************************************************************************"""
	def lampld(self):
		t = 'lampld'
		os.chdir(self.infolder+'/'+self.admix+'/lampld')
		mpath = self.infolder+'/'+self.admix+'/lampld/'
		for ch in self.CHR:
			"""ADMIXED Genotypes for either LAMPLD or WINPOP"""
			cpath = self.infolder+'/'+self.admix+'/lampld/'+str(ch)+'/'
			##PLINK bed to PLINK ped
			os.system(self.soft+'plink --bfile '+self.infolder+'/'+self.admix+'.'+str(ch)+ ' --recode --out '+ self.infolder+'/'+self.admix+'/'+self.admix+'.'+str(ch))
			os.system('rm -rf '+self.infolder+'/'+self.admix+'/'+self.admix+'.'+str(ch)+'.log')

			##PLINK ped to OXFORD config file 
			DataStruct = (self.infolder,self.admix,self.admix,ch, self.infolder,self.admix,self.admix,ch, self.infolder,self.admix,self.admix,ch, self.infolder,self.admix,t,self.admix,ch, self.infolder,self.admix,t,self.admix,ch, self.infolder,self.admix, t,self.admix)
			StrStruct = open(self.wrkdir+"/franc_util/configs/configPedGeno.txt").read()%DataStruct
			param = open(cpath+'parPED_OXFORD.'+str(ch)+'.PED','w')
			print(StrStruct, file = param)
			param.close()

			os.system(self.soft+'convertf -p '+cpath+'parPED_OXFORD.'+str(ch)+'.PED')
			##Inverting admixed genotypes
			fo=mpath+self.admix+'.'+str(ch)+'.geno'
			fw=open(cpath+'geno.'+str(ch)+'.txt','w')
			rows={}
			for line in fileinput.input(fo):
				data=line.split()
				rows[fileinput.lineno()]=data[0]
				sample=len(data[0])
			for i in range(sample):
				for des in rows:
					if rows[des][i] in ["2","1","0"]:					#if genotypes aren't missing
						fw.writelines(rows[des][i])						#no space for LAMPLD genotypes
					else:
						fw.writelines("?")								#missing values in LAMPLD
				fw.write("\n")
			fw.close()

			##Ancestral populations phasing
			for anc in self.anc_pop:
				os.system(self.soft+'eagle --bfile='+self.infolder+'/'+anc+'.'+str(ch)+' --geneticMapFile='+self.infolder+'/genetic_map_phase.txt.gz --chrom='+str(ch) +' --outPrefix='+mpath+anc+'.'+str(ch)+'.eagle --numThreads=4 --genoErrProb 0.003 --pbwtOnly 2>&1 | tee '+mpath+'lampld.'+str(ch)+'.log')

				##Removing 1st 5 rows of shapeit/eagle phased data & transposing from SNPs/rows to individuals/rows 
				fo=gzip.open(mpath+anc+'.'+str(ch)+'.eagle.haps.gz')
				fw1=open(mpath+anc+'.'+str(ch)+'.haps','w')

				os.system('zcat  '+mpath+anc+'.'+str(ch)+".eagle.haps.gz | awk '{for (i=1; i<=NF-5; i++) $i = $(i+5); NF-=5; print}' > "+ mpath+anc+'.'+str(ch)+'.haps')
				fw2=open(cpath+anc+'.'+str(ch)+'.haps_ref', 'w')
				with open(mpath+anc+'.'+str(ch)+'.haps') as f:
					lis = [x.split() for x in f]
					for x in zip(*lis):								#transpose sequence of iterables:zip(*original_list)
						for y in x:
							fw2.write(y)							#fw2.write(y+' ')#if space
						fw2.write('\n')
				fw2.close()
				fo.close()

			##POSITION file 
			os.system('zcat '+mpath+self.anc_pop[0]+'.'+str(ch)+".eagle.haps.gz |  awk '{print $3}' >  " +cpath+'chr'+str(ch)+'.pos')

			####**************************************************************************************
			##RUNNING LAMPLD
			####**************************************************************************************
			num_anc = len(self.anc_pop)
			if (2 or 3 or 5) != num_anc: 
				print("\nNumber of Ancestral Populations is %d. LAMPLD can only run for 2, 3 and 5 way admixture ONLY, please choose another tool and try again ..."%(num_anc))
				sys.exit(2)
			else:
				os.system('cp -r '+self.soft+t+'/bin '+cpath)
				pfil = []
				for p in self.anc_pop:
					pfil.append(cpath+p+'.'+str(ch)+'.haps_ref')
				F = " ".join(pfil)
				os.system('cp -r '+self.soft+'lampld/run_LAMPLD'+str(num_anc)+'way.pl '+cpath)
				os.chdir(cpath)
				os.chmod(cpath+'run_LAMPLD'+str(num_anc)+'way.pl',0o777)
				##Changing file permission
				os.chmod(cpath+'bin/unolanc2way',0o777);os.chmod(cpath+'bin/triolanc',0o777);os.chmod(cpath+'bin/unolanc5way',0o777)
				os.chmod(cpath+'bin/haplanc',0o777);os.chmod(cpath+'bin/gettriophase',0o777);os.chmod(cpath+'bin/trainhmm',0o777)

				##Running LAMPLD for 2, 3,5 WAY ADMIXTURE
				os.system('perl '+cpath+'run_LAMPLD'+str(num_anc)+'way.pl '+cpath+'chr'+str(ch)+'.pos '+F+ ' '+cpath+'geno.'+str(ch)+'.txt '+cpath+self.outfile+'.'+str(ch)+'.out > '+cpath+self.outfile+'.'+str(ch)+'.log')

			for anc in self.anc_pop:
				os.system('rm -rf '+mpath+anc+'.'+str(ch)+'.eagle.* '+mpath+anc+'.'+str(ch)+'.haps '+cpath+anc+'.'+str(ch)+'.haps_ref '+cpath+p+'.'+str(ch)+'.log')
			os.system('rm -rf '+ self.infolder+'/'+self.admix+'/'+self.admix+'.'+str(ch)+'.ped '+ self.infolder+'/'+self.admix+'/'+self.admix+'.'+str(ch)+'.map '+mpath+self.admix+'.* '+cpath+'chr'+str(ch)+'.pos '+cpath+'geno.'+str(ch)+'.txt '+ cpath+'parPED_OXFORD.19.PED '+ mpath+'lampld.'+str(ch)+'.log '+cpath+'bin '+cpath+'run_LAMPLD'+str(num_anc)+'way.pl '+ cpath+self.outfile+'.lampld.'+str(ch)+'.log')



	"""#####****************************************************************************************************************###
	###							ELAI										
	####********************************************************************************************************************"""
	def elai(self):
		mpath = self.infolder+'/'+self.admix+'/elai/'
		for ch in self.CHR:
			cpath = self.infolder+'/'+self.admix+'/elai/'+str(ch)+'/'
			pop = self.anc_pop+[self.admix]
			for p in pop:
				if self.phased == 'YES':
					##PLINK bed to VCF
					os.system(self.soft+'plink --bfile '+self.infolder+'/'+p+'.'+str(ch)+' --chr '+str(ch)+' --recode vcf --out '+mpath+p+'.'+str(ch))
					##Phasing ALL pop with EAGLE to get phased VCF 
					os.system(self.soft+'eagle --vcf='+mpath+p+'.'+str(ch)+'.vcf --geneticMapFile='+self.infolder+'/genetic_map_phase.txt.gz --chrom='+str(ch)+' --outPrefix='+mpath+p+'.'+str(ch)+'.eagle --numThreads=4 --genoErrProb 0.003 --pbwtOnly 2>&1 | tee '+mpath+'elai.'+str(ch)+'.log')
					##Unzip VCFs
					os.system('gzip -d '+mpath+p+'.'+str(ch)+'.eagle.vcf.gz')

					##VCF to PLINK bed
					os.system(self.soft+'plink --vcf '+mpath+p+'.'+str(ch)+'.eagle.vcf --double-id --vcf-require-gt --recode --make-bed --out '+mpath+p+'.'+str(ch)+'.eagle')
					##PLINK bed to BIMBAM
					os.system(self.soft+'plink --bfile '+mpath+p+'.'+str(ch)+'.eagle --recode bimbam --out '+cpath+p+'.'+str(ch))
					##Add an "=" to the first row
					os.system('sed -i "1 s|$|=|" '+cpath+p+'.'+str(ch)+'.recode.geno.txt')
					os.system('rm -rf '+ mpath+p+'.* '+cpath+p+'.'+str(ch)+'.nosex')

				else:
					os.system(self.soft+'plink --bfile '+self.infolder+'/'+p+'.'+str(ch)+' --recode bimbam --out '+cpath+p+'.'+str(ch))
			##RUNNING ELAI
			anc_pop = self.anc_pop
			geno = [10+anc_pop.index(i) for i in anc_pop]
			os.chdir(self.soft)
			os.system('unzip -qo elai.zip')
			pfil = [self.soft+'elai/elai-lin -g ']
			for p,g in zip(anc_pop,geno):
				pops = cpath+p+'.'+str(ch)+'.recode.geno.txt','-p',str(g),'-g '
				a = ' '.join(pops)
				pfil.append(a)
			F = ' '.join(pfil)
			os.chdir(cpath)

			os.system(F+' '+cpath+self.admix+'.'+str(ch)+'.recode.geno.txt -p 1 -pos '+cpath+ self.anc_pop[0]+'.'+str(ch)+'.recode.pos.txt -o '+self.outfile+'.'+str(ch)+' -C '+str(len(anc_pop))+' -c '+str(5*len(anc_pop))+ ' -mg '+ str(self.G)+' > '+cpath+self.outfile+'.'+str(ch)+'.log')

			##Moving the output file to the right directory for further processing
			os.system('mv '+cpath+'output/* ' + cpath)
			##Removing unneceassary files
			for p in pop:
				os.system('rm -rf '+cpath+p+'*.recode.* '+cpath+p+'*log') 
		os.system('rm -rf '+self.soft+'elai '+ cpath+'output')


	"""#####***********************************************************************************************************###
	###										PCADMIX												
	#####****************************************************************************************************************"""
	def pcadmix(self):
		os.chdir(self.infolder+'/'+self.admix+'/pcadmix')
		mpath = self.infolder+'/'+self.admix+'/pcadmix/'
		for ch in self.CHR:
			cpath = mpath+str(ch)+'/'
			pop = [self.admix] + self.anc_pop
			pop_index = [pop.index(i) for i in pop]
			for p in pop:
				##PLINK bed to VCF
				os.system(self.soft+'plink --bfile '+ self.infolder+'/'+p+'.'+str(ch)+' --chr '+ str(ch) +' --biallelic-only strict --recode vcf --out '+mpath+p+'.'+str(ch))
				##VCF to EAGLE VCFs
				os.system(self.soft+'eagle --vcf='+mpath+p+'.'+str(ch)+'.vcf --geneticMapFile='+self.infolder+'/genetic_map_phase.txt.gz --chrom='+str(ch)+' --outPrefix='+mpath+p+'.'+str(ch)+'.eagle --numThreads=4 --genoErrProb 0.003 --pbwtOnly 2>&1 | tee '+mpath+self.outfile+'.'+str(ch)+'.log')
				##Unzip VCFs for use in BEAGLE conversion
				os.system('gzip -d '+mpath+p+'.'+str(ch)+'.eagle.vcf.gz')
				##VCF to PLINK bed
				os.system(self.soft+'plink --vcf '+mpath+p+'.'+str(ch)+'.eagle.vcf --double-id --vcf-require-gt --recode --make-bed --out '+mpath+p+'.'+str(ch))
				##EAGLE VCFs & PLINK ped to BEAGLE VCFs
				os.system('java -Xmx18204m -jar '+self.soft+'beagle.jar gt='+mpath+p+'.'+str(ch)+'.eagle.vcf ped='+mpath+p+'.'+str(ch)+'.ped out='+mpath+p+'.bgl.'+str(ch))
				##Convert the BEAGLE vcf format to the BEAGLE v3 format
				os.system('zcat '+mpath+p+'.bgl.'+str(ch)+'.vcf.gz | java -Xmx18204m -jar '+self.soft+'vcf2beagle.jar  "?" '+cpath+p+'.'+str(ch))
				##Unzip BEAGLE VCFs 
				os.system('gzip -d '+cpath+p+'.'+str(ch)+'.bgl.gz')

				###Renaming
				os.system('mv '+cpath+p+'.'+str(ch)+'.bgl '+mpath+str(ch)+'/'+p+'.'+str(ch)+'.txt')

			####**********************************************************************************************
			####RUNNING PCADMIX
			####**********************************************************************************************

			pfil=[self.soft+'PCAdmix3_linux -anc ']
			lab = ' '.join(self.anc_pop)
			for anc in self.anc_pop:
				pfil.append(mpath+str(ch)+'/'+anc+'.'+str(ch)+'.txt')
			F = ' '.join(pfil)
			os.chdir(mpath+str(ch)+'/')

			os.system(F+' -adm '+mpath+str(ch)+'/'+self.admix+'.'+str(ch)+'.txt -map '+mpath+self.anc_pop[0]+'.'+str(ch)+'.map -rho ' + self.infolder+'/genetic_map_chr'+str(ch)+'.txt -lab '+lab+' '+ self.admix +' -w '+str(self.window)+' -prune 0 -o '+mpath+str(ch)+'/'+self.outfile+'.'+str(ch)+' > '+cpath+self.outfile+'.'+str(ch)+'.log')

			##Removing unneccessary files
			for p in pop:
				os.system('rm -rf '+cpath+p+'.'+str(ch)+'.markers '+cpath+p+'.'+str(ch)+'.int '+mpath+p+'.'+str(ch)+'.* '+mpath+p+'.bgl.'+str(ch)+'.vcf.gz '+mpath+p+'.bgl.'+str(ch)+'.log '+cpath+p+'.'+str(ch)+'.txt '+cpath+p+'.'+str(ch)+'.bgl.gz ')
			os.system('rm -rf '+mpath+self.outfile+'.'+str(ch)+'.log')

	"""#####******************************************************************************************************************###
	###															RFMIX														
	#####*******************************************************************************************************************###"""
	def rfmix(self):
		os.chdir(self.infolder+'/'+self.admix+'/rfmix')
		t = "rfmix"
		mpath = self.infolder+'/'+self.admix+'/'+t+'/'
		for ch in self.CHR:
			cpath = mpath+str(ch)+'/'
			pop = [self.admix] + self.anc_pop
			pop_index = [pop.index(i) for i in pop]
			for p in pop:
				##PLINK bed to VCF
				os.system(self.soft+'plink --bfile '+ self.infolder+'/'+p+'.'+str(ch)+' --chr '+ str(ch) +' --biallelic-only strict --recode vcf --out '+mpath+p+'.'+str(ch))
				##Phasing: VCF to EAGLE VCFs
				os.system(self.soft+'eagle --vcf='+mpath+p+'.'+str(ch)+'.vcf --geneticMapFile='+self.infolder+'/genetic_map_phase.txt.gz --chrom='+str(ch)+' --outPrefix='+mpath+p+'.'+str(ch)+'.eagle --numThreads=4 --genoErrProb 0.003 --pbwtOnly 2>&1 | tee '+mpath+'rfmix_'+self.admix+'.'+str(ch)+'.log')
				##Unzip VCFs for use in BEAGLE conversion
				os.system('gzip -d '+mpath+p+'.'+str(ch)+'.eagle.vcf.gz')
				##VCF to PLINK bed
				os.system(self.soft+'plink --vcf '+mpath+p+'.'+str(ch)+'.eagle.vcf --double-id --vcf-require-gt --recode --make-bed --out '+mpath+p+'.'+str(ch))
				##EAGLE VCFs & PLINK ped to BEAGLE VCFs
				os.system('java -Xmx18204m -jar '+self.soft+'beagle.jar gt='+mpath+p+'.'+str(ch)+'.eagle.vcf ped='+mpath+p+'.'+str(ch)+'.ped out='+mpath+p+'.bgl.'+str(ch))
				##Convert the BEAGLE vcf format to the BEAGLE v3 format
				os.system('zcat '+mpath+p+'.bgl.'+str(ch)+'.vcf.gz | java -Xmx18204m -jar '+self.soft+'vcf2beagle.jar  "?" '+cpath+ p+'.'+str(ch))
				##Unzip BEAGLE VCFs 
				os.system('gzip -d '+cpath+p+'.'+str(ch)+'.bgl.gz')
				##Renaming
				os.system('mv '+cpath+p+'.'+str(ch)+'.bgl '+cpath+p+'chr'+str(ch)+'_bgl.phased')

			##"""MarkersLocation"""
			os.system('tail -n +11 '+ mpath+self.admix+'.'+str(ch)+'.eagle.vcf > '+cpath+self.admix+'.'+str(ch)+'.markers.vcf')
			fo=open(cpath+self.admix+'.'+str(ch)+'.markers.vcf')
			fw =open(cpath+'markerLocationsChr'+str(ch),'wt')
			for line in fo:
				data = line.split()
				fw.writelines(str(float(data[1])/1000000.0)+"\n")
			fw.close()
			fo.close()

			##Creating the R code to process alleles file
			DataStruct = (ch, ch, self.admix, self.infolder, self.admix, t, ch, self.admix, ch, self.admix, self.admix, self.admix, self.admix, self.admix, self.admix)
			pdata = ""
			n = len(self.anc_pop)
			for i in range(n):
				pdata += '\n ##Ancestral population %s\n'%(self.anc_pop[i],)
				pdata += '%s=read.table(paste("%s/%s/%s/%d/%schr%d_bgl.phased",sep=""), header=T, stringsAsFactors=F)\n'%(self.anc_pop[i], self.infolder, self.admix, t, ch, self.anc_pop[i], ch)
				pdata += '%s= %s[, -c(1:2)]\n'%(self.anc_pop[i],self.anc_pop[i])
				pdata += ' %sT = t(%s)\n'%(self.anc_pop[i],self.anc_pop[i])
				pdata += ' num%s = nrow(%sT)\n'%(self.anc_pop[i],self.anc_pop[i])

			DataStruct += (pdata,)

			m= ['%sT'%(self.anc_pop[i],) for i in range(n)]
			DataStruct += (self.admix, ",".join(m), self.infolder, self.admix, t, ch, self.admix, ch)

			c = ["rep("+str(pi)+",num"+i+")" for (i,pi) in zip(pop,pop_index)]

			DataStruct += ('classes = c('+','.join(c)+')\n', self.infolder, self.admix, t)
			StrStruct = open(self.wrkdir+"/franc_util/configs/configRfmix.txt").read()%DataStruct

			param = open(cpath+'ConvertToRFMixFormat.R','wt')
			print(StrStruct, file=param)
			param.close()

			##Convert the BEAGLE file format to RFMIX
			os.system('Rscript '+cpath+'ConvertToRFMixFormat.R')

			os.chdir(self.soft)
			os.system('unzip -qo rfmix.zip')
			os.system('mv -f '+self.soft+'RFMix_v1.5.4 '+self.soft+'rfmix')
			os.system('rm -rf '+self.soft+'rfmix/TestData '+self.soft+'rfmix/RFMix_v1.5.4 ')

			####***************************************************************************************
			####				RUNNING RFMIX
			####***************************************************************************************

			if type(self.window) is int:
				w=(self.window)/100.0
			else:
				w=self.window
			if self.PopPhased == 'YES':
				##*******PopPhased**************
				os.chdir(self.soft+'rfmix/PopPhased')
				os.system('g++ -Wall -O3 -ftree-vectorize -fopenmp main.cpp getdata.cpp randomforest.cpp crfviterbi.cpp windowtosnp.cpp -o RFMix_PopPhased')
				os.chdir(self.soft+'rfmix/')
				rfmix_starttime = datetime.datetime.today()
				print("\n\nRFMIX PopPhased for %s population,a %s way admixture on chromosome %s started at %s"%(self.admix, str(len(self.anc_pop)), str(ch),str(rfmix_starttime)),file=open(self.infolder+'/ToolAdmixTimeRecord.txt','a'))
				os.system('python RunRFMix.py PopPhased '+cpath+self.admix+'_alleles.'+str(ch)+'.txt '+mpath+'classes.txt '+cpath+'markerLocationsChr'+str(ch)+' -w '+str(w) +' -n 5 -G '+str(self.G)+' -o '+cpath+self.outfile+'.'+str(ch)+' > '+cpath+self.outfile+'.'+str(ch)+'.log')

			else:
				##FOR TrioPhased **************
				os.chdir(self.soft+'rfmix/TrioPhased')
				os.system('g++ -Wall -O3 -ftree-vectorize -fopenmp main.cpp getdata.cpp randomforest.cpp crfviterbi.cpp windowtosnp.cpp -o RFMix_TrioPhased')
				os.chdir(self.soft+'rfmix/')

				os.system('python RunRFMix.py TrioPhased '+cpath+self.admix+'_alleles.'+str(ch)+'.txt '+mpath+'classes.txt '+cpath+'markerLocationsChr'+str(ch)+' -w '+str(w) +' -n 5 -G '+str(self.G)+' -o '+cpath+self.outfile+'.'+str(ch)+' > '+cpath+self.outfile+'.'+str(ch)+'.log')

			##Removing unneccessary files
			for p in pop:
				os.system('rm -rf ' + cpath+p+'chr19_bgl.phased '+cpath+p+'.'+str(ch)+'.markers '+cpath+p+'.'+str(ch)+'.int '+ mpath+p+'.* ' + cpath+p+'.'+str(ch)+'.bgl.gz')
			os.system('rm -rf '+ mpath+'rfmix_'+self.admix+'.'+str(ch)+'.log '+ cpath+'markerLocationsChr'+str(ch) +' ' +cpath+'ConvertToRFMixFormat.R '+cpath+self.admix+'.'+str(ch)+'.markers.vcf '+cpath+self.admix+'_alleles.'+str(ch)+'.txt '+cpath+self.outfile+'.'+str(ch)+'.log '+mpath+'classes.txt') 
		os.system('rm -rf '+self.soft+'rfmix '+self.soft+'__MACOSX')


	"""#####*********************************************************************************************************###
	###											MULTIMIX								
	#####***********************************************************************************************************###"""
	def multimix(self):
		t = "multimix"
		os.chdir(self.infolder+'/'+self.admix+'/multimix')
		mpath = self.infolder+'/'+self.admix+'/'+t+'/'
		pop = [self.admix] + self.anc_pop
		for ch in self.CHR:
			cpath = mpath + str(ch)+'/'
			dpath = cpath+'data/'
			if not os.path.exists(dpath):
				os.system('mkdir '+dpath)
			if not os.path.exists(dpath+'Samples/'):
				os.system('mkdir '+dpath+'Samples/')
			if not os.path.exists(dpath+'haplotypes/'):
				os.system('mkdir '+dpath+'haplotypes/')
			for anc in self.anc_pop:
				if not os.path.exists(dpath+'haplotypes/'+anc):
					os.system('mkdir '+dpath+'haplotypes/'+anc)
			if not os.path.exists(dpath+'haplotypes/genetic_maps'):
				os.system('mkdir '+dpath+'haplotypes/genetic_maps')
			if not os.path.exists(dpath+'haplotypes/legend_files'):
				os.system('mkdir '+dpath+'haplotypes/legend_files')

			pfil = []
			for p in self.anc_pop:
				pfil.append(p+'/')
			F = " ".join(pfil)
			pop = self.anc_pop +[self.admix]
			##Genetic mapfile
			os.system('cp '+self.infolder+'/genetic_map_chr'+str(ch)+'.txt '+dpath+'haplotypes/genetic_maps/chr'+str(ch)+'.map')

			##ALL PHASED, admix & anc_pop
			if not os.path.exists(dpath+'Samples/haplos/'):
				os.system('mkdir '+dpath+'Samples/haplos/')
			##Phasing admix & anc_pop
			for p in pop:
				os.system(self.soft+'eagle --bfile='+self.infolder+'/'+p+'.'+str(ch)+' --geneticMapFile='+self.infolder+'/genetic_map_phase.txt.gz --chrom='+str(ch) +' --outPrefix='+mpath+p+'.'+str(ch)+'.eagle --numThreads=4 --genoErrProb 0.003 --pbwtOnly 2>&1 | tee '+mpath+p+'_phase.'+str(ch)+'.log')
			##Phased Samples to MULTIMIX sample files
			mapp={'0':'0 0 1','1':'1 0 0'}													#'?':'9 9 9', '9':'9 9 9', 'NA':'9 9 9'}
			fo=gzip.open(mpath+self.admix+'.'+str(ch)+'.eagle.haps.gz' ,'r')				#phased haplotypes
			fw=open(dpath+'Samples/haplos/haplos_chr'+str(ch),'w')							#multimix file wth maped haps 
			for line in fo:
				tline=line.strip().split()
				tlin=[]
				for lin in tline[1:5]:
					tlin.append(str(lin))
				fw.write(' '.join(tlin))
				for g in tline[5:]:
					fw.write(' ' + ' '.join(map(lambda x: mapp[x], g.decode('UTF-8'))))
				fw.write('\n')
			fw.close()
			fo.close()
		
			##legend file from admixed phased data
			os.system('zcat '+mpath+self.anc_pop[0]+'.'+str(ch)+".eagle.haps.gz | awk '{print $2,$3,$4,$5 }' > "+dpath+'haplotypes/legend_files/chr'+str(ch)+'.legend')
			os.system('sed -i "1i rs position a0  a1" '+dpath+'haplotypes/legend_files/chr'+str(ch)+'.legend')

			##MULTIMIX ancestral haplotypes
			for anc in self.anc_pop:
				os.system('zcat  '+mpath+anc+'.'+str(ch)+".eagle.haps.gz | awk '{for (i=1; i<=NF-5; i++) $i = $(i+5); NF-=5; print}' > "+dpath+'haplotypes/'+anc+'/chr'+str(ch)+'.haps')

			##Delete unnecessary files
			for p in pop:
				os.system('rm -rf '+mpath+p+'.'+str(ch)+'.*')

			###***************************************************************************************
			###RUNNING MULTIMIX
			####**************************************************************************************

			os.chdir(cpath)
			if not os.path.exists(cpath+'out_'+self.model):
				os.system('mkdir '+cpath+'out_'+self.model)


			if self.resolve == "YES":
				if len(self.anc_pop)==2:
					os.system(self.soft+'MULTIMIX_'+self.model+'_v1.1.0  -dir ./data/ -o ./out_'+self.model+'/ -sourcepops '+F+'  -misfit 0.950 0.05 0.05 0.950 -chr.list '+str(ch)+ ' -M -resolve > '+cpath+'out_'+self.model+'/logfile.'+str(ch))

				elif len(self.anc_pop)==3:
					os.system(self.soft+'MULTIMIX_'+self.model+'_v1.1.0  -dir ./data/ -o ./out_'+self.model+'/ -sourcepops '+F+' -misfit 0.95 0.025 0.025 0.025 0.95 0.025 0.025 0.025 0.95 -chr.list '+str(ch)+ ' -lambda 0.1 -M -resolve > '+mpath+'out_'+self.model+'/logfile.'+str(ch))

				elif len(self.anc_pop)==5:
					os.system(self.soft+'MULTIMIX_'+self.model+'_v1.1.0  -dir ./data/ -o ./out_'+self.model+'/ -sourcepops '+F+' -misfit 0.950 0.01250 0.01250 0.01250 0.01250 0.01250 0.950 0.01250 0.01250 0.01250 0.01250 0.01250 0.950 0.01250 0.01250 0.01250 0.01250 0.01250 0.950 0.01250 0.01250 0.01250 0.01250 0.01250 0.950 -chr.list '+str(ch)+' -lambda 0.1 -M -resolve > '+cpath+'out_'+self.model+'/logfile.'+str(ch))
				else:
					os.system(self.soft+'MULTIMIX_'+self.model+'_v1.1.0  -dir ./data/ -o ./out_'+self.model+'/ -sourcepops '+F+' -misfit 0.95 0.0167 0.0167 0.0166 0.0167 0.95 0.0167 0.0166 0.0167 0.0167 0.95 0.0166 0.0167 0.0167 0.0166 0.95 -chr.list '+str(ch)+ ' -lambda 0.1 -M -resolve > '+cpath+'out_'+self.model+'/logfile.'+str(ch))

				##RESOLVING
				os.system(self.soft+'MULTIMIX_resolve_v1.1.0 -dir ./data/ -i ./out_'+self.model+'/ -o ./out_'+self.model+'/ -sourcepops '+F+' -chr.list '+str(ch)+ ' > '+ cpath+'out_'+self.model+'/logfile.'+str(ch))

			##MULTIMIX WITHOUT RESOLVING
			else:
				if len(self.anc_pop)==2:
					os.system(self.soft+'MULTIMIX_'+self.model+'_v1.1.0  -dir ./data/ -o ./out_'+self.model+'/ -sourcepops '+F+'  -misfit 0.950 0.05 0.05 0.950 -chr.list '+str(ch)+ ' > '+cpath+'out_'+self.model+'/logfile.'+str(ch))

				elif len(self.anc_pop)==3:
					os.system(self.soft+'MULTIMIX_'+self.model+'_v1.1.0  -dir ./data/ -o ./out_'+self.model+'/ -sourcepops '+F+' -misfit 0.95 0.025 0.025 0.025 0.95 0.025 0.025 0.025 0.95 -chr.list '+str(ch)+ ' -lambda 0.1 > '+mpath+'out_'+self.model+'/logfile.'+str(ch))

				elif len(self.anc_pop)==5:
					os.system(self.soft+'MULTIMIX_'+self.model+'_v1.1.0  -dir ./data/ -o ./out_'+self.model+'/ -sourcepops '+F+' -misfit 0.950 0.01250 0.01250 0.01250 0.01250 0.01250 0.950 0.01250 0.01250 0.01250 0.01250 0.01250 0.950 0.01250 0.01250 0.01250 0.01250 0.01250 0.950 0.01250 0.01250 0.01250 0.01250 0.01250 0.950 -chr.list '+str(ch)+' -lambda 0.1 -M -resolve > '+cpath+'out_'+self.model+'/logfile.'+str(ch))
				else:
					os.system(self.soft+'MULTIMIX_'+self.model+'_v1.1.0  -dir ./data/ -o ./out_'+self.model+'/ -sourcepops '+F+' -misfit 0.95 0.0167 0.0167 0.0166 0.0167 0.95 0.0167 0.0166 0.0167 0.0167 0.95 0.0166 0.0167 0.0167 0.0166 0.95 -chr.list '+str(ch)+ ' -lambda 0.1 > '+cpath+'out_'+self.model+'/logfile.'+str(ch))

			##Deleting unnecessary files
			os.system('rm -rf '+cpath+'data '+mpath+'*.log')


	"""#####**************************************************************************************************###
	###															CHROMOPAINTER							
	#####*****************************************************************************************************###"""
	def chromopainter(self):
		os.chdir(self.infolder+'/'+self.admix+'/chromopainter')
		t="chromopainter"
		mpath = self.infolder+'/'+self.admix+'/'+t+'/'
		pop = self.anc_pop + [self.admix]
		lcmds = []
		for ch in self.CHR:
			cpath = mpath+str(ch)+'/'
			##coping the files from working dir to chromopainter directories
			for p in pop:
				os.system('cp '+self.infolder+'/'+p+'.'+str(ch)+'* '+ mpath)
			##Merging all populations together
			for i in range(len(pop)-1):
				n1 = ''.join(pop[0:i+1])
				n2 = pop[i+1]
				n1n2 = ''.join(pop[0:i+2])
				os.system(self.soft+'plink --bfile '+ mpath+n1+'.'+str(ch)+' --bmerge ' + mpath+n2 +'.'+str(ch)+'.bed '+ mpath+n2 +'.'+str(ch)+'.bim '+ mpath+n2 +'.'+str(ch)+'.fam --recode --make-bed --out '+mpath+n1n2+'.'+str(ch))

			##Creating the file mantaining sample ordering 
			for p in pop:
				os.system('cat '+mpath+p+'.'+str(ch)+'.fam >> '+mpath+'chrohaplo'+str(ch)+'.sample')
			##Order the merged files
			os.system(self.soft+'plink --bfile '+mpath+n1n2+'.'+str(ch)+' --indiv-sort f '+mpath+'chrohaplo'+str(ch)+'.sample --make-bed --out '+mpath+'DonorRecipient.'+str(ch))
			os.system('rm -rf '+mpath+'chrohaplo'+str(ch)+'.sample')
			##Phasing the merged file
			os.system(self.soft+'eagle --bfile='+mpath+'DonorRecipient.'+str(ch)+' --geneticMapFile='+self.infolder+'/genetic_map_phase.txt.gz --chrom='+str(ch) +' --outPrefix='+mpath+'DonorRecipient.'+str(ch)+' --numThreads=4 --genoErrProb 0.003 --pbwtOnly 2>&1 | tee '+mpath+self.admix+'_phase.'+str(ch)+'.log')
			##Unzipping
			os.system('gzip -d '+mpath+'DonorRecipient.'+str(ch)+'.haps.gz')
			##Phased haps to PHASE v1 format
			os.system('perl '+self.soft+'impute2chromopainter.pl -v1 -f '+mpath+'DonorRecipient.'+str(ch)+'.haps '+mpath+'Recomrates.'+str(ch))
			os.system('tail -n +2 '+self.infolder+'/genetic_map_chr'+str(ch)+'.txt > '+self.infolder+'/'+self.admix+'/gen_maprem.'+str(ch)+'.txt')
			fo=open(self.infolder+'/'+self.admix+'/gen_maprem.'+str(ch)+'.txt')
			fw = open(mpath+'gen_map.'+str(ch)+'.txt','w')
			fw.write('Chromosome\t'+'Position(bp)\t'+'Rate(cM/Mb)\t'+ 'Map(cM)')
			for line in fo:
				fw.write('chr'+str(ch)+ '\t'+ line)
			fw.close()
			fo.close

			os.system('perl '+self.soft+'convertrecfile.pl -M hapmap '+mpath+'Recomrates.'+str(ch)+'.phase '+mpath+'gen_map.'+str(ch)+'.txt ' +cpath+'chr'+str(ch)+'.recomrates > '+cpath+'chr'+str(ch)+'.recolog')
			##Manipulate the .phase files for all anc_pop (donor) and admix (recipient) pops 
			os.system('perl '+self.soft+'impute2chromopainter.pl '+mpath+'DonorRecipient.'+str(ch)+'.haps ' +cpath+'haplotype_infile.'+str(ch))
			###Label_infile
			fw = open(mpath+'label_infile','w')
			for p in pop:
				fo=open(mpath+p+'.'+str(ch)+'.fam')
				for line in fo:
					data = line.split()
					fw.write(data[1] + '\t' + p +'\t1\n')
				fo.close()
			fw.close()

			##Population infile
			fw=open(mpath+'popn_list_infile','w')
			for p in pop:
				if p == self.admix:
					fw.write(p+ ' ' +'R'+'\n')
				else:
					fw.write(p + ' '+'D'+'\n')
			fw.close

			####********************************************************************************************
			##						RUNNING CHROMOPAINTER	
			####********************************************************************************************

			strcmd = self.soft+'ChromoPainterv2 -g '+cpath+'haplotype_infile.'+str(ch)+'.phase -r '+cpath+'chr'+str(ch)+'.recomrates -t '+mpath+'label_infile -f '+mpath+'popn_list_infile 0 0 -ip -b -o '+cpath+self.outfile+'.'+str(ch)+ ' > '+cpath+'logfile.'+str(ch)
			lcmds.append(strcmd)
			##Deleting unneccessary files
			"""Don't delete ****cpath+'chr'+str(ch)+'.recomrates**** in chromosome dir, we need it for standardizing outputs """
			for i in range(len(pop)-1):
				n1 = ''.join(pop[0:i+1])
				n2 = pop[i+1]
				n1n2 = ''.join(pop[0:i+2])
				os.system('rm -rf ' +mpath+n1n2+'.'+str(ch)+'.ped '+mpath+n1n2+'.'+str(ch)+'.map '+mpath+n1n2+'.'+str(ch)+'.bed '+mpath+n1n2+'.'+str(ch)+'.bim '+mpath+n1n2+'.'+str(ch)+'.fam ')
			os.system('rm -rf '+mpath+'DonorRecipient.'+str(ch)+'.* '+self.infolder+'/'+self.admix+'/gen_maprem.'+str(ch)+'.txt '+mpath+'Recomrates.'+str(ch)+'.phase '+mpath+'gen_map.'+str(ch)+'.txt ')
			for p in pop:
				os.system('rm -rf '+ mpath+p+'.'+str(ch)+'.* ')
		return lcmds

	"""#####******************************************************************************************************************###
	###															LOTER														
	#####*******************************************************************************************************************###"""
	def loter(self):
		os.chdir(self.infolder+'/'+self.admix+'/loter')
		try:
			import allel
		except ImportError:
			raise ImportError("The library allel under Scikit-allel is required for running Loter. \nPlease, install it and try again ...")
		try:
			import numpy as np
		except ImportError:
			raise ImportError("The library Numpy is required for some numerical processes. \nPlease, install it and try again ...")
		try:
			import sklearn
		except ImportError:
			raise ImportError("The library sklearn is required for running Loter. \nPlease, install it and try again ...")
		try:
			import matplotlib
			matplotlib.use('Agg')
			import matplotlib.pyplot as plt
		except ImportError:
			raise ImportError("The library matplotlib is required for plotting. \nPlease, install it and try again ...")

		mpath = self.infolder+'/'+self.admix+'/loter/'
		pop = self.anc_pop + [self.admix]
		for ch in self.CHR:
			cpath = mpath + str(ch) + '/'
			def vcf2npy(vcfpath):
				callset = allel.read_vcf(vcfpath)
				haplotypes_1 = callset['calldata/GT'][:,:,0]
				haplotypes_2 = callset['calldata/GT'][:,:,1]

				m, n = haplotypes_1.shape
				mat_haplo = np.empty((2*n, m))
				mat_haplo[::2] = haplotypes_1.T
				mat_haplo[1::2] = haplotypes_2.T
				return mat_haplo.astype(np.uint8)
			for p in pop:
				##PLINK bed to VCF
				os.system(self.soft+'plink --bfile '+ self.infolder+'/'+p+'.'+str(ch)+' --chr '+ str(ch) +' --biallelic-only strict --recode vcf --out '+mpath+p+'.'+str(ch))
				##Phase with EAGLE to VCFs
				os.system(self.soft+'eagle --vcf='+mpath+p+'.'+str(ch)+'.vcf --geneticMapFile='+self.infolder+'/genetic_map_phase.txt.gz --chrom='+str(ch)+' --outPrefix='+cpath+p+'.'+str(ch)+'.eagle --numThreads=4 --genoErrProb 0.003 --pbwtOnly 2>&1 | tee '+mpath+'loter_'+self.admix+'.'+str(ch)+'.log')
				##Unzip VCFs
				os.system('gzip -d '+cpath+p+'.'+str(ch)+'.eagle.vcf.gz')

			os.chdir(self.soft)
			os.system('unzip -qo loter.zip')
			os.chdir(self.soft+'Loter-master')
			if self.parallel == 'YES':
				os.system('make ')
			else:
				os.system('make no_omp=1')

			####********************************************************************************************
			####Running LOTER
			####********************************************************************************************

			sys.path.append(self.soft+'Loter-master/python-package')
			import loter.locanc.local_ancestry as lc

			pfile=[]
			for anc in self.anc_pop:
				anc = vcf2npy(cpath+anc+'.'+str(ch)+'.eagle.vcf')
				pfile.append(anc)
			##Admixed haplotypes
			admx = vcf2npy(cpath+self.admix+'.'+str(ch)+'.eagle.vcf')

			##Putting the matrices of anc_pop & admix pop together & running LOTER
			adm = lc.loter_smooth(pfile, admx)
			res_impute, res_no_impute = lc.loter_local_ancestry(pfile, admx)

			##Ploting the results
			plt.imshow(res_no_impute[0], interpolation='nearest', aspect='auto')
			plt.colorbar()
			plt.savefig(cpath+self.outfile+'.'+str(ch)+'.png', bbox_inches='tight')
			plt.close()
			##Saving the results
			np.savetxt(cpath+self.outfile+'.'+str(ch)+'.ancestry.txt', res_no_impute[0], fmt="%i ")		#local anc
			np.savetxt(cpath+self.outfile+'.'+str(ch)+'.gen_anc.txt', res_impute, fmt="%i ")			#genotypic anc
			np.savetxt(cpath+self.outfile+'.'+str(ch)+'.NoTimeAncPicked.txt', res_no_impute[1], fmt="%i")## of time anc is picked

			for p in pop:
				os.system('rm -rf '+ cpath+p+'.'+str(ch)+'.eagle.vcf.gz '+cpath+p+'.'+str(ch)+'.eagle.vcf '+mpath+p+'.'+str(ch)+'.vcf '+ mpath+p+'.'+str(ch)+'.log ')
		os.system('rm -rf '+self.soft+'Loter-master')


	def run(self):
		##Testing possibility of conversion
		t1 = set(["pcadmix","supportmix","chromopainter","lampld","multimix","rfmix","loter","elai"])	##convertible to outform1
		outform1 = set(["rfmix","loter","winpop","lait"])					#convertible to from t1
		t2 = set(["winpop"])
		outform2 = set(["lait"])
		ConvertDict = {"pcadmix": 3, "supportmix":4, "chromopainter": 5, "lampld":7, "elai":8, "multimix":9, "rfmix":1, "loter":2, "winpop":6, "lait":10}
		for t in self.tool:
			if t=='lampld': self.lampld()
			elif t=='winpop': self.winpop()
			elif t=='elai': self.elai()
			elif t=='supportmix': self.supportmix()
			elif t=='pcadmix': self.pcadmix()
			elif t=='multimix': self.multimix()
			elif t=='chromopainter': 
				lcmds = self.chromopainter()
				for strcmd in lcmds:
					os.system("bash -c \"%s\""%(strcmd,))
				for ch in self.CHR:
					cpath = self.infolder+'/'+self.admix+'/chromopainter/'+str(ch)+'/'
					os.system('rm -r '+cpath+'chr'+str(ch)+'.recomrates '+cpath+'chr'+str(ch)+'.recolog '+cpath+'haplotype_infile.'+str(ch)+'.phase '+cpath+self.outfile+'.'+str(ch)+'.chunkcounts.out '+cpath+self.outfile+'.'+str(ch)+'.chunklengths.out '+cpath+self.outfile+'.'+str(ch)+'.EMprobs.out '+cpath+self.outfile+'.'+str(ch)+'.mutationprobs.out '+cpath+self.outfile+'.'+str(ch)+'.regionchunkcounts.out '+cpath+self.outfile+'.'+str(ch)+'.regionsquaredchunkcounts.out')
				os.system('rm -rf '+self.infolder+'/'+self.admix+'/chromopainter/*.log '+self.infolder+'/'+self.admix+'/chromopainter/label_infile '+self.infolder+'/'+self.admix+'/chromopainter/popn_list_infile')

			elif t=='rfmix': self.rfmix()
			elif t=='loter': self.loter()

			##No conversion needed
			if t==self.outform:
				pass
			##Convert by calling relevant conversion function, if possible
			elif (t in t1 and self.outform in outform1) or (t in t2 and self.outform in outform2):
				runcvt = eval("convert_%d%d"%(ConvertDict[t], ConvertDict[self.outform]))
				runcvt(self.all)
			else:
				pass


def is_valid_file(parser, arg):
	"""
	Check if a valid path to a file has been passed.

	Parameters
	parser, file path
	return the file path or fail with an error"""
	arg = os.path.abspath(arg)
	return arg if os.path.exists(arg) else parser.error("The file %s does not exist!" % arg)

def main():
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(description = __doc__, formatter_class = ArgumentDefaultsHelpFormatter)
	parser.add_argument("-p", "--par", dest="parameter", type = lambda x: is_valid_file(parser, x), help="Parameter input file containing information to be processed ", required = True, metavar="FILE")
	parser.add_argument("-d", "--dir", dest="infolder", type = lambda x: is_valid_file(parser, x), default = os.getcwd(), help="Working directory ", metavar="FILE")
	parser.add_argument("-t", "--tool", dest="tool", type = str, help="Tool/s to be run ", required = True, metavar="FILE")
	parser.add_argument("-a", "--admix", dest="admix", type = str, help="Admixed population under analysis, name prefix ", required = True, metavar="FILE")
	parser.add_argument("-o", "--out", dest="outfile", type = str, help="Output file prefix ", metavar="FILE")
	parser.add_argument("-m", "--mode", dest="mode", type = str, default="local", help="Indicate whether you run it from local or server ")
	parser.add_argument("-f", "--outformat", dest="outformat", type = str, default="tool", help="Output format", metavar="FILE")
	args = parser.parse_args()

	## Providing system and platform requirements
	print('\n'+74*'*'+'\n*FRANC %s'%(__version__).center(71)+'\n'+'\n*Developed by %s'%(__author__).center(71)+'\n\n*Email %s'%(__email__).center(71))
	print('\n'+74*'*'+'\n*'+'These are requirements for running franc tool'.center(71)+'*\n'+74*'*')
	print('\n\nPlatform: Python version >= 2.7.x\n\nOS      : Linux')
	print('\nFor chromopainter: Please make sure "zlib" is installed\n(e.g. sudo apt-get install zlib1g-dev)')
	print('\n\nFor loter        : Please ensure "scikit-allel"" is installed\n')
	print('\nFor supportmix        : Please ensure "libpng12.so.0" is installed\n')
	print('\nIf "-m server", Please, load the following modules\npython >= 2.7, perl, R, and png, by adding them in the "franc.sh\n')
	print('\nOUTPUT FORMAT: franc provides a certain flexibility in the output format.\nwinpop is only converts to lait, while ,\n****pcadmix, supportmix, chromopainter, lampld, elai**** converts to \n****rfmix, loter, winpop, & lait**** in addition to their standard format. \nExcept output in their format, its impossible to convert to \n****pcadmix < supportmix < chromopainter < lampld < elai and multimix****.\n***loter & rfmix**** are only convertible to \n****loter/rfmix or rfmix/loter, winpop & lait**** \nPlease ensure your output format is possible when running a particular tool).\n')
	print(74*'*')
	print('\nEnter 1 to continue and 2 to exit')
	if sys.version_info.major < 3: raw_pass = raw_input
	else: raw_pass = input
	while True:
		a = raw_pass('> ')
		try: 
			a = int(a)
			if a in [1,2]: break
			print('Please enter 1 to continue and 2 to exit')
		except: print('Please enter 1 to continue and 2 to exit')
	if a==2:
		print("\nExiting FRANC, please install packages required \nand try again ...\n")
		sys.exit(2)

	###The parameters
	parameters = readuserinput(args.parameter)
	parameters['admix'] = args.admix; parameters['infolder'] = args.infolder
	parameters['tool'] = [t for t in args.tool.split(':')]
	parameters['mode'] = args.mode
	parameters['outformat'] = args.outformat
	##capturing cwd
	wd = os.getcwd()
	parameters['wrkdir'] = '/'.join(wd.split('/')[:-1])

	if not os.path.exists(args.infolder + '/'+ args.admix):
		os.system('mkdir '+args.infolder+'/'+args.admix)
	for t in args.tool.split(':'):
		if not os.path.exists(args.infolder + '/' + args.admix + '/' + t):
			os.system('mkdir ' + args.infolder + '/'+ args.admix + '/' + t)
		for ch in parameters['CHR']:
			if not os.path.isfile(args.infolder+'/'+args.admix+'.'+str(ch)+'.bed'): 
				print("\nProvide the genotype file.. "+ args.admix +'.'+str(ch)+'.bed')
				sys.exit(2)
			elif not os.path.isfile(args.infolder+'/'+args.admix+'.'+str(ch)+'.bim'):
				print("\nProvide the snps file i.e, "+ args.admix+'.'+str(ch)+'.bim')
				sys.exit(2)
			elif not os.path.isfile(args.infolder+'/'+args.admix+'.'+str(ch)+'.fam'):
				print("\n Please provide all the necessary files i.e, ... your admixed population family file")
				sys.exit(2)
			elif not os.path.exists(args.infolder + '/' + args.admix + '/' + t + '/'+str(ch)):
				os.system('mkdir '+args.infolder+ '/' +args.admix+'/'+t+'/'+str(ch))

	parameters['outfile'] = '.'.join([args.admix, args.tool]) if not args.outfile else args.outfile
	
	# Final check goes here
	if args.mode.lower()=='local':
		RunTool = franc(parameters)
		RunTool.run()
	elif args.mode.lower()=='server':
		fs = open(parameters['wrkdir']+"/franc_interface/ServerInput.txt",'w')
		for k in parameters: 
			fs.write('%s=%s\n'%(k, parameters[k]))
		fs.close()
		fr = open(parameters['wrkdir']+"/franc_interface/francserver.txt").read()%(parameters['wrkdir']+"/franc_interface", parameters['wrkdir']+"/franc_interface/ServerInput.txt")
		fw = open(parameters['wrkdir']+"/franc_interface/francserver.py",'w')
		fw.write(fr)
		fw.close()
		fr = open(parameters['wrkdir']+"/franc_interface/franc.sh").read()%(parameters['wrkdir']+"/franc_interface/",)
		fw = open(parameters['wrkdir']+'/franc_interface/FRANC.sh', 'w')
		fw.write(fr)
		fw.close()
		print('qsub '+parameters['wrkdir']+'/franc_interface/FRANC.sh')
		os.system('qsub '+parameters['wrkdir']+'/franc_interface/FRANC.sh')
	else:
		pass


if __name__=='__main__':
	main()

