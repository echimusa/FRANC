# FRANC
FRANC is a unified framework for multi-way local ancestry inference, FRANC, integrating eight existing state-of-the-art local ancestry deconvolution tools. FRANC is an adaptable, expandable and portable tool that manipulates tool-specific inputs, deconvolutes ancestry and standardizes tool-specific results. To facilitate both medical and population genetics studies, FRANC requires convenient and easy to manipulate input files and allows users to choose output formats to ease their use in further potential local ancestry deconvolution applications.

# Cite the below paper
Geza, E., Mulder, N.J., Chimusa, E.R. and Mazandu, G.K., 2020. FRANC: a unified framework for multi-way local ancestry deconvolution with high density SNP data. Briefings in Bioinformatics, 21(5), pp.1837-1845.

Tools included::
	WINPOP
	LAMPLD
	SUPPORTMIX
	RFMIX
	PCADMIX
	ELAI
	CHROMOPAINTER
	LOTER

Name: FRANC

Version: 19.1

Summary: A Python interface for inferring multi-way local 
         ancestry estimates with high density SNP data. It 
         integrates eight existing state-of-the-art ancestry 
         de convolution tools listed above. We refer the 
         interested reader to the PDF file, Appendix 1, for 
         a brief description of each tool, and Appendix 2 
         for a short explanation on the output standardiza-
         tion features implemented in FRANC.

Home-page: http://web.cbio.uct.ac.za/ITGOM/franc

Author: Gaston K Mazandu, Ephifania Geza, 
        Nicola J. Mulder, Emile R. Chimusa
        
Maintainer-email: ecimusa@gmail.com and kuzamunu@aims.ac.za 

License: Copyright (c) 2015 Building Genomics-Bioinformatics 
                            Tool Project@Authors

Users can find essential information about obtaining FRANC 
from its Home-page provided above. It is freely downloadable 
under GNU General Public License (GPL), pre-compiled for 
Linux version and protected by copyright laws. Users are free 
to copy, modify, merge, publish, distribute and display 
information contained in the package, provided that it is done 
with appropriate citation of the package and by including the 
permission notice in all copies or substantial portions of 
the module contained in this package.

Permission is hereby granted, free of charge, to any person 
obtaining a copy of this package and associated documentation 
files to deal in the package without restriction, including 
without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense and to permit persons to whom 
the Software is furnished to do so, subject to the following 
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

FRANC IS DISTRIBUTED AND PROVIDED "AS IS" IN THE HOPE THAT 
IT WILL BE USEFUL, BUT WITHOUT ANY WARRANTY OF ANY KIND, EXPRESS 
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE PACKAGE OR THE USE OR OTHER 
DEALINGS IN THE PACKAGE. 
See <https://www.gnu.org/licenses/gpl-3.0.en.html>.

Description:
	==============================
	The FRANC Python Interface
	==============================
	See the PDF documentation for more information on the use of
	the tool and different multi-way local ancestry deconvolution 
	tools with high density SNP data integrated in FRANC.

	FRANC is a unified framework integrating tools for estimating
	      multi-way local ancestry given high density SNP data of
	      admixed individuals. It contains two Python modules and 
	      sets of required files available in a specific folder 
	      for download. The two modules are described below and 
	      for more details about these and other files, please 
	      refer to the PDF documentation.
	      
         +========================================================+
            Files         Descriptions
         +--------------------------------------------------------+
          franc.py       The main Python code, which is the Python 
                         interface running different local ancestry 
                         deconvolution tools included.
         ----------------------------------------------------------
          convertall.py  The Python file containing functions that 
                         enable the standardization of a given tool 
                         outputs to different output format accor-
                         ding to user specifications.
         +=========================================================+

	Administration
	--------------

	FRANC v19.1 requires Linux operating system and Python(>= 2.7.x). 
	All other libraries and platforms required to deconvolve local 
	ancestry using the FRANC interface greatly depend on the choice 
	of the tool being executed. These libraries need to be installed 
	prior to the use of FRANC.
	
	* CHROMOPAINTER requires the system to have zlib, which is 
	installed using the following command:
	
	    sudo apt-get install zlib1g-dev
	
	* LOTER requires Scikit allel package,
	* SUPPORTMIX requires libpng12.so.0,
	* CHROMOPAINTER and LAMPLD requires Perl,
	* RFMIX requires R.

	Build status
	------------

	The main website for the FRANC interface is 
	http://web.cbio.uct.ac.za/ITGOM/franc where users can find 
	essential information about FRANC. It is freely downloadable
	under GNU General Public License 
	(GPL: https://www.gnu.org/licenses/gpl-3.0.en.html)
	pre-compiled for Linux version and protected by copyright 
	laws. Users are free to copy, modify, merge, publish, 
	distribute and display information contained in this package, 
	provided that it is done with appropriate citation of the 
	package and by including the permission notice in all copies
	or substantial portions of the module contained in this 
	package.

	The whole package is relatively large (around 110 Mb) and 
	contains two Python modules and sets of required files 
	available 	in a specific folder for download.
	It is currently maintained by one member of the core-development 
	team, Gaston K. Mazandu <kuzamunu@aims.ac.za>, who regularly 
	updates the information available in this package and makes every 
	effort to ensure the quality of this information.

	FRANC usage
	-----------
	
	To use FRANC, the user needs to download the `tar.gz' file 
	from the FRANC Home-page and extract all files as follows: 

        tar -xzf franc.tar.gz
	or
	git clone https://github.com/echimusa/FRANC.git
	
	It is worth mentioning that FRANC is a portable python package.
	So, the user is required to move in the franc_interface folder, 
	containing the franc.py python module to run different multi-
	way ancestry deconvolution tools. Thus, after uncompressing 
	franc.tar.gz, please move to the franc_interface directory, 
	using the following terminal command::
	     
	     	cd FRANC/Interface. 
	
	
	The FRANC interface is free to use under GNU General Public 
    License. You are free to copy, distribute and display information
	contained herein, provided that it is done with appropriate cita-
	tion of the tool. Thus, by using FRANC, it is assumed that you have
	read and accepted the agreement provided and that you agreed to be 
	bound to all terms and conditions of this agreement. Please, use 
	the following command line to see the package licence::

			python setup.py --licence

	 
	To get help on how to run FRANC using the following command::

		python franc.py -h

	The above command should produce the following output:
	
	usage: FRANC.py [-h] -p FILE [-d FILE] -t FILE -a FILE [-o FILE] [-m MODE] -f FILE 

	with different tags explained on your screen.

	The two commands above should be executed inside the franc_interface 
	folder.
	
	
	Running FRANC
	-------------
	As highlighted by the help option, FRANC is run using the following 
	one line command::

		python franc.py -p FrancPar.txt -d infolder -t tool -a admixed_prefix_name -o outfile -m local -f output_format

	When running more than one tool, these tools should be 
	separated by semicolon and provided after the tag t in 
	the command. For example, to run LAMPLD, WINPOP and RFMIX, 
	then -t lampld:winpop:rfmix. When FRANC	is run, an 
	interactive interface appears orienting the user of 
	systems, platforms or libraries requirements of the 
	framework. If the system satisfies all the requirements, 
	the user proceeds by entering 1, or alternatively the user 
	can enter 2 to exit when requirements are not satisfied.

	It is worth mentioning that FRANC can also be run from 
	server using the tag -m server. In this case, the user 
	should fill in the server details in the example.sh file. 
	Details include queue name, resources requested which 
	include loading modules required to run FRANC. 
	
	
	(a) Setting FRANC parameter file (FrancPar.txt)
	-----------------------------------------------
	Additional required parameters are provided in a parameter 
	file, FrancPar.txt, which is inside the franc_interface 
	directory. This parameter file also provides the structure 
	of any additional FRANC required parameters. Thus, it is 
	recommended that the user maintains the order of this file, 
	and edit it according to his/her specifications. The user 
	can also create a new parameter file with the same structure, 
	in which case, it is required to provide the path to this 
	parameter file after the tag 'p' or -p in the command line 
	running FRANC. Note that each tool requires this file to run 
	and more details on the contents of FrancPar.txt are given 
	in the PDF documentation and, for a given tool, not all 
	options are required.
	
	(b) FRANC file management
	-------------------------
	At the high level, FRANC system is composed of two main 
	folders as shown in Figure 1 in the PDF documentation: 
	franc_interface, franc_util and test folders. The test 
	folder contains an illustrative dataset for testing the 
	FRANC interface, as shown in in (d).
	
	- The franc_interface folder includes two main Python 
	  modules and one text file, namely franc.py and 
	  convertall.py described above, and FrancPar.txt, 
	  explained in (a). It also contains config file and bash 
	  and other python modules for server option.
	  
	- The franc_util folder contains two sub-folders: soft and
	  configs. soft contains source codes of all tools included 
	  in FRANC. configs contains following files required to 
	  execute tools included in FRANC: configPedGeno.txt, 
	  configsvm.txt, par_winpop.txt, configRfmix.txt,  
	  configGenoPed.txt and configGenoPed.py. 
	  The two files, configGenoPed.txt and configGenoPed.py 
	  are needed for converting the EIGENSTRAT/OXFORD to 
	  PLINK ped/map formats. 
	
	(c) FRANC input files
	---------------------
 	These are user input files specific to the population 
 	under consideration, which are in a user specified 
 	folder, and here referred to it as infolder. There are 
 	three input files: population data, genetic map and 
 	reference genetic maps. The user is required to provide 
 	the path to this infolder or alternatively copy or moved 
 	them into the franc_interface folder, which is the 
 	folder from which the user runs FRANC. In this case the 
 	infolder and the tag 'd' may not be used as 
 	franc_interface is considered to be the default folder 
 	that may contain these input files.
 	
 	(d) FRANC output files
 	----------------------
 	The FRANC output folder, named admix_prefix_name, is 
 	created within infolder after running franc.py to contain
 	the result outputs of executed tools (see Figure 1, 2 and
 	3 in the PDF documentation for more details about the 
 	structure of this folder). These outputs depend greatly 
 	on the user specifications. Each tool produces its output 
 	file format by default. In addition to their output format, 
 	PCADMIX,  	SUPPORTMIX, CHROMOPAINTER, LAMPLD, MULTIMIX for
 	phased and resolved and ELAI (phased option), can produce 
 	results in the format of RFMIX, LOTER, WINPOP and LAIT 
 	(refer to Appendix-2 in the PDF documentation for more 
 	details).

	It should be noted that, running any other tool different 
	from itself can not yield results in PCADMIX, SUPPORTMIX, 
	CHROMOPAINTER, LAMPLD and ELAI formats. Thus, RFMIX and 
	LOTER are convertible to LOTER/RFMIX, WINPOP and LAIT. 
	Although this works for two- and three-way admixtures only, 
	all output results from any tool within FRANC can be 
	standardized to the LAIT format. 
	
	Illustrating FRANC usage
	------------------------
	To illustrate how should a specific file be, please go to  
	http://web.cbio.uct.ac.za/ITGOM/franc/tests and download  
	or alternatively, go to the test_data sub-folder in the local 
	franc forder to see the illustrative data set  to  that can 
	be used to try and run FRANC. Refer to the PDF documentation
	for more details about this dataset.
	
	FRANC is run with the following command:

		python franc.py -p FrancPar.txt -d ../test_data/ -a SIM2 -t chromopainter -f rfmix -o SIM2 -m local
	
	The output results can be viewed under the folder test_data
	in the sub-folder SIM2.

#Version history
---------------

- 19.1: Initial FRANC release in January 2019.

# Maintainer
----------

Emile R. Chimusa
Email: echimusa@gmail.com 

#Contributors
------------

Ephifania Geza, Nicola J Mulder, Emile R Chimusa, 
Gaston K Mazandu

        
Classifier: License :: GPL (>= 3)
Classifier: Operating System :: OS Independent, 
            but tested only on Linux (Ubuntu)
Classifier: Programming Language :: Python (>=2.7.x), but 
            not tested on the version less than 2.7
Classifier: Topic :: Software Development :: Libraries
