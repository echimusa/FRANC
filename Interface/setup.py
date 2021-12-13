#!/usr/bin/env python

from distutils.core import setup

LICENSE = open("LICENSE").read()

LONG_DESCRIPTION = open("PKG-INFO").read().replace("`_", "`")

setup(name='FRANC',
   version='19.1',
   description='\nFRANC: A framework that allows a user to estimate the ancestry of each chromosomal segment of an admixed individual, \ncommonly known as local ancestry inference or local ancestry deconvolution.\n',
   long_description=LONG_DESCRIPTION,
   author='Gaston K. Mazandu et al.',
   author_email='kuzamunu@aims.ac.za, ephie@aims.ac.za, nicola.mulder@uct.ac.za, emile.chimusa@uct.ac.za',
   maintainer = 'Gaston K. Mazandu',
   maintainer_email = 'gmazandu@gmail.com, kuzamunu@aims.ac.za',
   url='http://web.cbio.uct.ac.za/ITGOM/franc',
   license=LICENSE,
   classifiers= [ "\nDevelopment Status :: 4 - Beta",
                  "License :: OSI Approved :: GNU General Public License",
                  "Operating System :: OS Independent, but tested only on Linux",
                  "Programming Language :: Python :: Not tested on the version less than 2.7.x",
                  "Topic :: Software Development :: Libraries",
                  "Following programming languages need to be installed prior to using FRANC:",
                   "\t::Python (>=2.7.x)\n\t::R\n\t::Perl\n"],
   platforms = 'FRANC is OS Independent, but it was tested only on Linux',
)
