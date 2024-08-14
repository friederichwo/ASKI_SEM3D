#!/usr/bin/env python3
""" --------------------------------------------------------------------------
   Copyright 2013 Wolfgang Friederich

   This file is part of Gemini II.

   Gemini II is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   any later version.

   Gemini II is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.

   Python module to create dependencies from use and include
   statements in Fortran programs

   Usage: makeDepFromUseInclude [[dir] [dir]....]
   Search all f90 files in . and given directories and create a list of basename.o: deps.o lines
   to be included into a Makefile
"""
import glob
import re
import os.path
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
# ------------------------------------------------------------------
#    Functions
# ---------------------------------------------------------------


def search_text_file(file, pattern):
    """ Search each line of a text file for a pattern.
    Start searching at the beginning of the line.
    Returns the group(1) match.
    """
    mlist = []
    with open(file, 'r') as f:
        for line in f:
            m = re.match(pattern, line)
            if m is not None:
                mlist.append(m.group(1))
    return mlist


def build_makefile_entry(file, mlist, depext):
    """ Build a makefile entry.
    Form: datei_without_path_and_extension.o: dependencies.
    mlist: list of dependencies
    depext: extension to be appended to the dependencies
    """
    target = str.split(os.path.basename(file), '.')[0] + '.o:'
    dep = target
    for el in set(mlist):
        dep = dep + ' ' + el + depext
    return dep


# --------------------------------------------------------------------
#   Start of main script
# --------------------------------------------------------------------
ap = ArgumentParser(description="Generate a dependency file for use with make",
                    formatter_class=ArgumentDefaultsHelpFormatter)
ap.add_argument('--dirs', help='List of directories to be searched', default='f90')
ap.add_argument('--depfile', help='Name of file to which dependencies are written', default='dependencies')
args = ap.parse_args()

fdep = open(args.depfile, 'w')
for folder in args.dirs.split():

    #   search for use statements in f90 files
    #
    for datei in glob.glob(folder+'/*.f90'):
        matchlist = search_text_file(datei, r'\tuse ([a-zA-Z0-9_]+)')
        if len(matchlist) != 0:
            entry = build_makefile_entry(datei, matchlist, '.o')
            fdep.write(entry+'\n')
        matchlist = search_text_file(datei, r'^ {2,4}use ([a-zA-Z0-9_]+)')
        try:
            matchlist.pop(matchlist.index('hdf5'))
        except ValueError:
            pass
        try:
            matchlist.pop(matchlist.index('mpi'))
        except ValueError:
            pass
        if len(matchlist) != 0:
            entry = build_makefile_entry(datei, matchlist, '.o')
            fdep.write(entry+'\n')

    #   search for include statements in f-files, either tab or four blanks
    #
    for datei in glob.glob(folder+'/*.f'):
        matchlist = search_text_file(datei, r'\tinclude \'([a-zA-Z0-9_]+)')
        if len(matchlist) != 0:
            entry = build_makefile_entry(datei, matchlist, '.h')
            fdep.write(entry+'\n')
        matchlist = search_text_file(datei, r'    include \'([a-zA-Z0-9_]+)')
        if len(matchlist) != 0:
            entry = build_makefile_entry(datei, matchlist, '.h')
            fdep.write(entry+'\n')

    #   search for include statements in .f.m4-files
    #
    for datei in glob.glob(folder+'/*.f.m4'):
        matchlist = search_text_file(datei, r'\tinclude \'([a-zA-Z0-9_]+)')
        if len(matchlist) != 0:
            entry = build_makefile_entry(datei, matchlist, '.h')
            fdep.write(entry+'\n')

fdep.close()
