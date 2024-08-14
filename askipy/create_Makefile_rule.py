#!/usr/bin/env python3
"""
   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)

   This file is part of ASKI version 1.2.

   ASKI version 1.2 is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   ASKI version 1.2 is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.

   Python script that creates a Makefile rule for an ASKI executable given by
   sys.argv[1], also accounting for library requirements of its dependencies.
"""
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
# ------------------------------------------------------------
#   Functions
# -------------------------------------------------------------


def line_break(line):
    """
    break lines after max_char_line characters (in make line break syntax)
    """
    max_char_line = 100
    words = str.split(line, " ")
    broken_line = ""
    nchar = 0
    for i, w in enumerate(words):
        nchar += len(w)
        if nchar < max_char_line:
            broken_line += w + ' '
        else:
            broken_line += '\\\n\t' + w + ' '
            nchar = len(w)
    return broken_line


# -----------------------------------------------------------
#   START MAIN SCRIPT HERE
# -----------------------------------------------------------
ap = ArgumentParser(description="Generate a Makefile rule for a given target",
                    formatter_class=ArgumentDefaultsHelpFormatter)
ap.add_argument('target', help='Target of Makefile rule')
ap.add_argument('--depfile', help='Name of file with dependencies', default='dependencies')
args = ap.parse_args()
target = args.target
depfile = args.depfile

depobs = []  # list of all objects on which target depends
#
#   read lines of dependencies file
#
with open(depfile, 'r') as fdep:
    dep_lines = fdep.readlines()
#
#   search for target in dependencies (should be first word in line)
#   add depending objects to depobs
#
for line in dep_lines:
    objects = str.split(line)
    if target in objects[0]:
        depobs += objects[1:]
        break
#
#   now add dependencies of each of the depending objects to depobs
#
for obj in depobs:
    for line in dep_lines:
        objects = str.split(line)
        if obj in objects[0]:
            depobs += objects[1:]
#
#   sort and make unique list of depending objects and
#   print (broken) line(s) of target dependencies
#
depobs_unique = sorted(set(depobs))
makerule = line_break(target+': %: %.o '+" ".join(depobs_unique))
print(makerule)
