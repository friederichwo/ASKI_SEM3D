# ----------------------------------------------------------------------------
#   Copyright 2024 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of ASKI version 1.2.
#
#   ASKI version 1.2 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   ASKI version 1.2 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
#---------------------------------------------------------------------------------
#   Class handling an SGE python job
#
from askipy.basejob import Basejob


class Pyjob(Basejob):
    """
    Derived class for handling non-parallel Python jobs.
    Here, Python specific resources and options are added to the basic ones.
    The specific command and command options for Python jobs are defined here.
    Allows to run different Python jobs using the submit()-function.
    """
    def __init__(self, **kwargs):
        """
        Keyword arguments:
            pypath:      Colon separated list of paths to Python code
            condapath:   Path to conda executable
            further arguments: see base class
        """
        pypath = "/home/wolle/work_git/ASKI_SEM3D_CAR/"
        condapath = "/opt/anaconda3/bin/"
        if 'pypath' in kwargs: pypath = kwargs['pypath']
        if 'condapath' in kwargs: condapath = kwargs['condapath']

        # initiate base class variables
        super().__init__(**kwargs)

        #  options specific for MPI jobs
        self.sge_options += "-b yes -v OMP_NUM_THREADS=1 " + "-v PYTHONPATH={:s} ".format(pypath)

        #self.command = condapath + "conda"
        #self.cmd_options += "-n wolle38 -p "


    def submit(self, execpath, execfile, jobargs, logfile, **kwargs):
        """
        Submit function which calls base class submit
        :param execpath: path to executable
        :param execfile: name of executable (without path)
        :param jobargs: arguments to executable program
        :param logfile: file for logging output of program
        :param kwargs: further keyword arguments passed on to Basejob
        :return: a CompletedProcess object
        Keywords:
            :key memory: Amount of memory requested for job per process (default = '2G'), passed to Basejob
        """
        cp = super().submit(execpath, execfile, jobargs, logfile, **kwargs)
        return cp


