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
#   Class handling an SGE parallel job
#
from askipy.mpijob import Mpijob


class Gpujob(Mpijob):
    """
    Derived class for handling GPU jobs under MPI.
    Here, GPU specific resources and options are added to the MPI ones.
    Allows to run different GPU jobs with different number of processes
    and GPUs using the submit()-function.
    """
    def __init__(self, **kwargs):
        """
        Keyword arguments: see base class sgejobBase
        """
        # initiate base class variables
        super().__init__(**kwargs)

        # add CUDA library path to options
        ldlibpath = "/opt/nvidia/hpc_sdk/Linux_x86_64/latest/cuda/lib64"
        self.sge_options += "-v LD_LIBRARY_PATH={:s} ".format(ldlibpath)


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
            :key ngpu:         Number of GPUs (default = 0)
            :key nproc:        Number of processes (default = 1)
            :key memory:       Amount of memory requested for job per process (default = '2G'), passed to Basejob
        """
        ngpu = 0
        if 'ngpu' in kwargs: ngpu = kwargs['ngpu']

        # store original resources
        resources = self.resources

        # add GPUs to resources
        self.resources += "-l gpu={:d} ".format(ngpu)

        #  call submit function of parent class
        cp = super().submit(execpath, execfile, jobargs, logfile, **kwargs)

        #  restore original resources
        self.resources = resources

        return cp




