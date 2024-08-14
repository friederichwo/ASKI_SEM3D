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
from askipy.basejob import Basejob


class Mpijob(Basejob):
    """
    Derived class for handling MPI jobs.
    Here, MPI specific resources and options are added to the basic ones.
    The specific command and command options for MPI jobs are defined here.
    Allows to run different MPI jobs with different number of processes
    using the submit()-function.
    """
    def __init__(self, **kwargs):
        """
        Keyword arguments: see base class
        """
        # initiate base class variables
        super().__init__(**kwargs)

        #  options specific for MPI jobs
        self.sge_options += "-binding pe linear -v OMP_NUM_THREADS=1 "

        #  mpirun command and options
        mca_par_1 = "--mca fs_ufs_lock_algorithm 1 "           # makes parallel MPI IO fast
        mca_par_2 = "--mca btl openib,self "                   # enables Infiniband (only required if MPI crosses machines)
        mca_par_3 = "--mca btl_openib_allow_ib 1 "             # allows use of IB host channel adapter
        self.command = "mpirun -n "
        self.cmd_options += mca_par_1+mca_par_2+mca_par_3


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
            :key nproc:        Number of processes (default = 1)
            :key memory:       Amount of memory requested for job per process (default = '2G'), passed to Basejob
        """
        nproc = 1
        if 'nproc' in kwargs: nproc = kwargs['nproc']

        #  store original values
        command = self.command
        sge_options = self.sge_options

        # adjust values
        self.command += "{:d}".format(nproc)
        self.sge_options += "-pe mpi {:d} ".format(nproc)

        #  call base class submit function
        cp = super().submit(execpath, execfile, jobargs, logfile, **kwargs)

        # restore original values
        self.command = command
        self.sge_options = sge_options

        return cp




