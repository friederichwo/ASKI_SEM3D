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
import subprocess
from os import path as os_path
from os import system as os_system


class sgejob(object):
    """
    class holds properties of an SGE job and has functions to run jobs
    """
    def __init__(self, **kwargs):
        """
        Keywords arguments:
            queue:         queue (either low or normal), (def = low)
            hostname:      name of host where job is executed (def = minos15)
            memory:        memory per process in GB (def = 2G)
            ngpu:          number of GPUs to be used (def = 0)
            free_scratch:  required free space on scratch (def = 10G)
            sge_options:   SGE specific options (string)
            command:       Command to be executed by SGE (string) (def = mpirun -n nproc)
            cmd_options:   Options of command (string)
            qsubpath:      Path to qsub command (def = /usr/bin)
        """
        ldlibpath = "/opt/nvidia/hpc_sdk/Linux_x86_64/latest/cuda/lib64"
        mca_par_1 = "--mca fs_ufs_lock_algorithm 1 "          # makes parallel MPI IO fast
        mca_par_2 = "--mca btl openib,self "                  # enables Infiniband (only required if MPI crosses machines)
        mca_par_3 = "--mca btl_openib_allow_ib 1 "            # allows use of IB host channel adapter

        self.queue = 'low'
        self.hostname = 'minos15'
        self.memory = '2G'
        self.gpu = 0
        self.free_scratch = '10G'
        self.sge_options = "-cwd -b yes -pe mpi {nproc:4d} -binding pe linear -o {logfile:s} -sync yes -v OMP_NUM_THREADS=1 -v LD_LIBRARY_PATH=" + ldlibpath
#        self.sge_options = "-cwd -b yes -pe mpi {nproc:4d} -o {logfile:s} -sync yes -v LD_LIBRARY_PATH=" + ldlibpath
        self.command =  "mpirun -n {nproc:4d}"
        self.cmd_options = mca_par_1 + mca_par_2 + mca_par_3
        self.qsubpath = '/usr/bin/'

        if 'queue' in kwargs: self.queue = kwargs['queue']
        if 'hostname' in kwargs: self.hostname = kwargs['hostname']
        if 'memory' in kwargs: self.memory = kwargs['memory']
        if 'ngpu' in kwargs: self.ngpu = kwargs['ngpu']
        if 'free_scratch' in kwargs: self.free_scratch = kwargs['free_scratch']
        if 'sge_options' in kwargs: self.sge_options = kwargs['sge_options']
        if 'command' in kwargs: self.command = kwargs['command']
        if 'cmd_options' in kwargs: self.cmd_options = kwargs['cmd_options']
        if 'qsubpath' in kwargs: self.qsubpath = kwargs['qsubpath']


    def submit(self, nproc, execpath, executable, jobargs, logfile):
        """
        create lists from job properties and pass to subprocess
        """
        command = self.command.format(nproc=nproc).split()
        resources = ["-l", self.queue, "-l", "hostname=" + self.hostname, "-l", "memory=" + self.memory,
                     "-l", "h_vmem=" + self.memory, "-l", "scratch_free=" + self.free_scratch,
                     "-l", "gpu=" + str(self.ngpu)]
        sge_options = self.sge_options.format(nproc=nproc, logfile=logfile).split()
        cmd_options = self.cmd_options.split()

        #  remove logfile if it exists
        if os_path.exists(logfile):
            os_system('rm -f ' + logfile)
            print('rm -f ' + logfile + "\n")

        print([self.qsubpath + "qsub"] + resources + sge_options + command +
                              cmd_options + [execpath+executable] + jobargs.split())
        cp = subprocess.run([self.qsubpath + "qsub"] + resources + sge_options + command +
                              cmd_options + [execpath+executable] + jobargs.split())
        cp.check_returncode()


    def setMemory(self, memory):
        self.memory = memory


    def setNgpu(self, ngpu):
        self.ngpu = ngpu
