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
#
import subprocess
import logging

logger = logging.getLogger(__name__)


class Basejob(object):
    """
    Base class for handling binary SGE jobs.
    This class holds common attributes  of an SGE job such as
    basic resources and the path to the qsub command.
    Further attributes such as the command, the command options
    are provided by derived classes designed for specific job types
    like an MPI job, a script job or a GPU job.
    The job instance can be used to submit several related jobs
    with different executables whose names and arguments and log files
    are specified with the submit()-function of this class.
    """
    def __init__(self, **kwargs):
        """
            queue:         queue (either low or normal)
            hostname:      name of host where job is executed
            memory:        memory per process in GB
            free_scratch:  required free space on scratch
            qsubpath:      Path to qsub command
        """
        #  defaults
        queue = 'low'; hostname = 'minos15'; memory = '2G'; free_scratch = '10G'
        self.qsubpath = '/usr/bin/'
        self.sge_options =  "-b yes -cwd -sync yes "

        #  set variables from keyword arguments
        if 'queue' in kwargs: queue = kwargs['queue']
        if 'hostname' in kwargs: hostname = kwargs['hostname']
        if 'memory' in kwargs: memory = kwargs['memory']
        if 'free_scratch' in kwargs: free_scratch = kwargs['free_scratch']
        if 'qsubpath' in kwargs: self.qsubpath = kwargs['qsubpath']

        # members set by derived classes
        self.command = ''
        self.cmd_options = ''

        # resources constructed from variables (may be extended by derived classes)
        self.resources = "-l {:s} -l hostname={:s} -l memory={:s} -l h_vmem={:s} -l scratch_free={:s} ".\
            format(queue, hostname, memory, memory, free_scratch)


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
            :key memory:       Amount of memory requested for job per process (default = '2G')
        """
        executable = execpath + execfile
        sge_options = self.sge_options + "-o {:s}".format(logfile)

        resources = self.resources
        if 'memory' in kwargs: 
            resources = self.updateMemory(kwargs['memory'])

        # log start of job
        logger.info(execfile + " submitted")

        # print submission command to joblog
        logger.info("%s qsub %s %s %s %s %s %s",
                    self.qsubpath, resources, sge_options, self.command, self.cmd_options,
                    executable, jobargs)

        # do submit to SGE (arguments must be converted to lists and concatenated)
        cp = subprocess.run((self.qsubpath + "qsub").split() + resources.split() + sge_options.split() +
                            self.command.split() + self.cmd_options.split() + executable.split() + jobargs.split(),
                            capture_output=True, encoding="utf-8")

        # write error output of subprocess to joblog and log end time
        logger.info(execfile + " ended")
        logger.info("Return code: %d", cp.returncode)
        logger.error(cp.stderr)

        return cp


    def updateMemory(self, memory):
        """
        replace the values of memory and h_vmem in resources by new ones
        """
        k1 = self.resources.index('memory')
        k2 = self.resources.index('-l scratch')
        tail = self.resources[k2-1:]
        front = self.resources[0:k1]
        new = "memory=" + memory + " -l h_vmem=" + memory
        return front + new + tail
        








