from os import path as os_path
from sys import stdout
from os import getcwd
from time import time as ttime
from time import ctime as tctime


def check_file(filename):
    """Check existence of a file and produce error message in stdout"""
    if not os_path.exists(filename):
        stdout.write("### ERROR : The file '" + filename + "' does not exist \n\n")
        raise Exception("Missing file")


def check_directory(dirname):
    """Check existence of a directory and produce error message in stdout"""
    if not os_path.isdir(dirname):
        stdout.write("### ERROR : The directory '" + dirname + "' does not exist \n\n")
        raise Exception("Missing directory")


def check_keys(nokeys):
    if len(nokeys) > 0:
        stdout.write("### ERROR! the following keywords are required:".join(nokeys) + "\n")
        raise Exception("Missing parameters")


def check_free(pathname):
    """Raise error if a path exists"""
    if os_path.exists(pathname):
        raise Exception("Path " + pathname + " already exists")


def check_in_mainpath(mainpath):
    if mainpath != getcwd() + '/':
        raise Exception("Run script in main inversion directory")


def get_time_string():
    return tctime(ttime())
