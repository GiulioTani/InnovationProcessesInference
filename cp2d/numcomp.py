#!/usr/bin/env python3

# CP2D -- Constrained Probability Poisson-Dirichlet
# Copyright (C) 2023  Giulio Tani Raffaelli

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# @file numcomp.py
# @author Giulio Tani Raffaelli (tani@cs.cas.cz)
# @brief Functions launching the numeric computation to extract LZ77 dictionaries and the conditional probabilities.
# @version 1.0
# @date 2023-06-27
# 
# @copyright Copyright (c) 2023
#

# @package numcomp
#  Functions launching the numeric computation to extract LZ77 dictionaries and the conditional probabilities.
#
#  Although is possible to use C code directly from Python is quite annoying so I prefer to call external programs.
#  All such programs are called and managed from di module.
#
#  The main aim of this mondule is to manage the core computation of the conditional probabilities of the fragments over the authors.
#  This is done by fragprob that is quite resource intensive: this module includes also a workload manager that in case other users are on the same machine reduces the resources used.

import sys
import os
from datetime import timedelta
from subprocess import Popen, PIPE
import multiprocessing as mp
import psutil
import signal
import time
from select import select
import configparser
import warnings
from multiprocessing import Pool
import logging
from glob import glob

logger = logging.getLogger(__name__)
verboseDebugging = False
# @var numCpus
#  Total number of cores available. Hard coded default value overwritten from initialization file.
#

# @var buffer
#  Number of always free cores to avoid frequent changes in the number of running threads. Hard coded default value overwritten from initialization file.

# @var limit
#  Maximum number of cores left to other users. Hard coded default value overwritten from initialization file.
#

# @var timeout
#  Interval between stdin readouts. Hard coded default value overwritten from initialization file.
#

# @var cpuThreshold
#  Minimum cpu percentage to consider a process. Hard coded default value overwritten from initialization file.
#
numCpus = 40
buffer = 3
limit = 27  # 7#27
timeout = 0.5
cpuThreshold = 10
el_path = os.path.dirname(__file__)


def init(iniFile):
    """
    Reads config file.

    Reads the config file to get the values of the parameters controlling the job manager and the numeric computations.

    In case of badly written config file tries to deduce some of the parameters from the system or the values of others, if fails falls back to the default value.

    Also performs some checks to avoid strange behaviours when the number of cores to leave free is higher than the total number of cores and similar things.
    In case of incompatible balues issues a warning and falls back to the less resource consuming valid configuration.

    Args:
        iniFile:    path to the config file.
    """
    global numCpus, buffer, limit, timeout, cpuThreshold
    if iniFile:
        if not os.path.isfile(iniFile):
            raise FileNotFoundError(f"Config file {iniFile} not found.")
        config = configparser.ConfigParser()
        config.read(iniFile)
        if "numcomp" in config:
            if config["numcomp"].getint("numCpus") == None:
                logger.warning(f"Missing 'numCpus' value in config file")
                print(
                    "Missing 'numCpus' value in config file. Trying to deduce from system.", file=sys.stderr)
                try:
                    numCpus = mp.cpu_count()  # relying on multiprocessing
                except:
                    logger.warning(f"Failed to get number of cpus")
                    # we shouldn't get here
                    print(
                        "Failed to get number of cpus. Going on with default value, fingers crossed.", file=sys.stderr)
            else:
                numCpus = config["numcomp"].getint('numCpus')

            if config["numcomp"].getint("buffer") == None:
                print(
                    "Missing 'buffer' value in config file. Going on with default value.", file=sys.stderr)
                if "elaboration" in config and "workers" in config["elaboration"] and config["elaboration"].getint("workers") < numCpus:
                    # this is complementary to the analogous deduction in elaboration.py
                    buffer = numCpus - config["elaboration"].getint("workers")
                else:
                    logger.warning(f"Failed to deduce 'buffer' size")
                    print(
                        "Failed to deduce 'buffer' size. Going on with default value.", file=sys.stderr)
            else:
                buffer = config["numcomp"].getint('buffer')

            if config["numcomp"].getint("limit") == None:
                logger.warning(f"Missing 'limit' value in config file")
                print(
                    "Missing 'limit' value in config file. Going on with default value.", file=sys.stderr)
            else:
                limit = config["numcomp"].getint('limit')

            if config["numcomp"].getfloat("timeout") == None:
                logger.warning(f"Missing 'timeout' value in config file")
                print(
                    "Missing 'timeout' value in config file. Going on with default value.", file=sys.stderr)
            else:
                timeout = config["numcomp"].getfloat('timeout')

            if config["numcomp"].getfloat("cpuThreshold") == None:
                logger.warning(f"Missing 'cpuThreshold' value in config file")
                print(
                    "Missing 'cpuThreshold' value in config file. Going on with default value.", file=sys.stderr)
            else:
                cpuThreshold = config["numcomp"].getfloat('cpuThreshold')

            # checks for invalid numbers
            if limit > numCpus:
                logger.warning(f"Limit higher than numCpus, resetting")
                warnings.warn("Limit higher than numCpus, resetting.")
                limit = numCpus
            if buffer > numCpus:
                logger.warning(f"Buffer higher than numCpus")
                warnings.warn(
                    "Buffer higher than numCpus. I assume you want the work done, changing to numCpus-1.")
                buffer = numCpus-1
            if buffer > limit:
                logger.warning(
                    f"Buffer higher than limit, job manager will have no effect")
                warnings.warn(
                    "Buffer higher than limit, job manager will have no effect.")
        else:
            logger.warning(f"Missing 'numcomp' section in config file")
            print(
                "Missing 'numcomp' section in config file. Going on with default values.", file=sys.stderr)
    else:
        logger.warning(f"No config file provided")
        print("No config file provided. Going on with default values.", file=sys.stderr)


def process_controller():
    """
    Self-implemented workload manager. Checks if other users are using the cluster.

    If up to 10 processes are detected mid intensity usage is assumed and the same amount of cores are freed. If more than 10 processes are fund restricts the anlisys to 18 cores.
    If the other users' processes terminate the cores are populated again.
    """
    nomi = ["betta", "claudioc", "londei", "napoletano", "servedio",
            "ubaldi", "campanelli", "gravinno", "monechi", "schueller", "guest"]
    now = 0
    # reads once to initialize cpu_percent
    process = [p for p in psutil.process_iter(['username', 'cpu_percent'])]
    try:
        logger.info(f"Starting process controller loop")
        my_proc = []
        while True:
            time.sleep(2)
            count = 0
            target = 0
            try:
                # looks every process
                for process in psutil.process_iter(['username', 'cpu_percent', 'exe'], "0.0"):
                    # to avoid shrinking from shells
                    if process.info['username'] in nomi and float(process.info['cpu_percent']) > cpuThreshold:
                        count += 1
                    if not my_proc and 'fragprob' in process.info['exe']:
                        my_proc = [psutil.Process(
                            A.id) for A in process.threads() if A.id != process.pid]
                        my_proc = [p for p in my_proc if float(
                            process.info['cpu_percent']) > cpuThreshold]
                        my_proc.sort(key=lambda x: x.pid)
                        logger.info(f"Found {len(my_proc)} threads")
                        old_r = my_proc[:]
                        old_s = []
                if not my_proc:
                    continue
                tmp_buff = numCpus-len(my_proc)
                if count > tmp_buff and limit > tmp_buff:
                    target = (count-tmp_buff if count <
                              limit else limit-tmp_buff)
                else:
                    target = 0
                my_proc = old_r[:]+old_s[:]
                new_s = my_proc[:target]
                new_r = my_proc[target:]
                if target != now:
                    logger.info(f"Paused processes from {now} to {target}")
                    print(f"Paused processes from {now} to {target}.")
                    now = target
                stopping = 0
                for pr in new_s:
                    if pr in old_r or pr.status() == psutil.STATUS_RUNNING:
                        stopping += 1
                        # don't know why I can't effectively catch SIGSTOP, stops everything
                        pr.send_signal(signal.SIGUSR1)
                if stopping:
                    logger.info(f"Sent {stopping} stop signals")
                starting = 0
                for pr in new_r:
                    if pr.status() != psutil.STATUS_RUNNING:
                        pr.resume()
                        if pr in old_s:
                            starting += 1
                if starting:
                    logger.info(
                        f"Sent {starting} start signals (first: {new_r[0].pid})")
                old_s = new_s[:]
                old_r = new_r[:]
            except psutil.NoSuchProcess:
                old_s = [pr for pr in old_s if pr.is_running()]
                old_r = [pr for pr in old_r if pr.is_running()]
                logger.info(
                    f"Now {len(old_s)+len(old_r)} processes left running")
                if len(old_s)+len(old_r) == 0:
                    break
                #print("Subprocess terminated while assessing . . . ¯\_ (ツ)_/¯")
    except KeyboardInterrupt:
        pass
    except Exception as ecc:
        logger.error("Uncaught exception in process controller")
        logger.error(ecc)
        logger.error(type(ecc))
        raise


def calprob(dataDir, resultsDir, outFile="", slicesize=0, ngram=0, fragment=0, authsize=0, P0file: str="", dumpP0: bool=False):
    """
    Organizes the computation of the probabilities.

    Args:
        dataDir:    Directory containing the input data.
        probFile:   Name of the file containing probability.
        outFile:    Name of the file where to store the results.
    """
    if outFile:
        outName = outFile
    else:
        name = os.path.split(dataDir)
        if not name[1]:
            name = os.path.split(dataDir)
        outName = os.path.join(dataDir, name[1])
    if len(glob(os.path.join(resultsDir, outName+"_*.res"))):
        return
    if len(glob(os.path.join(resultsDir, "tmpres_*.res"))) and not len(glob(os.path.join(resultsDir, "tmpres_uo_*.res"))):
        for fname in glob(os.path.join(resultsDir, "tmpres_*.res")):
            nnam=os.path.basename(fname).split("_")[1]
            os.rename(fname, os.path.join(resultsDir, "_".join([outName,nnam])))
        return
    t_start = time.time()
    # \n\tl : approximate number of jobs remaining
    help_msg = "Press <?> + return for info on the execution:\n\tt : Time since the beginning of this stage\n\tr : number of currently Running processes\n\th : this help"
    help_short = "Press [t,r,h] + return for info on the execution."
    print(help_short)
    try:
        logger.info(f"Starting process controller")
        p = [mp.Process(target=process_controller, args=()), ]
        #p[0].start()
        logger.info(f"Starting probability computation")
        p.append(mp.Process(target=prob_compute, args=(
            dataDir, resultsDir, slicesize, ngram, fragment, authsize, P0file, dumpP0)))
        p[1].start()
        logger.info(f"Starting interactive loop")
        while p[1].is_alive():
            rlist, _, _ = select([sys.stdin], [], [], timeout)
            if rlist and sys.stdin.isatty():
                inp = sys.stdin.readline().strip()
                if inp == "t":
                    print(
                        f"{timedelta(seconds=int(time.time()-t_start))} running this step.")
                elif inp == "r":
                    now_running = 0
                    try:
                        process = [pr for pr in psutil.process_iter(
                            ['cpu_percent', ])]
                        time.sleep(0.1)
                        # looks every process
                        for process in psutil.process_iter(['cpu_percent', 'exe'], "0.0"):
                            if "fragprob" in process.info['exe'] and float(process.info['cpu_percent']) > cpuThreshold:
                                now_running += process.num_threads()
                    except:
                        logger.warning("Error retrieving processes info")
                        print("Error retrieving processes info\n",
                              process, process.info)
                    print(f"{now_running-1} processes running.")
                    if now_running <= 0:
                        logger.info(f"Exiting interactive loop")
                        print("Exiting loop.")
                        psutil.Process(p[1].pid).send_signal(signal.SIGINT)
                        break
                elif inp == "h":
                    print(help_msg)
                else:
                    print("Unknown option.", help_msg, sep="\n")
        logger.info(f"Closing computation processes")
        p[1].join()
        #p[0].terminate()
        #p[0].join()
        if p[1].exitcode:
            raise mp.ProcessError(p[1].exitcode)
        for proc in p:
            proc.close()
        logger.info(f"Renaming and moving output file")
        print(os.path.join(resultsDir, "tmpres_ids.res"), outName+"_ids.res")
        for fname in glob(os.path.join(resultsDir, "tmpres_*.res")):
            nnam=os.path.basename(fname).split("_")[1]
            os.rename(fname, os.path.join(resultsDir, "_".join([outName,nnam])))
        print(f"Elapsed: {timedelta(seconds=int(time.time()-t_start))}.")
    except:
        logger.error("Exception has occurred")
        print("Exception has occurred. Don't trust the output file.")
        raise


def prob_compute(dataDir, resultsDir, slicesize, ngram, fragment, authsize, P0file: str="", dumpP0: bool=False):
    """
    Launches the computation of conditional probabilities

    Args:
        dataDir:    Directory containing the input data.
        probFile:   Name of the file containing probability.

    Raises:
        RuntimeError
    """
    args = [os.path.join(el_path, "../bin/fragprob"), "-f", os.path.join(dataDir, "seq"),
            "-o", resultsDir, "-F", str(fragment), "-t", str(numCpus-buffer)]
    if slicesize:
        args.extend(['-s', str(slicesize)])
    else:
        args.extend(['-s', "1"])
    if ngram:
        args.extend(['-n', str(ngram)])
    if authsize:
        args.extend(['-d', str(authsize)])
    if P0file:
        args.extend(['-p', P0file])
    elif dumpP0:
        args.append('-P')
    logger.info(f"Launching C++ code")
    print(" ".join(args))
    sys.stdout.flush()
    process = Popen(args, stderr=PIPE)  # gdb --args
    (_, err) = process.communicate()
    exit_code = process.wait()
    if exit_code:
        logger.error(f"Error in probabilities computation")
        logger.error(err.decode())
        print(f"\n\nError in probabilities computation.\n",
              err.decode(), flush=True)
        raise RuntimeError("Failed computation.")
    else:
        logger.warn(err.decode())


if __name__ == '__main__':
    print("This script is intended only as a support for workflow.py. Use that.", file=sys.stderr)
    exit(1)
