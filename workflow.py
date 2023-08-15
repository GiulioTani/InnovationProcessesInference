#!/home/giulio/miniconda3/bin/python3

# @file workflow.py
# @author Giulio Tani (giulio.tani@uniroma1.it)
# @brief Complete workflow for the attribution of books to authors.
# @version 0.1
# @date 2020-07-05
#
# @copyright Copyright (c) 2020
#
# @bug if the suffix contains the strings "seq" or "wnt" the subdirectory structure may fail.

import baa.elaboration as el
import baa.numcomp as nc
import sys
import os
from getopt import getopt
import signal
import shutil
import logging
logging.basicConfig(filename='../common_workflow.log', filemode='w', level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


def signal_handler(sig, frame):
    if sig == signal.SIGHUP:
        print("Ignoring SIGHUP.")
    if sig == signal.SIGINT:
        print('Received SIGINT.')
        sys.exit(0)


logging.info(" ".join(sys.argv))

logging.info("Configure signal handlers")
signal.signal(signal.SIGHUP, signal_handler)
signal.signal(signal.SIGINT, signal_handler)

# @var read
#  Options read from command line.
#
#  Possible options are:
#  -# 'z' to use LZ77 dictionaries;
#  -# 'n' size of the overlapping space free n-grams (incompatible with 'z');
#  -# 'f' [fragment length] size of the book fragments;
#  -# 'w' [window size] size of the window for LZ77;
#  -# 's' suffix to folder name;
#  -# 'c' config file if not named 'config.ini';
#  -# 'l' Size of the slices used to split the corpus (default 1: one book at time);
#  -# 'S' Always create new Slicing:
#  -# 'R' Always Retrieve information from other folders;
#  -# 'O' Always overwrite folder (cancel any R option).
#  -# 'o' Always keep folder (cancel any O option, protects R).
#  -# 'K' Keep seq and wnt folders after computation.
#  -# 'd' Directory where of the corpus database (e.g. the output folder chosen with data preprocessing).
#  -# 'b' Directory where of the base corpus if using processed versions.
#  -# 'D' Delta value for the correction of the base probability. If not passed defaults to 1.
#
# @sa getopt

# @var extra
#  Extra options from command line.
logging.info("Read command line arguments")
read, extra = getopt(sys.argv[1:], "SMPROoKzUf:F:w:s:d:b:B:l:n:c:D:")

# @var param
#  Dictionary containing the parameters for the computation.
logging.info("Calling init function")
param = el.init(read)
logging.info("Initialize numcomp")
nc.init(param['c'])

print("Extracting fragments. . .", end="\r")
# (inputDir, outputDir, baseFragLength, LZ77, ngrams=0, window=-1, useEncoding="", baseCorpus="")
logging.info("Extracting fragments")
el.fragments(param['d'], param['dataDir'], param['f'], param['z'],
             param['n'], param['w'], param['b'], param["U"], param["B"])

logging.info("Preparing support files")
print("Preparing support files. . .", end="\r")
# computes parameters per author and produces the list of fragments for the analisys
el.support_files(param['dataDir'], param['M'])

logging.info("Computing probabilities")
print("Computing probabilities. . .", end="\r")
nc.calprob(param['dataDir'], os.path.join(param['dataDir'], "probFile.prob"),
           slicesize=param['l'], ngram=param['n'], authsize=param['F'])

logging.info("Computing results")
el.results(param['dataDir'], param['D'], margOut=True, allowPartial=param['P'])
if not param['K']:
    print("Cleaning . . .")
    if os.path.isfile(os.path.join(param['dataDir'], "probFile.prob")):
        logging.info("Removing probFile")
        os.unlink(os.path.join(param['dataDir'], "probFile.prob"))
    if os.path.isdir(os.path.join(param['dataDir'], "wnt")):
        logging.info("Cleaning wnt directory")
        shutil.rmtree(os.path.join(param['dataDir'], "wnt"))
logging.info("Done")
