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

import argparse
import os
from . import support, elaboration


def main ():
    parser = argparse.ArgumentParser(
        prog="CP2D", description="Authorship attribution using the CP2D  approach. For each subcommand help, use <subcommand> -h.", epilog="Giulio Tani Raffaelli, 2023")
    parser.add_argument('--version', action='version', version='%(prog)s 1.1')
    subparsers = parser.add_subparsers(required=True, dest='command')

    run_parser = subparsers.add_parser('attribute', help="Attributes the the documents in a given corpus to their authors.")
    run_parser.add_argument('database', type=str, help="Path to the database folder (e.g. the output folder chosen with data preprocessing subcommand 'prepare').", metavar='corpus path')

    group_run1 = run_parser.add_mutually_exclusive_group()
    group_run1.add_argument('-n', type=int, dest="ngramSize", metavar='N', default=0,
                            help="Size of the overlapping space free n-grams (incompatible with 'z' and 'U')")
    group_run1.add_argument('-U', action='store_true',
                        help="Keep underscores and punctuation when using dictionary method (incompatible with 'z' and 'n')", dest="keepNonAlfabetical")
    group_run1.add_argument('-z', action='store_true',
                        help="Use LZ77 dictionaries (incompatible with 'n' and 'U')", dest="LZ77")
    run_parser.add_argument('-w', type=int, help="Length of the LZ77 window (default: use full document)."
                            "A length<=0 searches on the whole document. Usually long windows perform better. "
                            "The window length is limited by the document size, avoid using very different window "
                            "lengths across the documents. Ignored if '-z' not provided.", metavar='length', default=0, dest='window')
    run_parser.add_argument('-f', type=int, help="Size of the document fragments (default: use full documents).",
                            metavar='size', default=0, dest='fragment')
    run_parser.add_argument('-F', type=int, help="Size of the author slices (default: use full authors). Usually this length should be at least ten times the fragments length.",
                            metavar='size', default=0, dest='authorLength')
    run_parser.add_argument('-s', type=str, metavar='Name', dest="suffix",
                            help="Suffix to folder name to keep order in the results. If left blank a random unique 5 char string is added during output directory creation.")
    run_parser.add_argument('-c', type=str, help="Path to the config file if not using a global configuration.", metavar='path', dest='configFile', default="")
    run_parser.add_argument('-p', type=str, help="Path to the file storing P0.", metavar='path', dest='P0file')
    run_parser.add_argument('-D', type=float, nargs='+', dest="delta", default=1, help="Values of delta for attribution.")
    group_run2 = run_parser.add_mutually_exclusive_group()
    group_run2.add_argument('-l', type=int, nargs='?', const=1, metavar='leave-out', default=0, dest='leaveNout',
                            help="Size of the slices used to split the corpus (default 1: leave-one-out).")
    group_run2.add_argument('-L', type=int, nargs='?', const=10, metavar='folds', default=10, dest='folds',
                            help="Number of folds for cross-validation (default: 10 also if both 'l' and 'L' missing).")
    group_run3 = run_parser.add_mutually_exclusive_group()
    group_run3.add_argument('-S', action='store_false', default=True, dest='mantainSlicing',
                            help="Always create new slicing of the database in folds.")
    group_run3.add_argument('-M', action='store_true', default=True, dest='mantainSlicing',
                            help="Mantain previous slicing of the database in folds.")
    group_run4 = run_parser.add_mutually_exclusive_group()
    group_run4.add_argument('-o', action='store_false', default=False, dest='overwriteOutputs',
                            help="Keep all the useful files in the results folder.")
    group_run4.add_argument('-O', action='store_true', default=False, dest='overwriteOutputs',
                            help="Always overwrite results folder.")
    run_parser.add_argument('-K', action='store_true', dest="keepTemporary",
                        help="Keep temporary files in the after computation. Useful if they are expensive and you want to retrieve '-R' them later.")
    run_parser.add_argument('-R', action='store_true', dest="retrieve",
                        help="Retrieve useful files from other folders.")
    run_parser.add_argument('-P', action='store_true', dest="allowPartial",
                        help="Allow partial attributions.")
    run_parser.add_argument('-A', action='store_true', dest="authorFiles",
                        help="Each input file contains all the documents of an author.")
    run_parser.add_argument('-V', action='store_true', dest="margOut",
                        help="Dump file with detailed results.")
    run_parser.add_argument('-d', action='store_true', dest="dumpP0",
                        help="Dump binary file with token hashes and multiplicities.")
    run_parser.add_argument('-v', action='store_true', dest="sliceSeparated",
                        help="Print separate statistics for every validation fold.")
    run_parser.add_argument('-a', type=str, help="Path to the author association file.", metavar='path', dest='association', default={})
    run_parser.set_defaults(func=elaboration.from_command_line)

    prepare_parser = subparsers.add_parser(
        'prepare', help="Prepares a properly formatted corpus for further processing.", description="The document file name should be in the format UNIQUEAUTHORNAME_futile_fields_BOOK_NAME.ext\nThe separator '_' can be replaced with others using option '-s', the extension can be changed with option '-e'.\nThe number of futile fields (plus the author name) can be specified throug option 't'.")
    prepare_parser.add_argument('--input', '-i', type=str, dest='inFolder',
                                help="Input folder containing the corpus.", default='', metavar='path')
    prepare_parser.add_argument('--output', '-o', type=str, dest='outFolder',
                                help="Output folder for the processed corpus.", default='', metavar='path')
    prepare_parser.add_argument('-A', action='store_true', dest="authorFiles",
                                help="If set, stores all the documents of each author in a single file.")
    prepare_parser.add_argument('--sep', '-s', type=str, dest='sep', default='_',
                                help="Field separator in the documents' file name (default: '_'). More than one character is accepted.")
    prepare_parser.add_argument('--ext', '-e', type=str, dest='extension', default='txt',
                                help="Extension of the documents' files (default: 'txt').")
    prepare_parser.add_argument('--title-part', '-t', type=int, dest='bookStart', default=2, metavar='N',
                                help="First field of the book title (counting from 0, deafult:2).")
    prepare_parser.set_defaults(func=support.data_preprocessing)

    config_parser = subparsers.add_parser(
        'config', help="Create a configuration file.")
    group_conf = config_parser.add_mutually_exclusive_group()
    group_conf.add_argument('--global', '-g', action='store_true',
                            help="Create a global configuration file.")
    group_conf.add_argument('--output', '-o', type=str, help="Output folder for the config file.",
                            default='', metavar='path', dest='outFolder')
    config_parser.set_defaults(func=support.build_config,
                            globalDir=os.getcwd())

    shuffle_parser = subparsers.add_parser(
        'shuffle', help="Shuffle a corpus document-author association for noisy test.")
    shuffle_parser.add_argument('--input', '-i', type=str,
                                help="Input folder containing the corpus.", default='', metavar='path', dest='inFolder')
    group_shuff = shuffle_parser.add_mutually_exclusive_group()
    group_shuff.add_argument('--preserve', '-p', action='store_true',
                            help="Preserve the distribution of books per author.", default=None, dest='preserveDist')
    group_shuff.add_argument('--no-preserve', '-n', action='store_false',
                            help="Do not preserve the distribution of books per author.", default=None, dest='preserveDist')
    shuffle_parser.set_defaults(func=support.corpus_shuffler)

    args = parser.parse_args()
    args.func(**vars(args))