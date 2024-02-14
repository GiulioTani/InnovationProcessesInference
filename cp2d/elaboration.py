#!/home/giulio/miniconda3/bin/python3

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

# @file elaboration.py
# @author Giulio Tani Raffaelli (tani@cs.cas.cz)
# @brief Contains the python elaboration for the attribution of books to authors.
# @version 1.0
# @date 2023-06-27
#
# @copyright Copyright (c) 2023
#

# @package elaboration
# @brief Contains the python elaboration for the attribution of books to authors.
#
# All the functions in this module, except for data_preprocessing and build_config, are meant to be called exclusively from another script, i.e. the workflow.
# These functions perform all the high level elaboration on the files: management of rough texts, standard optimization of parameters and so on up to the extraction of the results from the computed conditional probabilities.
#


import json
import logging
import multiprocessing as mp
import os
import shutil
import sys
from glob import glob
import numpy as np
import pandas as pd
import typing as tp
import time

from . import support, logic, numcomp


logging.basicConfig(filename=os.path.join(os.path.dirname(__file__), 'cp2dExperiment.log'), filemode='w', level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger()


class cp2dExperiment ():
    __noauth = int(0x02)
    __nofrag = int(0x01)
    __auth_t = 'u2'
    __slice_t = 'u4'
    __book_t = 'u4'
    __frag_t = 'u4'
    __tok_t = 'u4'

    def __init__(self, LZ77: bool = False, fragment: int = -1, authorLength: int = 0, window: int = -1, suffix: str = None, configFile: str = "", database: str = "", ngramSize: int = 0, leaveNout: int = 0, mantainSlicing: bool = True, allowPartial: bool = False, retrieve: bool = False, keepNonAlfabetical: bool = False, overwriteOutputs: bool = False, keepTemporary: bool = False, delta: float = 1, folds: int = 10, authorFiles: bool = False, goodSlices = [], P0file:str="", dumpP0:bool=False, **kwargs) -> None:
        """
        Initializes the parameters for the analysis from command line options.

        Parses the command line arguments to build all the parameters needed. Created the name of the output folder and the folder if missing.
        If the output folder exists prompt for overwriting or resuming the previous analysis.
        Possible options are:
            -# 'z' to use LZ77 dictionaries (incompatible with 'n' and 'U');
            -# 'f' [fragment length] size of the book fragments;
            -# 'F' [author length] size of the virtual authors to compare on, usually >~10 times fragment length;
            -# 'w' [window size] size of the window for LZ77;
            -# 's' suffix to folder name;
            -# 'c' config file name if different from \c config.ini
            -# 'd' path to database folder (e.g. the output folder chosen with data preprocessing)
            -# 'n' size of the overlapping space free n-grams (incompatible with 'z' and 'U');
            -# 'l' Size of the slices used to split the corpus (default 1: one book at time);
            -# 'S' Always create new Slicing;
            -# 'M' Mantain previous slicing (cancels any S option);
            -# 'P' Allow partial attributions;
            -# 'R' Always Retrieve information from other folders;
            -# 'U' Keep underscores and punctuation when using dictionary method (incompatible with 'z' and 'n')
            -# 'O' Always overwrite folder (cancel any R option).
            -# 'K' Keep probFile and wnt folder after computation.
            -# 'b' Directory where of the base corpus if using processed versions.
            -# 'B' Base size of the fragments, only useful for LZ77 processing.
        The path to the config file may be absolute or relative and is searched in the cwd, then in data, then in res and finally in res/{dataDir}.
        In case of badly written config file tries to deduce some of the parameters from the system or the values of others, if fails falls back to the default value.

        Args:
            read: A list of (option, value) pairs in the format returned by getopt.

        Returns:
            param: A dictionary with the parameters needed to carry on the analisys.

        @sa getopt, workflow.read
        """
        # global logic._Q_association, logic._Q_excluded, logic._Q_which_slice, logic._Q_apt, logic._Q_fra, logic._Q_partList
        parameters = {key: value for key,
                      value in locals().items() if key != "self"}

        self.__param = support.lockable_dict(dict=parameters)
        self.__attributions = None
        self.__unk = None
        if hasattr(logic, "_Q_lock") and logic._Q_lock is not None:
            raise RuntimeError(
                "Only one instance of exp at the time is allowed.")
        logic._Q_lock = True
        try:
            logic._Q_association = None
            logic._Q_excluded = None
            logic._Q_which_slice = None
            logic._Q_apt = None
            logic._Q_fra = None
            logic._Q_partList = None
            logic._Q_F = None
            self.__numSlices = None
            logic.check_parameters(self.__param)
            self.__param['database']=os.path.abspath(self.__param['database'])
            if self.__param['configFile']:
                self.__param['configFile']=os.path.abspath(self.__param['configFile'])
            logic._Q_F = self.__param['fragment']
            self.workers, self.Encoding, self.__param['configFile']= logic.read_config(
                configPath=self.__param['configFile'], workers=4, Encoding="latin1")

            if os.path.isfile(os.path.join(self.__param['database'], "encoding.dat")):
                with open(os.path.join(self.__param['database'], "encoding.dat")) as fp:
                    self.Encoding = fp.readline()
            self.__foundResults = False
            if not self.__param['LZ77']:  # without the LZ77 is always CopyFrag
                self.__param['window'] = -1

            self.__param['databaseCtime'] = os.path.getctime(
                self.__param['database'])

            self.el_path = os.path.abspath(
                os.path.join(os.path.dirname(__file__), '..'))
            self.__param["sourceInfo"] = logic.check_sources(self.el_path)

            if folds and not leaveNout:
                self.__param["leaveNout"] = logic.foldsize(
                    database, folds, authorFiles)
            self.__param["baseDir"], self.__param["authDir"], self.__param["fragDir"] = logic.dirnames(
                self.__param)
            self.__param["resultsName"] = "-".join(
                [self.__param["baseDir"], self.__param["authDir"], self.__param["fragDir"]])
            self.__param['resultsDir'] = os.path.abspath(os.path.join(
                self.el_path, "../res", self.__param["baseDir"], self.__param["authDir"], self.__param["fragDir"]))
            self.__param['dataDir'] = os.path.abspath(os.path.join(
                self.el_path, "../res", self.__param["baseDir"]))

            self.__param["TS"] = time.time()
            self._create_experiment_directory()
            with open(os.path.join(self.__param['dataDir'], "parameters.json"), "a") as fp:
                print(json.dumps({k: v for k, v in self.__param.items()
                           if k != "kwargs"}, default=support.numpy_json), end="\n", file=fp)

            if self.__param['retrieve'] and not self.__foundResults:
                self.retrieve()
        except:
            logic._Q_lock = None
            raise

    def __del__(self):
        # global logic._Q_association, logic._Q_excluded, logic._Q_which_slice, logic._Q_apt, logic._Q_fra, logic._Q_partList
        logic._Q_association = None
        logic._Q_excluded = None
        logic._Q_which_slice = None
        logic._Q_apt = None
        logic._Q_fra = None
        logic._Q_partList = None
        logic._Q_F = None
        logic._Q_lock = None

    def run(self):
        if not self.__foundResults:
            logging.info("Extracting fragments")
            self.fragments()
            logging.info("Preparing support files")
            self.support_files()
            logging.info("Computing probabilities")
            print("Computing probabilities. . .", end="\r")
            self.compute_probabilities()

    def runall(self):
        self.run()
        logging.info("Computing results")
        self.results(self.__param['delta'], self.__param['allowPartial'], nonAttri=[
        ], excluded=[], association={}, margOut=True)
        self.clean()

    @property
    def resultsDir(self):
        return self.__param['resultsDir']

    @resultsDir.setter
    def resultsDir(self, __):
        raise AttributeError("Can't change directory structure by hand!!")

    @property
    def resultsName(self):
        return self.__param['resultsName']

    @resultsName.setter
    def resultsName(self, __):
        raise AttributeError("Can't change directory structure by hand!!")

    @property
    def dataDir(self):
        return self.__param['dataDir']

    @dataDir.setter
    def dataDir(self, value):
        try:
            if os.path.isdir(value):
                raise FileExistsError(f"Directory {value} already existing.")
            oldDir = self.__param['dataDir']
            shutil.move(os.path.join(oldDir), os.path.join(value))
            self.__param['dataDir'] = value
        except ValueError as e:
            print(e)

    @property
    def Encoding(self):
        return self.__encoding

    @Encoding.setter
    def Encoding(self, value):
        try:
            "a".encode(value)
        except LookupError as e:
            print(e)
            self.__encoding = None
            return
        self.__encoding = value

    @property
    def unknowns(self):
        return self.__unk

    @unknowns.setter
    def unknowns(self, __: tp.Any) -> None:
        raise AttributeError("Can't change attributions by hand!!")

    @property
    def attributions(self):
        return self.__attributions

    @attributions.setter
    def attributions(self, __: tp.Any) -> None:
        raise AttributeError("Can't change attributions by hand!!")

    @property
    def whichSlice(self) -> tp.Mapping[tp.Sequence[int], int]:
        if logic._Q_which_slice is None:
            print("Loading slices")
            self.__numSlices = logic.loadSlices(os.path.join(self.__param['resultsDir'], "slices.bin"), self.__slice_t, self.__auth_t, self.__book_t)
        return logic._Q_which_slice

    @whichSlice.setter
    def whichSlice(self, slices: tp.Mapping[tp.Sequence[int], int]) -> None:
        logic._Q_which_slice = slices
        self.__numSlices = len(set(logic._Q_which_slice.values()))
        raise DeprecationWarning(
            "Manually setting the slicing is prone to errors and therefore deprecated.")

    @property
    def sliceNum(self) -> int:
        return self.__numSlices

    @sliceNum.setter
    def sliceNum(self, __: tp.Any) -> None:
        raise AttributeError("Can't change number of slices by hand!!")

    def _create_experiment_directory(self):
        if os.path.isdir(self.__param['dataDir']):
            logger.info("Found existing folder")
            old_different = False
            if os.path.isfile(os.path.join(self.__param['dataDir'], "parameters.json")):
                with open(os.path.join(self.__param['dataDir'], "parameters.json")) as fp:
                    atLeastOneLine = False
                    for line in fp:
                        atLeastOneLine = True
                if atLeastOneLine:
                    try:
                        old_param = json.loads(line)
                        comp_string = []
                        cause = []
                        if old_param['database'] != self.__param['database']:
                            comp_string.append(
                                f"database now: {self.__param['database']}, was: {old_param['database']}")
                            old_different = True
                            cause.append("different path")
                        if 'databaseCtime' not in old_param or old_param['databaseCtime'] != self.__param['databaseCtime']:
                            comp_string.append(
                                "Databases may have been edited since results were computed.")
                            old_different = True
                            cause.append("edited corpus")
                        if 'sourceInfo' in old_param:
                            if self.__param['sourceInfo']["warnSourceNewer"] or old_param['sourceInfo']["warnSourceNewer"]:
                                old_different = True
                                cause.append("source newer")
                                comp_string.append(
                                    "Impossible to check source at compile time, source is newer.")
                            parts = []
                            for k in self.__param['sourceInfo']:
                                if self.__param['sourceInfo'][k] != old_param['sourceInfo'][k]:
                                    old_different = True
                                    parts.append(k)
                                    cause.append("different sourceflags")
                            comp_string.append(
                                "Different sources on: "+", ".join(parts))
                    except:
                        with open(os.path.join(self.__param['dataDir'], "parameters.json"), "a") as fp:
                            print(file=fp)
                        old_param = None
                        old_different = True
                        cause = ["missing parameters"]
                else:
                    old_param = None
                    old_different = True
                    cause = ["missing parameters"]
            else:
                old_param = None
                old_different = True
                cause = ["missing parameters"]

            if old_different:
                if not support.get_input("A folder with the same name is present but corpora seem to be different.\n"
                                                   f"{'The parameters used are missing.' if old_param is None else ' '.join(comp_string)}\n"
                                                   f"Do you want to IGNORE this difference and treat as good those results?"):
                    raise FileExistsError(
                        f"A directory with the same name and different parameters ({' - '.join(cause)}) is already existing.\nTry changing the suffix.")

            if os.path.isdir(self.__param['resultsDir']):
                if self.__param['overwriteOutputs']:
                    logger.info("Erasing existing folder")
                    print("Erasing. . .")
                    shutil.rmtree(self.__param['resultsDir'])
                    os.mkdir(self.__param['resultsDir'])
                    logger.info("Created folder")
                else:
                    inName = os.path.join(
                        self.__param['resultsDir'], self.__param["resultsName"])
                    if len(glob(os.path.join(inName+"_*.res"))) or (os.path.isfile(inName+"_ids.res") and os.path.isfile(inName+"_pro.res")):
                        logger.info("Found existing results")
                        self.__foundResults = True
                        if not self.__param.islocked():
                            self.__param.lock()
            else:
                os.makedirs(self.__param['resultsDir'])
                logger.info("Created folder")
        else:
            # raise FileNotFoundError(f"This experiment does not exist! {self.__param['resultsDir']}")
            os.makedirs(self.__param['resultsDir'])
            logger.info("Created folder")

    def _author_processor(self, inputDir, authName):
        extension = "seq"
        anum = int(authName[1:])
        booklist = []
        outputFileName, outputFileTmpName = logic.output_file_names(
            self.__param['dataDir'], authName, extension)

        if os.path.isfile(outputFileName):
            logger.debug(f"Author {authName}, already there")
            return booklist

        logger.debug(f"Author {authName} started")
        with open(outputFileTmpName, "w", encoding=self.Encoding) as f_out:
            for bn, text_b in support.auth_books_iterator(os.path.join(inputDir, authName), self.Encoding, self.__param["authorFiles"]):
                booklist.append((anum, int(bn)))
                fragproc = logic.seq_to_wnt(
                    fragment=text_b, compression=self.__param["LZ77"], window=self.__param[
                        "window"], Encoding=self.Encoding, translation=self.__translation)
                logic.write_sequence(
                    bookNumber=bn, sequence=fragproc, file=f_out)

        logger.debug(f"Author {authName} done")

        if os.path.isfile(outputFileTmpName):
            os.rename(outputFileTmpName, outputFileName)

        return booklist

    def fragments(self):
        """
        Splits texts in fragments.

        This is Margherita's function, the only edit I made were to avoid
        useless computation and remove prints. Is not safe but it will fail unlikely.

        If baseFragLength is greater than zero all the fragments generated will have the same length, shorter parts are discarded.

        Converts sequence files to wnt.

        The fragments are converted to wnt and encoding unified to utf-8.
        The original length of the file (in bytes) is lost but shall give
        less problems.

        Args:
            dataDir:    Path to the directory containing the fragments to be converted.
            useEncoding:    Encoding of the sequence files, not needed if already specified calling other functions from this module.


        Args:
            inputDir:   Directory to look for texts to split.
            outputDir:  Directory that will contain the fragments.
            baseFragLength: If set to 0 only a single fragment containing the whole book is created. Otherwise zero or more fragments of length baseFragLength are created.
            useEncoding:    Encoding of the sequence files, not needed if already specified calling other functions from this module.

        Raises:
            ValueError
        """
        if not self.__param.islocked():
            self.__param.lock()

        logger.info("Extracting fragments.")
        # prepare the right translation
        self.__translation = logic.derive_translation(
            LZ77=self.__param["LZ77"], ngramSize=self.__param["ngramSize"], keepUnder=self.__param["keepNonAlfabetical"])

        outputDir = logic.create_tmpfile_dir(
            dataDir=self.__param["dataDir"], tmpKind="seq")

        if self.__param["LZ77"]:
            if os.path.isfile(os.path.join(self.__param["database"], "encoding.dat")):
                shutil.copy(os.path.join(self.__param["database"], "encoding.dat"),
                            os.path.join(outputDir, "encoding.dat"))
        if os.path.isfile(os.path.join(self.__param["database"], "slices.bin")) and self.__param['mantainSlicing']:
            shutil.copy(os.path.join(self.__param["database"], "slices.bin"),
                        os.path.join(outputDir, "slices.bin"))

        self._parallel_extraction(self.__param["database"])

    def _parallel_extraction(self, inputDir):
        # prepare extraction
        extraction_pool = mp.Pool(self.workers)
        extraction_jobs = []
        booklist = []
        fileList = [f for f in os.listdir(inputDir) if f.endswith(".txt")]
        logger.info(f"Begin extraction from {len(fileList)} books")
        if self.__param["authorFiles"]:
            authList = [os.path.splitext(f)[0] for f in fileList]
        else:
            authList = set(f.split("B")[0] for f in fileList)
        # extract wnt or seq files
        for authName in authList:
            extraction_jobs.append(extraction_pool.apply_async(self._author_processor, [
                inputDir, authName], callback=booklist.extend, error_callback=lambda x: print(x)))

        support.progress_wait(extraction_jobs)
        logger.info(
            f"Extracted {len(fileList)} books.")
        print(
            f"\rExtracted {len(fileList)} books.")
        print(f"Done {len(extraction_jobs)} authors")

        extraction_pool.close()
        extraction_pool.join()

    def support_files(self):
        """
        Generates the support files used by fragprob.

        Produces the fragments and params support files.
            -# Fragments contains a list of books with the relative author and the number of fragments of which each book is composed.
            -# Params contains the \c Alpha and \c Theta parameters for all the authors.

        It is possible to manually change the Fragments file to exclude books from the analysis or limit the number of fragments taken into consideration.
        Even if an author is missing or removed from the analysis her line must remain in the params file.

        In case of error in the computation of the parameters for an author gives a warning an then uses the average from the other authors as parameters.
        It doesn't raise an actual warning as the default behaviour is to display them only once. To be modified in the future.

        Args:
            dataDir:    Path to a directory containing the subdirectories seq and wnt (files of the fragment for the analysis with utf-8 encoding).
            mantainSlice: True keep slicing, False generate a new one, None (default) ask.
        """
        if not self.__param.islocked():
            self.__param.lock()

        if os.path.isfile(os.path.join(self.__param['dataDir'], "seq", "slices.bin")):
            logger.info(f"Found slicing")
            print(
                "Found slicing of the corpus (possibly with related parameters).")
            if not self.__param['mantainSlicing']:
                logger.info(f"Removing slicing")
                os.remove(os.path.join(
                    self.__param['dataDir'], "seq", "slices.bin"))
                for auth in os.listdir(self.__param['dataDir']):
                    if os.path.isdir(os.path.join(self.__param['dataDir'], auth)):
                        if os.path.isfile(os.path.join(self.__param['dataDir'], auth, "params.bin")):
                            os.remove(os.path.join(
                                self.__param['dataDir'], auth, "params.bin"))
                print("Slicing removed.")

    def compute_probabilities(self):
        if not self.__param.islocked():
            self.__param.lock()
        try:
            numcomp.init(self.__param['configFile'])
            # prepare common files
            if os.path.isfile(os.path.join(self.__param['dataDir'], self.__param['authDir'], "params.bin")):
                shutil.copy2(os.path.join(self.__param['dataDir'], self.__param['authDir'], "params.bin"), os.path.join(
                    self.__param['dataDir'], "seq"))

            # compute
            numcomp.calprob(self.__param['dataDir'], self.__param['resultsDir'], outFile=self.__param["resultsName"], slicesize=self.__param['leaveNout'],
                            ngram=self.__param['ngramSize'], fragment=self.__param['fragment'], authsize=self.__param['authorLength'], goodSlices=self.__param['goodSlices'], P0file=self.__param['P0file'], dumpP0=self.__param['dumpP0'])
        except RuntimeError:
            sys.exit(1)
        finally:
            # clean common files
            # add params to authDir
            if not os.path.isfile(os.path.join(self.__param['dataDir'], self.__param['authDir'], "params.bin")) and os.path.isfile(os.path.join(self.__param['resultsDir'], "params.bin")):
                shutil.copy2(os.path.join(self.__param['resultsDir'], "params.bin"), os.path.join(
                    self.__param['dataDir'], self.__param['authDir']))
            # remove params from seq
            if os.path.isfile(os.path.join(self.__param['dataDir'], "seq", "params.bin")):
                os.remove(os.path.join(
                    self.__param['dataDir'], "seq", "params.bin"))
            # add slices to seq
            if not os.path.isfile(os.path.join(self.__param['dataDir'], "seq", "slices.bin")) and os.path.isfile(os.path.join(self.__param['resultsDir'], "slices.bin")):
                shutil.copy2(os.path.join(self.__param['resultsDir'], "slices.bin"), os.path.join(
                    self.__param['dataDir'], "seq"))

        self.__foundResults = True

    def juicer(self, excluded=None, association=None):
        if excluded is None or association is None:
            partial = logic.availableResults(self.__param['resultsDir'], self.__param['database'])
            if excluded is not None:
                partial = [(exc, asso) for exc, asso in partial if exc==frozenset(excluded)]
            if association is not None:
                partial =[(exc, asso) for exc, asso in partial if asso==tuple(sorted(association.items()))]
            print("juicing:", partial)
            final = [(list(exc), {k:v for k, v in asso}) for exc, asso in partial]
            for exc, asso in final:
                self.__juicer(exc, asso)
        else:
         self.__juicer(excluded, association)

        with open(os.path.join(self.__param['resultsDir'], self.__param["resultsName"]+"_RESULTSGUARD.res"), "w") as fp:
            pass

    def __juicer(self, excluded=[], association={}):
        """
        Computes the results of the attribution in terms of success percentages.

        From the file of the conditional probabilities of every fragment against every author computes the percentage of successful attributions.
        It also computes the attributions excluding one author in case a single one biases the analysis.

        Args:
            excluded:   Names in the format int(author_number) of the authors that must completely ignored (no comparison nor attribution).
            association:Dictionary whose keys and values are author numbers. This allows to consider some author as a different
                        one (even if the parameters used in the estimate of probabilities are those of the formerly selected) both when assigning 
                        its books and when considering it in assigning others. Association is performed prior to exclusion thus this allows also
                        the exclusion of selected authors in group by assigning them to a dummy author then excluded. Each author can be assigned
                        only once.
        """
        # global logic._Q_association, logic._Q_excluded, logic._Q_which_slice, logic._Q_apt, logic._Q_fra, logic._Q_partList, logic._Q_proFile, logic._Q_topRes
        if not hasattr(excluded, '__iter__'):
            logic._Q_excluded = [excluded, ]
        else:
            logic._Q_excluded = excluded
        logic._Q_topRes = -2

        logic._Q_association = logic.loadAssociation(association, self.__param['database'])
        if "_association_" in logic._Q_association:
            raise logic._Q_association["_association_"]

        i = 0
        out_name = os.path.join(self.__param['resultsDir'], self.__param["resultsName"]+f"_JUICED{i:02}.csv")
        while os.path.isfile(out_name):
            i+=1
            out_name = os.path.join(self.__param['resultsDir'], self.__param["resultsName"]+f"_JUICED{i:02}.csv")

        logger.info(f"Assigning books.")


        ret, self.__unk = logic.previous_results(
            self.__param['resultsDir'], [],  [], False, allDeltas=True)
        self._logDelta = logic.selectDelta(ret)

        logic._Q_missingDelta = {self._logDelta[0],0}

        logic._Q_proFile = os.path.join(
            self.__param['resultsDir'], self.__param["resultsName"]+"_pro")

        try:
            self.prepare_for_new_results(True, True)
        except FileNotFoundError as e:
            raise FileNotFoundError(*e.args, "While preparing results")
        logger.info("Computing Attributions")
        chunksize = min(
            int(len(logic._Q_partList)/(self.workers*2)+1), 2000)
        with mp.Pool(self.workers) as p:
            # clearing it beforehand gives a huge speedup on second call to results
            self.__attributions = None
            self.__attributions = {ab: tmp for ab, tmp in support.tqdm(p.imap_unordered(logic.attributor, logic._Q_partList, chunksize=chunksize), total=len(
                logic._Q_partList), desc="Attribution ", dynamic_ncols=True, disable=None, file=sys.stderr, leave=True)}
        logger.info("Assigning texts")
        assigned = logic.assign(self.__attributions )

        if (self.__flagbyte & self.__noauth):
            assigned.drop(columns=["TOP", "WP", "TMR", "WMR", "FRA_TMR", "FRA_WMR"],level=1, inplace=True)

        
        with open(out_name, 'w') as fp:
            print('#',json.dumps({"excluded":excluded, "association":association}), file=fp)
        assigned.to_csv(out_name,mode='a')



    def results(self, delta=None, allowPartial=False, nonAttri=[], excluded=[], association={}, margOut=False, machine=False, PANStyle="strict", groundTruth=None, othResult=None, forceRecompute=False, sliceSeparated=False, topres=2, saveUnk=True):
        """
        Computes the results of the attribution in terms of success percentages.

        From the file of the conditional probabilities of every fragment against every author computes the percentage of successful attributions.
        It also computes the attributions excluding one author in case a single one biases the analysis.

        Args:
            dataDir:    Path to the directory containing the results.
            nonAttri:   Names in the format int(author_number) of the authors whose books are not attributed (but used for comparisons).
            excluded:   Names in the format int(author_number) of the authors that must completely ignored (no comparison nor attribution).
            association:Dictionary whose keys and values are author numbers. This allows to consider some author as a different
                        one (even if the parameters used in the estimate of probabilities are those of the formerly selected) both when assigning 
                        its books and when considering it in assigning others. Association is performed prior to exclusion thus this allows also
                        the exclusion of selected authors in group by assigning them to a dummy author then excluded. Each author can be assigned
                        only once.
            resName:    The name of the file containing the probabilities if not named as <dataDir>.res.
            margOut:    Print to file detailed information for the attribution af every book and author. Default is False.
            machine:    Whether to silence outputs and provide results as a dictionary for automated machine use.
            allowPart:  Allows partial attribution of books if more tha one author is proposed as first candidate. Relevant only for MR and mostly
                        when books have few fragments.
        """
        # global logic._Q_association, logic._Q_excluded, logic._Q_which_slice, logic._Q_apt, logic._Q_fra, logic._Q_partList, logic._Q_proFile, logic._Q_topRes
        if delta is None:
            delta = self.__param["delta"]
        if not hasattr(delta, '__iter__'):
            delta = [delta, ]
        if not hasattr(nonAttri, '__iter__'):
            self.__nonAttri = [nonAttri, ]
        else:
            self.__nonAttri = nonAttri
        if not hasattr(excluded, '__iter__'):
            logic._Q_excluded = [excluded, ]
        else:
            logic._Q_excluded = excluded
        self._logDelta = np.log10(delta)
        if not np.isfinite(self._logDelta).all():
            raise ValueError(f"Invalid delta: {delta}")
        if topres is not None:
            assert int(topres) > 0, "At least one result must be required"
            logic._Q_topRes = -1-int(topres)
        else:
            logic._Q_topRes = None
        self.__unk = {}
        logic._Q_association = logic.loadAssociation(association, self.__param['database'])
        if "_association_" in logic._Q_association:
            raise logic._Q_association["_association_"]

        slSep = None
        logger.info(f"Assigning books.")
        ret, self.__unk = logic.previous_results(
            self.__param['resultsDir'], self.__nonAttri,  self._logDelta, allowPartial)
        logic._Q_missingDelta = [
            delta for delta in self._logDelta if forceRecompute or delta not in ret]
        logic._Q_proFile = os.path.join(
            self.__param['resultsDir'], self.__param["resultsName"]+"_pro")
        if logic._Q_missingDelta:
            self.prepare_for_new_results(machine, sliceSeparated)
            logger.info("Computing Attributions")
            chunksize = min(
                int(len(logic._Q_partList)/(self.workers*2)+1), 2000)
            with mp.Pool(self.workers) as p:
                # clearing it beforehand gives a huge speedup on second call to results
                self.__attributions = None
                self.__attributions = {ab: tmp for ab, tmp in support.tqdm(p.imap_unordered(logic.attributor, logic._Q_partList, chunksize=chunksize), total=len(
                    logic._Q_partList), desc="Attribution ", dynamic_ncols=True, disable=None, file=sys.stderr, leave=True)}
            self.authors = list(
                set(ab[0] for ab in self.__attributions if ab[0] > 0))
            logger.info("Assigning texts")
            tret, slSep = logic.return_dict(
                self.__attributions, self.__nonAttri, allowPartial, self.authors, sliceSeparated=self.__numSlices if (sliceSeparated and len(self.__param['goodSlices'])!=1) else False, goodSlices=self.__param['goodSlices'])
            scoreDelta = [(d,tret[d]['all']['FNN']["weigh"]["R"]) for d in tret]
            bestDelta = max(scoreDelta, key=lambda x: x[1])[0]
            self.__unk.update(logic.assign_unknown(self.__attributions, bestDelta))
            for nowDelta in logic._Q_missingDelta:
                if margOut:
                    logger.info(f"Creating margOut files")
                    logic.print_margout(
                        self.__param['resultsDir'], self.__attributions, self.__nonAttri, nowDelta, allowPartial, tret, nowDelta)

                if not machine and not self.__param['goodSlices']:
                    print(f"#########\n     {10**nowDelta}\n#########")
                    support.PAN11({ab: self.__attributions[ab][nowDelta] for ab in self.__attributions}, PANStyle,
                                  groundTruth, othResult, self.__numSlices if sliceSeparated else False, logic._Q_which_slice, self.authors)
                if nowDelta not in ret:
                    logic.add_results_line(self.__param["dataDir"], tret[nowDelta], self.__nonAttri, allowPartial, nowDelta,
                                           self.__unk[nowDelta], saveUnk=saveUnk, F=self.__param['fragment'], A=self.__param['authorLength'], goodSlices=self.__param['goodSlices'])
                    logic.add_results_line(
                        self.__param["resultsDir"], tret[nowDelta], self.__nonAttri, allowPartial, nowDelta, self.__unk[nowDelta], saveUnk=saveUnk, goodSlices=self.__param['goodSlices'])
            ret.update(tret)
        for nowDelta in self._logDelta:
            if not machine:
                print(f"#########\n     {10**nowDelta}\n#########")
                print("   ", *list(ret[nowDelta].keys()), sep="\t")
                for which in ret[nowDelta]:
                    print(f'\n\n{which}:\n\t\tWeigh\t\t\tMacro\n\tP\tR\tF1\tP\tR\tF1\t', end="\n")
                    for r in ret[nowDelta][which].items():
                        if r[0] == "FRA":
                            print(f"{r[0]}\t     \t{r[1]*100:.4}%", end="\n")
                        else:
                            print(r[0], *(f"{r[1]['weigh'][measure]*100:.4}%" for measure in ['P', 'R', 'F']),
                                    *(f"{r[1]['macro'][measure]*100:.4}%" for measure in ['P', 'R', 'F']), sep="\t")
                if sliceSeparated and slSep is not None:
                    for avg in ["weigh", "macro"]:
                        for measure in ["P", "R", "F"]:
                            print(avg, measure)
                            print("S",*[typ for typ in slSep[nowDelta][0]["all"]], sep="\t")
                            for i, slres in enumerate(slSep[nowDelta]):
                                if self.__param['goodSlices']:
                                    print(self.__param['goodSlices'][i], end='\t')
                                else:
                                    print(i, end='\t')
                                print(*[round(slres['all'][typ][avg][measure], 4) if typ!="FRA" else (round(slres['all'][typ], 4) if (avg, measure)==("weigh", "R") else "---")
                                    for typ in slres['all']], sep="\t")

        returnRet = None
        if machine:
            if sliceSeparated:
                returnRet = (ret, slSep)
            else:
                returnRet = ret
        return returnRet

    def prepare_for_new_results(self, machine: bool, sliceSeparated: bool = False):
        # global logic._Q_association, logic._Q_excluded, logic._Q_which_slice, logic._Q_apt, logic._Q_fra, logic._Q_partList
        inName = os.path.join(
            self.__param['resultsDir'], self.__param['resultsName'])
        if not machine:
            print("Loading results. . .", end="\r")

        if logic._Q_fra is None:
            with open(inName+"_fra.res", "rb") as fin:
                byte = fin.read(1)
                self.__flagbyte = int.from_bytes(byte, "little")
                dtype = [('a', self.__auth_t), ('b', self.__book_t)]
                if self.__flagbyte & self.__nofrag:
                    dtype.append(('f_siz', self.__tok_t))
                    logic._Q_fra = pd.DataFrame(np.fromfile(fin, dtype=dtype)).set_index(
                        ["a", "b"]).sort_index(level='a')
                else:
                    fra_rows = []
                    dtype.append(('fn', self.__frag_t))
                    try:
                        while True:
                            fra_rows.append(
                                list(np.fromfile(fin, dtype=dtype, count=1)[0]))
                            fra_rows[-1].append(np.fromfile(fin,
                                                dtype=self.__tok_t, count=fra_rows[-1][-1]))
                    except IndexError as e:
                        if fin.read() != b'':
                            raise e
                    logic._Q_fra = pd.DataFrame(fra_rows, columns=["a", "b", "fn", "f_siz"]).set_index(
                        ["a", "b"]).sort_index(level='a')
                bookSet = logic._Q_fra.index.tolist()

        if logic._Q_apt is None:
            if not (self.__flagbyte & self.__noauth):
                logic._Q_apt = {}
                if os.path.isfile(inName+"_aut.res"):
                    with open(inName+"_aut.res", "rb") as fin:
                        N_auth = np.fromfile(fin, dtype=self.__auth_t, count=1)[
                            0]  # number of authors
                        for __ in range(N_auth):
                            # author id, number of books, default size
                            AuNbDs = np.fromfile(fin, dtype=[
                                                 ('Au', self.__auth_t), ('Nb', self.__book_t), ('Ds', self.__book_t)], count=1)
                            Asizes = np.fromfile(
                                fin, dtype=[('s', self.__slice_t), ('fn', self.__book_t)], count=AuNbDs['Nb'][0])
                            logic._Q_apt[AuNbDs['Au'][0]] = [
                                {a['s']:a['fn'] for a in Asizes}, AuNbDs['Ds'][0]]
                else:
                    raise FileNotFoundError(f"Missing file {inName}_aut.res")
            else:
                logic._Q_apt = {a: [{}, 1]
                                for a in logic._Q_fra.index.get_level_values("a").unique()}
                if 0 in logic._Q_apt:
                    del logic._Q_apt[0]

        if logic._Q_which_slice is None and (sliceSeparated or not (self.__flagbyte & self.__noauth) or self.__param['goodSlices']):
            self.__numSlices = logic.loadSlices(os.path.join(self.__param['resultsDir'], "slices.bin"), self.__slice_t, self.__auth_t, self.__book_t)
                
        if (self.__flagbyte & self.__noauth) and not sliceSeparated:
            logic._Q_which_slice_is_None = True
        else:
            logic._Q_which_slice_is_None = False

        if logic._Q_partList is None:
            if len(glob(inName+"_ids*res")):
                logger.info(f"Loading results")
                # load ids
                if os.path.isfile(inName+"_ids.res"):
                    with open(inName+"_ids.res", "rb") as fin:
                        ids = np.fromfile(
                            fin, dtype=[('a', self.__auth_t), ('b', self.__book_t), ('tb', 'u4')], count=len(logic._Q_fra))
                else:
                    ids = []
                    for fids in sorted(glob(inName+"_ids*res")):
                        with open(fids, "rb") as fin:
                            ids.append(np.fromfile(
                                fin, dtype=[('a', self.__auth_t), ('b', self.__book_t), ('tb', 'u4')], count=-1))
                if len(ids)==0:
                    raise RuntimeError("Missing ids in file!")
            else:
                logger.error(f"Missing {inName}_* unable to compute results.")
                raise FileNotFoundError(f"Missing {inName}_*")

            booksAssociation = logic.apply_association(
                bookSet, logic._Q_association)

            if booksAssociation:
                newFraIndex = []
                for ind in logic._Q_fra.index:
                    if ind in booksAssociation:
                        newFraIndex.append(booksAssociation[ind])
                    else:
                        newFraIndex.append(ind)
                logic._Q_fra.index = pd.MultiIndex.from_tuples(
                    newFraIndex, names=["a", "b"])
                BAName = logic.sequential_file_name(
                    self.__param['dataDir'], "booksAssociation", "json")
                with open(BAName, "w") as fp:
                    baj = {str(ab): v for ab, v in booksAssociation.items()}
                    json.dump([logic._Q_association, baj], fp)

            logic._Q_partList = logic.fragments_numbers(
                ids, booksAssociation, logic._Q_excluded, logic._Q_fra)

    def clean(self):
        if not self.__param['keepTemporary']:
            print("Cleaning . . .")
            if os.path.isdir(os.path.join(self.__param['dataDir'], "seq")):
                logging.info("Cleaning seq directory")
                for seqfile in os.listdir(os.path.join(self.__param['dataDir'], "seq")):
                    if seqfile.endswith("seq") or seqfile in ["slices.bin", "params.bin"]:
                        os.unlink(os.path.join(
                            self.__param['dataDir'], "seq", seqfile))

    def delete_cached_results(self):
        for fname in glob(os.path.join(self.__param['resultsDir'], "*.res")):
            os.unlink(fname)
        self.__foundResults = False

    def retrieve(self):
        logger.info(f"Retrieving")
        print("Retrieving . . .", end="")
        # first try to get all at once
        if not os.path.isdir(os.path.join(self.__param['dataDir'], "seq")) or not len([a for a in os.listdir(os.path.join(self.__param['dataDir'], "seq")) if a.endswith("seq")]):
            self._retrieve_directoryes()

        if os.path.isdir(os.path.join(self.__param['dataDir'], "seq")):
            print(" done.")
        else:
            print(" failed. (no candidates)")

    def _retrieve_directoryes(self):
        pieces = self.__param['baseDir'].split("-")
        goodones = {}
        dirList = glob(os.path.join(self.el_path, f"../res/{pieces[0]}*"))
        dirList = [a for a in dirList if os.path.isdir(a)]
        if pieces[0] == "LZ77":
            dirList = [a for a in dirList if pieces[1] in a]
        for name in dirList:
            if not os.path.samefile(name, self.__param['dataDir']):
                if not os.path.isfile(os.path.join(name, "parameters.json")):
                    continue
                with open(os.path.join(name, "parameters.json")) as fp:
                    atLeastOneLine = False
                    for line in fp:
                        atLeastOneLine = True
                if atLeastOneLine:
                    try:
                        old_param = json.loads(line)
                    except json.decoder.JSONDecodeError:
                        continue
                    if self._is_good_candidate(old_param):
                        if os.path.isdir(os.path.join(name, "seq")):
                            goodones[name] = 2
                        if pieces[3] in name:
                            goodones[name] = goodones.get(name, 0) + 1
        if goodones:
            # selecting the most complete candidate
            chosen = sorted(goodones.keys(), key=lambda x: goodones[x])[-1]
            if os.path.isdir(os.path.join(chosen, "seq")):
                if os.path.isdir(os.path.join(self.__param['dataDir'], "seq")):
                    shutil.rmtree(os.path.join(self.__param['dataDir'], "seq"))
                print(f" from {os.path.basename(chosen)}", end="")
                shutil.copytree(os.path.join(chosen, "seq"),
                                os.path.join(self.__param['dataDir'], "seq"))
                if pieces[3] not in chosen and os.path.isfile(os.path.join(self.__param['dataDir'], "seq", "slices.bin")):
                    os.unlink(os.path.join(
                        self.__param['dataDir'], "seq", "slices.bin"))

    def _is_good_candidate(self, old_param):
        for key in ['database', 'databaseCtime', 'window']:
            try:
                if old_param[key] != self.__param[key]:
                    return False
            except KeyError:
                return False
        return True


def from_command_line(margOut=False, sliceSeparated=False, association={}, **kwargs):
    logging.info("Read command line arguments")

    print("\n       _.%%%%%%%%%%%%%\n"
          "      /- _%%%%%%%%%%%%%\n"
          "     (_ %\|%%% CP2D %%%\n"
          "        %%%$$$$$$$$$$$%\n"
          "          S%S%%%*%%%%S\n"
          "      ,,,,# #,,,,,,,##,,,\n")

    exp = cp2dExperiment(**kwargs)
    exp.run()
    exp.results(margOut=margOut, machine=False,
                PANStyle="strict", sliceSeparated=sliceSeparated, association=association)
