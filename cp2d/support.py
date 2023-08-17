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

import json
import sys
import os
import numpy as np
import shutil
import configparser
from multiprocessing import cpu_count
from collections import MutableMapping
from ctypes import Structure, addressof, byref, c_char, c_char_p, c_int, c_long, c_uint16, c_uint32, c_ushort, c_void_p, cdll, c_ulong, POINTER, cast, c_char_p, c_double
import logging
import time
import pandas as pd
import typing as tp
from glob import glob
el_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

unicpunt = u"\u00B4\u02B9\u02BC\u02C8\u0301\u2018\u2019\u201B\u2032\u2034\u2037"+u"\u00AB\u00BB\u02BA\u030B\u030E\u201C\u201D\u201E\u201F\u2033\u2036\u3003\u301D\u301E"+u"\u00AD\u2010\u2011\u2012\u2013\u2014\u2212\u2015"+u"\u01C3\u2762"+u"\u266F"+u"\u066A\u2052"+u"\u066D\u204E\u2217\u2731\u00D7"+u"\u201A\uFE51\uFF64\u3001" + \
    u"\u00F7\u0338\u2044\u2215"+u"\u0589\u05C3\u2236"+u"\u203D"+u"\u27E6"+u"\u20E5\u2216"+u"\u301B"+u"\u02C4\u02C6\u0302\u2038\u2303"+u"\u02CB\u0300\u2035"+u"\u2983" + \
    u"\u01C0\u05C0\u2223\u2758\u00A6"+u"\u02DC\u0303\u2053\u223C\u301C"+u"\u2039\u2329\u27E8\u3008" + \
    u"\u203A\u232A\u27E9\u3009"+u"\u02CD\u0331\u0332\u2017\u0333" + \
    u"\u200B\u2060\uFEFF\u00a0 "
# the characters u"\u200B\u2060\uFEFF\u00a0 " must be treated differently depending on the method
equivalen = "'"*11+'"'*14+"-"*8+"!"*2+"#"+"%"*2+"*"*5+","*4+"/"*4+":" * \
    3+"?"+"["+"\\"*2+"]"+"^"*5+"`"*3+"{"+"|"*5+"~"*5+"<"*4+">"*4+"_"*5+" "*5
control = ''.join([chr(i) for i in range(int('a', 16))]+[chr(i) for i in range(
    int('e', 16), int('20', 16))]+[chr(i) for i in range(int('7f', 16), int('a0', 16))])

nonAlfa = "!\"#$%&'()*+,-./:;<=>?@[\\]^`{_|}~«»§°1234567890"


class cyclecount:
    _max = 0
    _now = 0

    def __init__(self, m, n=0):
        if m <= 0 or n < 0:
            raise ValueError()
        if isinstance(m, int) and isinstance(n, int):
            self._max = m
            self._now = n
        else:
            raise TypeError()

    @property
    def now(self):
        return self._now

    @now.setter
    def now(self, o):
        if isinstance(o, int):
            self._now = o % self._max
        else:
            raise TypeError()

    def __add__(self, o):
        if isinstance(o, int):
            return cyclecount(self._max, (self.now+o) % self._max)
        else:
            raise TypeError()

    def __mul__(self, o):
        if isinstance(o, int):
            return cyclecount(self._max, (self.now*o) % self._max)
        else:
            raise TypeError()

    __rmul__ = __mul__

    def __eq__(self, o):
        return self.now == o

    def __ge__(self, o):
        return self.now >= o

    def __le__(self, o):
        return self.now <= o

    def __ne__(self, o):
        return self.now != o

    def __lt__(self, o):
        return self.now < o

    def __gt__(self, o):
        return self.now > o

    def __str__(self):
        return str(self.now)

    def __int__(self):
        return self.now

    def __float__(self):
        return float(self.now)

    def __floordiv__(self, o):
        return self.now//o

    def __truediv__(self, o):
        return self.now/o

    def __rtruediv__(self, o):
        return o/self.now

    def __mod__(self, o):
        return self.now % o


class dataLoad (Structure):
    _fields_ = [("fname", c_char_p), ("bucket", c_uint16), ("N", c_long), ("offs", c_long), ("ncids", c_long),
                ("cids", POINTER(c_ushort)), ("goods", POINTER(c_char)), ("fn", c_long), ("lengths", POINTER(c_uint32)), ("F", c_uint32)]


_libRes = cdll.LoadLibrary(os.path.join(el_path, 'bin/libattributor.so'))
_libRes.load.argtypes = (dataLoad, POINTER(POINTER(c_double)))
_libRes.load.restype = c_int
_libRes.attribute.argtypes = (
    POINTER(c_double), c_ulong, c_double, POINTER(POINTER(c_char)))
_libRes.attribute.restype = c_int
_libRes.clean.argtypes = (c_void_p,)


class pointerManager ():
    """This pointer manager is intended to be the only owner of the pointer.
    If you copy the raw pointer around, expect double free or corruption."""

    def __init__(self, ptr=None) -> None:
        self.__valid = False
        self.__set_ptr(ptr)

    def __set_ptr(self, ptr):
        if ptr is None or isinstance(ptr, (POINTER(c_char), POINTER(c_double))):
            self.clean()
            self.__ptr = ptr
            if self.__ptr is not None and self.__ptr != POINTER(c_char)() and self.__ptr != POINTER(c_double)():
                self.__valid = True
            else:
                self.__valid = False
        else:
            self.__valid = False
            raise ValueError("ptr may be only ctypes.POINTER() or None.")

    def clean(self):
        if self.__valid:
            if addressof(self.__ptr.contents) > 10:
                _libRes.clean(cast(self.__ptr, c_void_p))
            self.__ptr = POINTER(c_double)()
            self.__valid = False

    @property
    def valid(self):
        return self.__valid

    @valid.setter
    def valid(self, *args, **kwargs):
        raise ValueError("You should never change pointers validity by hand.")

    @property
    def ptr(self):
        if not self.__valid:
            raise RuntimeWarning("Using invalid pointer.")
        return self.__ptr

    @ptr.setter
    def ptr(self, ptr):
        self.__set_ptr(ptr)

    def __del__(self):
        self.clean()


class resultsExtractor ():
    def __init__(self, proFile, bucket, offs, N, ncids, cids, good, fn, lengths, F) -> None:
        self.data = dataLoad(proFile.encode(), bucket, N, offs, ncids, cids.ctypes.data_as(
            POINTER(c_ushort)), good.ctypes.data_as(POINTER(c_char)), fn, lengths.ctypes.data_as(
            POINTER(c_uint32)), F)
        self.__res = pointerManager()
        self.__proc = pointerManager()
        self.__oldDelta = None
        self.__oldAttri = None

    def __enter__(self):
        return self.load()

    def __exit__(self, exc_type, exc_value, exc_tracebac):
        self.__res.clean()
        self.__proc.clean()

    def __check_result(self, val, msg=""):
        errors = {1: "Cannot open file", 2: "File ended too soon",
                  3: "Unknown I/O error", 4: "Unknown fatal error"}
        if val in errors:
            raise RuntimeError(f"Failed {msg}. {errors[val]}.")
        return True

    def load(self):
        tmp_p = POINTER(c_double)()
        ret = _libRes.load(
            self.data, byref(tmp_p))
        if self.__check_result(ret, "load "+self.data.fname.decode()):
            self.__res.ptr = tmp_p

        return self

    def attribute(self, delta):
        if delta == self.__oldDelta and self.__oldAttri is not None:
            return self.__oldAttri

        self.__oldDelta = delta
        deltac = c_double(delta)
        if not self.__res.valid:
            raise RuntimeError(str(resultsExtractor) +
                               " not initialised (use 'load').")
        tmp_p = POINTER(c_char)()
        res = _libRes.attribute(
            self.__res.ptr, self.data.fn, deltac, byref(tmp_p))
        if self.__check_result(res, "attribution"):
            self.__proc.ptr = tmp_p

        reslen = -res
        self.__oldAttri = np.frombuffer(cast(self.__proc.ptr, POINTER(
            c_char*(reslen*38))).contents, dtype='u2, f8, f8, f8, i4, i4, i4', count=reslen)
        self.__oldAttri.dtype.names = (
            "cid", "FNN", "TOP", "WP", "MR", "TMR", "WMR")
        return self.__oldAttri

    @property
    def raw_results_p(self):
        return self.__res

    @raw_results_p.setter
    def raw_results_p(self, *args, **kwargs):
        raise ValueError("You should never change pointers by hand.")

    @raw_results_p.deleter
    def raw_results_p(self):
        self.__res.clean()

    @property
    def processed_results_p(self):
        return self.__proc

    @processed_results_p.setter
    def processed_results_p(self, *args, **kwargs):
        raise ValueError("You should never change pointers by hand.")

    @processed_results_p.deleter
    def processed_results_p(self):
        self.__proc.clean()

    @property
    def cached(self):
        return self.__oldDelta, self.__oldAttri

    @cached.setter
    def cached(self, *args, **kwargs):
        raise ValueError("You should never set results by hand.")

    @cached.deleter
    def cached(self):
        self.__oldDelta, self.__oldAttri = None, None


class ConflictingAssociation(BaseException):
    def __init__(self, message, base_message=None, *args):

        self.message = message
        self.base_message = base_message
        super(ConflictingAssociation, self).__init__()

    def __str__(self):
        if self.base_message is None:
            return self.message

        return self.message + " '" + str(self.base_message) + "'"


def get_input(msg, dir=False, num=False):
    """Promprts message and reads input.

    There are three working modes:
        -# \b yes/no when the user has to answer a simple question;
        -# \b directory when a path (absolute or relative) to a directory is needed.
        -# \b number when a number is needed (always a float)
    The function takes care to accept only valid input and prompt again if needed.

    Args:
        msg:    The message to be prompted.
        dir:    If \b True enables the directory mode (default \b False)
        num:    If \b True enables the number mode (default \b False)

    Returns:
        ch: [bool] The choice made in \b yes/no mode.
        dirName: [str] A valid path to an existing directory in \b directory mode.
        number: [float] A number if in  \b number mode.

    Raises:
        RuntimeError if bot directory and number mode are activated.
    """
    if dir and num:
        raise RuntimeError(
            "Cant get a directory and a number at the same time.\nOnly one of dir and num can be True.")

    if dir:
        if not sys.stdout.isatty() or not sys.stdin.isatty():
            try:
                ipy_str = str(type(get_ipython()))
                dirName = ""
                while os.path.isdir(dirName):  # stupid way to avoid existing default name
                    dirName += "w"
                while not os.path.isdir(dirName):
                    dirName = input(
                        msg+(" Directory must exists. " if dirName and dirName[0] != 'w' else " "))
                return dirName
            except:
                raise RuntimeError("Input needed but no tty output: "+msg)
        dirName = ""
        while os.path.isdir(dirName):  # stupid way to avoid existing default name
            dirName += "w"
        while not os.path.isdir(dirName):
            dirName = input(
                msg+(" Directory must exists. " if dirName and dirName[0] != 'w' else " "))
        return dirName
    elif num:
        if not sys.stdout.isatty() or not sys.stdin.isatty():
            try:
                ipy_str = str(type(get_ipython()))
                number = np.nan
                while np.isnan(number):
                    try:
                        number = float(input(msg+" "))
                    except ValueError:
                        print("Invalid. Insert a number: ")
                return number
            except:
                raise RuntimeError("Input needed but no tty output")
        number = np.nan
        while np.isnan(number):
            try:
                number = float(input(msg+" "))
            except ValueError:
                print("Invalid. Insert a number: ")
        return number
    else:
        if not sys.stdout.isatty() or not sys.stdin.isatty():
            try:
                ipy_str = str(type(get_ipython()))
                ch = -1
                while ch == -1:
                    tmp = input(msg+" [y/n] ")
                    if tmp.lower() in ["1", "y", "yes"]:
                        ch = True
                    elif tmp.lower() in ["0", "n", "no"]:
                        ch = False
                    else:
                        print("Invalid. Yes or no?")
                return ch
            except:
                return False
        ch = -1
        while ch == -1:
            tmp = input(msg+" [y/n] ")
            if tmp.lower() in ["1", "y", "yes"]:
                ch = True
            elif tmp.lower() in ["0", "n", "no"]:
                ch = False
            else:
                print("Invalid. Yes or no?")
        return ch


try:
    try:
        ipy_str = str(type(get_ipython()))
        if 'zmqshell' in ipy_str:
            from tqdm.notebook import tqdm
        if 'terminal' in ipy_str:
            from tqdm import tqdm
    except:
        if sys.stderr.isatty():
            from tqdm import tqdm
        else:
            def tqdm(iterable=None, desc=None, leave=True, *args, **kwargs):
                if desc:
                    print(desc, file=sys.stderr, end="\n" if leave else "\r")
                if iterable is not None:
                    return iterable.__iter__()
                else:
                    class FakeT ():
                        def __init__(self) -> None:
                            pass
                        def write (self, *args, **kwargs):
                            print(*args, **kwargs)
                        def update (self, *args, **kwargs):
                            pass
                    return FakeT()
except:
    def tqdm(iterable, desc=None, leave=True, *args, **kwargs):
        if desc:
            print(desc, file=sys.stderr, end="\n" if leave else "\r")
            sys.stderr.flush()
        return iterable.__iter__()


def build_config(globalDir, outputDir="", useGlobal=False, **kwargs):
    """
    Builds interactively a config file.

    Giudes the user in the choice of the parameters for the execution. These are technical parameters that do not interfere with the anlisys.
    The procedure tries to avoid all the most stupid choices the user can make (like try to run without using any processor).

    The file generated with this procedure has the standard name config.ini and any other file with the same name in the folder is renamed
    config_old.ini or config_old<progressive number>.ini.

    If no destination folder is passed as argument prompts the user for one.

    Args:
        dir:    The directory where the file will be created. 
    """
    if useGlobal:
        dir = globalDir
    elif outputDir and os.path.isdir(outputDir):
        dir = outputDir
    elif get_input("You want the configuration to be global?"):
        dir = globalDir
    else:
        dir = get_input(
            "Were do you want to create the config file?", True)

    config = configparser.ConfigParser()
    numCpus = cpu_count()
    if not get_input(f"It seems you have {numCpus} CPUs, go on with this number?"):
        numCpus = -1
        while numCpus <= 0:
            numCpus = int(
                get_input("Insert the desired number of CPUs to work with:", num=True))
            if numCpus <= 0:
                print(f"You want at least one CPU working.")
    config["numcomp"] = {"numCpus": str(numCpus)}
    buffer = numCpus+1
    while buffer >= numCpus or buffer < 0:
        buffer = int(
            get_input("How many processors you want to leave always free?", num=True))
        if buffer >= numCpus or buffer < 0:
            print(f"You want at least one of your {numCpus} CPUs working.")
    config["numcomp"]["buffer"] = str(buffer)
    if get_input(f"You want to use the same numbrer of {numCpus-buffer} CPUs in every parallelized section?"):
        config["elaboration"] = {"workers": str(numCpus-buffer)}
    else:
        workers = numCpus+1
        while workers > numCpus or workers <= 0:
            workers = int(
                get_input("How many processors you want to use?", num=True))
            if workers > numCpus or workers <= 0:
                print(f"You want at least one of your {numCpus} CPUs working.")
        config["elaboration"] = {"workers": str(workers)}

    limit = -1
    while limit < 0 or limit > numCpus:
        limit = int(
            get_input("Set the maximum number of CPUs to leave to others:", num=True))
        if limit < 0 or limit > numCpus:
            print(f"Can't stop more CPUs than you have.")
        elif limit == numCpus:
            if get_input("Are you sure you want to stop completely if someone else joins?"):
                break
            else:
                limit = -1
        elif limit < buffer:
            if get_input("A limit smaller than the buffer will disable the job manager, are you sure?"):
                break
            else:
                limit = -1
    config["numcomp"]["limit"] = str(limit)

    timeout = 0.5
    if not get_input(f"I suggest {timeout} s between checks if the main computation has finished, accept?"):
        timeout = -1
        while timeout < 0:
            timeout = get_input("Set the timeout:", num=True)
            if timeout == 0:
                if get_input("A null timeout may break the interactive monitoring ad be resource intensive, are you sure?"):
                    break
                else:
                    timeout = -1
    config["numcomp"]["timeout"] = str(timeout)

    cpuThreshold = 10
    if not get_input(f"I suggest a {cpuThreshold}% CPU threshold to consider a process from another user, accept?"):
        cpuThreshold = -1
        while cpuThreshold < 0 or cpuThreshold > 100:
            cpuThreshold = get_input("Set the cpuThreshold:", num=True)
            if cpuThreshold == 0:
                if get_input("A null cpuThreshold will stop your processes for every shell opened by other users, are you sure?"):
                    break
                else:
                    cpuThreshold = -1
    config["numcomp"]["cpuThreshold"] = str(cpuThreshold)

    Encoding = "utf-8"
    if get_input("Do you know the encoding of ypur input files?"):
        again = True
        while again:
            Encoding = input("Insert a valid encoding name: ")
            try:
                b'qq'.decode(Encoding)
            except LookupError:
                continue
            again = False
    else:
        print("Assuming 'utf-8'.")
    config["elaboration"]["Encoding"] = Encoding

    maxPreloadResultsMB = 750
    print("With big corpora and small fragments result files may grow large (>1GB).\nLoading results a bit at the time is slower but uses less memory.")
    if not get_input(f"I suggest to switch mode for files larger than  {maxPreloadResultsMB} MB, accept?\n\t(Increase [decrease] if your RAM is huge [tiny])"):
        maxPreloadResultsMB = -1
        while maxPreloadResultsMB < 0:
            maxPreloadResultsMB = get_input(
                "Set the maxPreloadResultsMB:", num=True)
            if maxPreloadResultsMB == 0:
                if get_input("A null maxPreloadResultsMB will disable preloading altogether, are you sure?"):
                    break
                else:
                    maxPreloadResultsMB = -1
    config["elaboration"]["maxPreloadResultsMB"] = str(maxPreloadResultsMB)

    if os.path.isfile(os.path.join(dir, "config.ini")):
        if not os.path.isfile(os.path.join(dir, "config_old.ini")):
            shutil.move(os.path.join(dir, "config.ini"),
                        os.path.join(dir, "config_old.ini"))
            print("A config file was already there, renamed config_old.ini.")
        else:
            i = 0
            while os.path.isfile(os.path.join(dir, f"config_old{i}.ini")):
                i += 1
            shutil.move(os.path.join(dir, "config.ini"),
                        os.path.join(dir, f"config_old{i}.ini"))
            print(
                f"A config file was already there, renamed config_old{i}.ini.")
    with open(os.path.join(dir, "config.ini"), "w") as configfile:
        config.write(configfile)
    print("Configuration complete.")
    return os.path.join(dir, "config.ini")


def data_preprocessing(inFolder="", outFolder="", authorFiles: bool = False, sep: str = '_', bookStart: int = 2, extension: str = "txt", **kwargs):
    """
    Prepares raw .txt books for the analisys.

    The format of the filenames must be <AuthorName>_<year>_<title>.txt where spaces in the author's name must be suppressed, year can be anything but must not contain underscores and spaces in the book title must be replaced with underscores.

    Books are renamed and copied in the folder chosen by the user, the new names follow the format A<author id>B<book id>.txt. In the folder are created also files containing the original name of the authors ("authorNames.dat") and of the books ("bookNames.dat") for future reference.

    This is the only function that may be accessed by calling directly the module from commandline with subcommand dataprep:

        ./elaboration.py dataprep
    """
    if not os.path.isdir(inFolder):
        inFolder = ""
    if not os.path.isdir(outFolder):
        outFolder = ""
    if not (outFolder and inFolder):
        print("Folders must be valid paths to existing folders.")
        if not inFolder:
            inFolder = get_input("Insert input folder: ", True)
        if not outFolder:
            outFolder = get_input("Insert output folder: ", True)

    if os.path.isfile(os.path.join(inFolder, "encoding.dat")):
        shutil.copy(os.path.join(inFolder, "encoding.dat"),
                    os.path.join(outFolder, "encoding.dat"))
        if authorFiles:
            with open(os.path.join(inFolder, "encoding.dat")) as fp:
                encoding = fp.read()
            for file in os.listdir(outFolder):
                if file.endswith("txt"):
                    os.unlink(os.path.join(outFolder, file))

    elif authorFiles:
        raise FileNotFoundError(os.path.join(
            inFolder, "encoding.dat")+" is missing, can't create author files.")
    authors = []
    books = {}
    for name in sorted(os.listdir(inFolder)):
        if name[0] != '.' and name.endswith(extension):
            baseName = os.path.splitext(name)[0]
            parts = baseName.split(sep)
            if not parts[0] in authors:
                authors.append(parts[0])
            try:
                books[authors.index(parts[0])].append(
                    sep.join(parts[bookStart:]))
            except:
                books[authors.index(parts[0])] = [
                    sep.join(parts[bookStart:]), ]
            if authorFiles:
                with open(os.path.join(outFolder, f"A{authors.index(parts[0])+1}.txt"), "a", encoding=encoding) as fout, open(os.path.join(inFolder, name), encoding=encoding) as fin:
                    text_b = fin.read().replace("\n", " ")
                    if text_b[0] == '#':
                        text_b = "".join(["_added_", text_b])
                    print(
                        f"# {len(books[authors.index(parts[0])])}", text_b, sep="\n", file=fout)
            else:
                newName = f"A{authors.index(parts[0])+1}B{len(books[authors.index(parts[0])])}.txt"
                shutil.copy(os.path.join(inFolder, name),
                            os.path.join(outFolder, newName))
    with open(os.path.join(outFolder, "authorNames.dat"), "w") as fp:
        for n, auth in enumerate(authors, start=1):
            print(n, auth, file=fp)

    with open(os.path.join(outFolder, "bookNames.dat"), "w") as fp:
        for A in books:
            for n, book in enumerate(books[A], start=1):
                print(A+1, n, book, file=fp)


def corpus_shuffler(inFolder="", preserveDist=None, **kwargs):
    if not os.path.isdir(inFolder):
        print("Folder's path must be a valid path to an existing and prepared corpus.")
        inFolder = get_input("Insert input folder: ", True)
    if preserveDist is None:
        preserveDist = get_input(
            "Do you want to preserve the distribution of authors' size? ")
    books = [a.split(".")[0]
             for a in os.listdir(inFolder) if a.endswith("txt")]
    shelf = {}
    for book in books:
        aut = book.split("B")[0]
        try:
            shelf[aut].append(book)
        except:
            shelf[aut] = [book, ]
    order = [a for a in sorted(
        shelf, reverse=True, key=lambda x: len(shelf[x]))]
    if preserveDist:
        shelfsize = [(i, len(shelf[a])) for i, a in enumerate(shelf)]
        shelfsize.sort(reverse=True, key=lambda x: x[1])
    else:
        tot = sum([len(shelf[a]) for a in shelf])
        shelfsize = [(i, tot//len(shelf)+(1 if i < tot % len(shelf) else 0))
                     for i in range(len(shelf))]
    newshelf = [[] for _ in shelf]
    autN = cyclecount(len(shelf))
    for aut in order:
        for book in shelf[aut]:
            while len(newshelf[int(autN)]) >= shelfsize[int(autN)][1]:
                autN += 1
            newshelf[int(autN)].append(book)
            autN += 1
    outdir = os.path.dirname(os.path.join(inFolder, ""))+"-shuf0"
    i = 1
    while os.path.exists(outdir):
        outdir = outdir[:-1]+str(i)
        i += 1
    os.mkdir(outdir)
    print("Created: ", outdir)
    for n, autlist in enumerate(newshelf, start=1):
        for b, tit in enumerate(autlist, start=1):
            shutil.copyfile(os.path.join(
                inFolder, f"{tit}.txt"), os.path.join(outdir, f"A{n}B{b}.txt"))
    if os.path.isfile(os.path.join(inFolder, f"encoding.dat")):
        shutil.copyfile(os.path.join(inFolder, f"encoding.dat"),
                        os.path.join(outdir, f"encoding.dat"))


class lockable_dict(dict, MutableMapping):
    def __new__(cls, *args, **kwargs):
        obj = super().__new__(cls)
        obj.__lock = False
        obj.__checks = {}
        return obj

    def __init__(self, dict={}, locked: bool = False, checks={}):
        super(lockable_dict, self).__init__(dict)
        self.__lock = locked
        if checks:
            self.register_checks()

    def __getitem__(self, key):
        return dict.__getitem__(self, key)

    def __setitem__(self, key, value):
        if self.__lock:
            raise ValueError(
                "Can't change parameters when the experiment has already begun.")
        if key in self.__checks:
            if self.__checks[key](value):
                dict.__setitem__(self, key, value)
            else:
                raise ValueError(f"Value {value} for element {key} invalid.")
        else:
            dict.__setitem__(self, key, value)

    def __delitem__(self, key):
        dict.__delitem__(self, key)

    def __iter__(self):
        return dict.__iter__(self)

    def __len__(self):
        return dict.__len__(self)

    def __contains__(self, x):
        return dict.__contains__(self, x)

    def lock(self):
        if not self.__lock:
            self.__lock = True
        else:
            raise UserWarning("{v} already locked.".format(
                v=self.__class__.__name__))

    def unlock(self):
        if self.__lock:
            self.__lock = False
        else:
            raise UserWarning("{v} already unlocked.".format(
                v=self.__class__.__name__))

    def islocked(self):
        return self.__lock

    def register_checks(self, checkDict: dict = None, **kwargs):
        if checkDict and type(checkDict) == dict:
            for key, check in checkDict.items():
                if callable(check):
                    self.__checks[key] = check
                else:
                    raise ValueError(
                        f"Checks must be callable. Check the check for {key}.")
        for key, check in kwargs.items():
            if callable(check):
                self.__checks[key] = check
            else:
                raise ValueError(
                    f"Checks must be callable. Check the check for {key}.")

    def run_checks(self):
        for key, value in self.items():
            if key in self.__checks:
                if not self.__checks[key](value):
                    raise ValueError(
                        f"Value \"{value}\" of type {type(value)} for element {key} invalid.")

    def __getstate__(self):
        tmpdict = {k: v for k, v in self.__dict__.items()}
        tmpdict['_lockable_dict__checks'] = {}
        return tmpdict

    def __setstate__(self, state):
        self.__dict__.update(state)


def numpy_json(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    else:
        raise TypeError


def auth_books_iterator(path, encoding, auth_file: bool):
    if auth_file:
        with open(path+".txt", encoding=encoding) as f_in:
            bn = -1
            known = set()
            for line in f_in:
                if line[0] == '#':
                    if bn > 0:
                        yield bn, ""
                    bn = int(line[1:].strip())
                    if bn in known:
                        raise FileExistsError(
                            f"Already found a document {bn} in author {path}.")
                    known.add(bn)
                elif line.strip():
                    if bn < 0:
                        raise ValueError("Wrong author file format.")
                    if line.startswith("_added_#"):
                        line = line[7:]
                    yield bn, line.strip()
                    bn = -1
    else:
        for fileName in glob(path+"B*.txt"):
            bn = os.path.basename(fileName).split("B")[1].split(".")[0]
            with open(fileName, encoding=encoding) as f_in:
                text_b = f_in.read()
            yield bn, text_b


class fakeGet:
    def __init__(self, n: int = 0) -> None:
        self.__n = n

    def __getitem__(self, __: tp.Any, *args, **kwargs) -> int:
        return self.__n


def progress_wait(extraction_jobs):
    logger = logging.getLogger(__name__)
    remaining_jobs = [j for j in extraction_jobs if not j.ready()]
    switch_log = 0
    print(f"\r                                          ", end="")
    done = -1
    while remaining_jobs:
        new_done = len(extraction_jobs)-len(remaining_jobs)
        if new_done != done:
            switch_log = (switch_log+1) % 5
            done = new_done
            if not switch_log:
                logger.info(
                    f"Extracted {done} files")
            else:
                logger.debug(
                    f"Extracted {done} files")
            print(f"\r{(1-len(remaining_jobs)/len(extraction_jobs))*100:6.6}%  \t({done:{len(str(len(extraction_jobs)+1))}})", end="", flush=True)
        time.sleep(1)
        remaining_jobs = [j for j in remaining_jobs if not j.ready()]


def load_book_len(dataDir, Encoding="latin1", output="var"):
    if output not in ["file", "both", "var"]:
        raise ValueError("output must be one of file, both, var")

    if not os.path.isdir(os.path.join(dataDir, "wnt")):
        raise FileNotFoundError(
            "The wnt directory is missing. Maybe not created or already cleaned.")

    booktok = {}
    for file in os.listdir(os.path.join(dataDir, "wnt")):
        A, B = file[1:].split(".")[0].split("B")
        N = 0
        with open(os.path.join(dataDir, "wnt", file), encoding=Encoding) as fp:
            for line in fp:
                if line[0] == "#":
                    continue
                __, c = line.split()
                N += int(c)
        booktok[(int(A), int(B))] = N

    if output in ["file", "both"]:
        with open(os.path.join(dataDir, "bookLen.dat"), "w") as fp:
            for ab in booktok:
                print(*ab, booktok[ab], sep="\t", file=fp)

    if output in ["var", "both"]:
        return booktok


class associator:
    def __init__(self, asso: dict = None) -> None:
        assert isinstance(asso, dict)
        self.__asso = asso
        if asso:
            def getitem(self, it):
                try:
                    return self.__asso[it]
                except:
                    return it
        else:
            def getitem(self, it):
                return it
        self.getitem = getitem

    def __getitem__(self, it):
        if hasattr(it, '__iter__'):
            return np.array([self.getitem(self, i) for i in it])
        return self.getitem(self, it)


def PAN11_stats(attribution, mode="strict", groundTruth=None, reject=None, shallow=False, sliceSeparated: tp.Optional[int] = None, slices: tp.Optional[tp.Mapping[tp.Sequence[int], int]] = None, authors: tp.Sequence[int] = None):
    """Attribution is a dictionary as returned from cp2d.elaboration.cp2dExperiment.attributions.
        Mode is one of strict, optimist, partial describing what to do when more than one author
        is proposed for the same text.
        groundTruth is a dictionary whose keys are the numbers of the books with unknown author
        and the values are the number of the true authors.
        reject is a callable that receives a tuple with (author number, book number) and the
        attribution info with number of the first and second classified authors and their score."""
    if not reject:
        def reject(ab, r): return False
    if groundTruth:
        if isinstance(groundTruth, str):
            with open(groundTruth) as fp:
                groTru = {int(k): v for k, v in json.load(fp).items()}
        elif isinstance(groundTruth, dict):
            groTru = groundTruth
        else:
            raise TypeError(
                "groundTruth should be a dictionary(int->int) or a (relative)path to a json containing such dictionary.")
    if authors is None:
        authors = set((ab[0] for ab in attribution))

    if groundTruth or not sliceSeparated:
        slices = fakeGet()
        sliceSeparated = 1

    KA = np.array([(a > 0) for a in authors], dtype=bool)  # known author
    if groundTruth:
        KA = ~KA
    AO = {a: i for i, a in enumerate(authors)}  # author order
    TB = np.ma.zeros((len(authors), sliceSeparated))  # author books in slice
    for ab in attribution:
        if groundTruth and not ab[0]:
            tau = groTru[ab[1]]
        elif not groundTruth and ab[0]:
            tau = ab[0]
        else:
            continue
        if tau < 0:
            continue
        TB[AO[tau], slices[ab]] += 1
    if mode not in ["strict", "optimist", "partial"]:
        raise ValueError("Mode must be one of strict, optimist, partial.")
    # for a in authors:
    #     if a not in stats:
    #         stats[a] = [{"co": {"FNN": 0, "TOP": 0, "WP": 0, "MR": 0, "TMR": 0, "WMR": 0},
    #                     "at": {"FNN": 0, "TOP": 0, "WP": 0, "MR": 0, "TMR": 0, "WMR": 0}} for __ in range(sliceSeparated)]
    for ab in attribution:
        sab = ab
        break
    micro = [{} for __ in range(sliceSeparated)]
    macro = [{} for __ in range(sliceSeparated)]
    for delta in [1] if shallow else attribution[sab]:
        stats = [{a: {"co": {"FNN": 0, "TOP": 0, "WP": 0, "MR": 0, "TMR": 0, "WMR": 0},
                      "at": {"FNN": 0, "TOP": 0, "WP": 0, "MR": 0, "TMR": 0, "WMR": 0}} for a in authors} for __ in range(sliceSeparated)]
        for ab in attribution:
            if groundTruth and not ab[0]:
                tau = groTru[ab[1]]
            elif not groundTruth and ab[0]:
                tau = ab[0]
            else:
                continue
            if tau < 0:
                continue
            for t, r in attribution[ab].items() if shallow else attribution[ab][delta].items():
                if t != "FRA":
                    try:
                        winners = r[0][0]
                        if winners == [-1] or reject(ab, r):
                            continue
                        if tau in winners:
                            if mode == "optimist" or len(winners) == 1:
                                stats[slices[ab]][tau]["co"][t] += 1
                            elif mode == "partial":
                                stats[slices[ab]][tau]["co"][t] += 1/len(winners)
                        if mode == "partial":
                            for A in winners:
                                stats[slices[ab]][A]["at"][t] += 1/len(winners)
                        else:
                            for A in winners:
                                stats[slices[ab]][A]["at"][t] += 1
                    except:
                        print(t, r)
                        print(attribution[ab])
                        raise
        PRaF = [{} for __ in range(sliceSeparated)]
        for i, s in enumerate(stats):
            for a in s:
                if not TB[AO[a], i]:
                    continue
                PRaF[i][a] = {"P": {}, "R": {}, "F": {}, "a": TB[AO[a], i]}
                for t in s[a]['co']:
                    if not s[a]['co'][t]:
                        PRaF[i][a]['P'][t] = 0
                        PRaF[i][a]['F'][t] = 0
                        PRaF[i][a]['R'][t] = 0
                    else:
                        PRaF[i][a]['P'][t] = s[a]['co'][t]/s[a]['at'][t]
                        PRaF[i][a]['R'][t] = s[a]['co'][t]/TB[AO[a], i]
                        PRaF[i][a]['F'][t] = 2*s[a]['co'][t] / \
                            (TB[AO[a], i]+s[a]['at'][t])

        for i, prf in enumerate(PRaF):
            micro[i][delta] = {}
            macro[i][delta] = {}
            for S in ['P', 'F', 'R']:
                micro[i][delta][S] = {}
                macro[i][delta][S] = {}
                for t in ["FNN", "TOP", "WP", "MR", "TMR", "WMR"]:
                    micro[i][delta][S][t] = np.average([prf[a][S][t] for a in prf], weights=[
                        prf[a]['a'] for a in prf])
                    macro[i][delta][S][t] = np.average(
                        [prf[a][S][t] for a in prf])
    if shallow:
        micro = [micro[i][1] for i in range(sliceSeparated)]
        macro = [macro[i][1] for i in range(sliceSeparated)]
    if sliceSeparated == 1:
        micro = micro[0]
        macro = macro[0]
    return micro, macro


def print_PAN11(_micro, _macro, results=None, shallow=False, sliceSeparated: tp.Optional[int] = None):
    """micro and macro averages results are those produced by PAN11_stats.
        results is a (path to a) pandas DataFrame containing the results obtained in the
        PAN11 competition for comparison. The seven columns of interest are: RankSum and
        Prec, Recall and F1 prefixed with m_ or M_ for micro and macro averages."""
    if results:
        if isinstance(results, str):
            with open(results) as fp:
                PANres = pd.read_csv(results)
        elif isinstance(results, pd.DataFrame):
            PANres = results
        else:
            raise TypeError(
                "groundTruth should be a DataFrame or a (relative)path to a csv containing such DataFrame.")
        microPAN = {}
        macroPAN = {}
        colnames = {'P': "Prec", 'R': "Recall", 'F': "F1"}
        for S in colnames:
            microPAN[S] = PANres[f"m_{colnames[S]}"].sort_values().to_numpy()
            macroPAN[S] = PANres[f"M_{colnames[S]}"].sort_values().to_numpy()
        PANtop = PANres.RankSum.min()

    for micro, macro in zip(_micro, _macro) if sliceSeparated else [(_micro, _macro)]:
        for delta in [1] if shallow else micro:
            print(('' if shallow else f'Delta: {10**delta:3}\n') +
                  f"\t\tMicro\t\t\tMacro\n\tP\tR\tF1\tP\tR\tF1\t{''if not results else f'(Min Rank:{PANtop})'}")
            tmic = micro if shallow else micro[delta]
            tmac = macro if shallow else macro[delta]
            for t in tmic['P']:
                print(t, *(round(tmic[S][t], 4) for S in ['P', 'R', 'F']),
                      *(round(tmac[S][t], 4) for S in ['P', 'R', 'F']), sep="\t", end="\t")
                if results:
                    rankSum = 0
                    for S in tmic:
                        rankSum += microPAN[S].size+1 - \
                            np.searchsorted(microPAN[S], tmic[S][t])
                        rankSum += macroPAN[S].size+1 - \
                            np.searchsorted(macroPAN[S], tmac[S][t])
                    print(rankSum)
                else:
                    print()


def PAN11(attributions, PANStyle="strict", groundTruth=None, othResult=None, sliceSeparated: tp.Optional[int] = None, slices: tp.Optional[tp.Mapping[tp.Sequence[int], int]] = None, authors: tp.Sequence[int] = None, shallow=True):
    if groundTruth:
        micro, macro = PAN11_stats(
            attributions, PANStyle, groundTruth, shallow=shallow)
        print("Unknown:")
    else:
        micro, macro = PAN11_stats(attributions, PANStyle, shallow=shallow,
                                   sliceSeparated=sliceSeparated, slices=slices, authors=authors)
        print("Known:")
    print_PAN11(micro, macro, othResult, shallow=shallow,
                sliceSeparated=sliceSeparated)
    return micro, macro
