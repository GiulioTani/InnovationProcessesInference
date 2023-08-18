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

from . import support
import os
from ctypes import c_char_p, c_int, cdll, create_string_buffer, c_char_p
import logging
import numpy as np
import sys
import configparser
import json
import pandas as pd
import typing as tp
from numpy.lib.recfunctions import structured_to_unstructured
logger = logging.getLogger(__name__)
el_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

libLZ = cdll.LoadLibrary(os.path.join(el_path, 'bin/libLZ77dict.so'))


def LZ77dict(text, window, Encoding):
    # strict upper limit is 4/3 aka 12/9, this accounts for rounding and stray bytes
    s = create_string_buffer(b'\000' * int(len(text)*14/10))
    if libLZ.LZ77dict(c_char_p(text.encode(Encoding)), c_int(window), s):
        raise RuntimeError("Error in compression")
    return s.value.decode(Encoding)  # text,c_int(len(text)),len(text)


def output_file_names(dataDir, filename, extension):
    outputFileName = os.path.join(
        dataDir, extension, f"{filename}."+extension)
    outputFileTmpName = os.path.join(
        dataDir, extension, f"tmp{filename}."+extension)
    return outputFileName, outputFileTmpName


def write_sequence(bookNumber, sequence, file) -> None:
    print(f"# {bookNumber}", file=file)
    print(" ".join(sequence), file=file)


def get_association_from_file(dataDir, association):
    if os.path.isfile(os.path.join(dataDir, association)):
        with open(os.path.join(dataDir, association)) as fp:
            tmp_association = json.load(fp)
        association = {int(a): b for a, b in tmp_association.items()}
        return association
    else:
        raise FileNotFoundError(
            f"Missing association file '{os.path.join(dataDir, association)}'.")


def apply_association(bookSet, association):
    # since associations are authorwise i have to create a mapping for books given associations and books found
    booksAssociation = {}
    if association:
        autass = {a: a for a in set(ab[0] for ab in bookSet) if a > 0}
        autass.update(association)
        for a in set(autass.values()):
            mali = [b for b in autass if autass[b] == a]
            if len(mali) == 1 and mali[0] == a:
                continue
            if a in mali:
                next = max(ab[1] for ab in bookSet if ab[0] == a)+1
                mali.remove(a)
            else:
                next = 1
            for a2 in sorted(mali):
                for b in sorted(ab[1] for ab in bookSet if ab[0] == a2):
                    booksAssociation[(a2, b)] = (a, next)
                    next += 1
    return booksAssociation


def fragments_numbers(ids: tp.Union[np.array, tp.Sequence[np.array]], booksAssociation: dict, excluded: list, fra: pd.DataFrame):
    bookAss = {}
    bookAss = {tuple(ab): tuple(ab) for ab in fra.index.tolist()}
    bookAss.update(booksAssociation)
    partsList = []
    offset = 0

    if isinstance(ids, list):
        for n, buck in support.tqdm(enumerate(ids), total=len(ids), desc="Load file", dynamic_ncols=True, disable=None, file=sys.stderr, leave=False):
            offset = 0
            for i in buck:
                b2 = (np.uint32(i['a']), np.uint32(i['b']))
                b2 = bookAss[b2]
                if b2[0] not in excluded:
                    partsList.append((b2, (offset, i['tb'], n)))
                offset += i['tb']
    else:
        for i in support.tqdm(ids, desc="Load fragments", dynamic_ncols=True, disable=None, file=sys.stderr, leave=False):
            b2 = (np.uint32(i['a']), np.uint32(i['b']))
            b2 = bookAss[b2]
            if b2[0] not in excluded:
                partsList.append((b2, (offset, i['tb'], 9999)))
            offset += i['tb']

    return partsList


def loadSlices(path, __slice_t, __auth_t, __book_t):
    global _Q_which_slice
    _Q_which_slice = {}
    if os.path.isfile(path):
        with open(path, "rb") as fin:
            N_slice = np.fromfile(fin, dtype=__slice_t, count=1)[
                0]  # number of slices
            for i in range(N_slice):
                N_auth = np.fromfile(fin, dtype=__auth_t, count=1)[
                    0]  # number of authors in slice
                for __ in range(N_auth):
                    # author id, number of books
                    AuNb = np.fromfile(
                        fin, dtype=[('Au', __auth_t), ('Nb', __book_t)], count=1)
                    _Q_which_slice.update({(AuNb['Au'][0], a): i for a in np.fromfile(
                        fin, dtype=__book_t, count=AuNb['Nb'][0])})
        return N_slice
    else:
        raise FileNotFoundError(f"Missing file slices.bin")

def selectDelta (ret):
    oldres_df=pd.DataFrame([{"delta":d, "FNN":v["all"]["FNN"], "MR":v["all"]["MR"]} for d,v in ret.items()])
    maxFNN = oldres_df.FNN.max()
    maxMR = oldres_df.MR.max()
    d_maxFNN = set(oldres_df[oldres_df.FNN==maxFNN].delta.unique())
    d_maxMR = set(oldres_df[oldres_df.MR==maxMR].delta.unique())
    d_max = d_maxFNN & d_maxMR
    if len(d_max)==0:
        d_max = d_maxFNN if maxFNN >= maxMR else d_maxMR
    return sorted(list(d_max),key=lambda x: np.abs(x))

def preproc_pro(ab):
    global _Q_association, _Q_excluded, _Q_which_slice, _Q_apt, _Q_fra, _Q_which_slice_is_None
    if _Q_association:
        def goodness(a): return (a not in _Q_association and a not in _Q_excluded) or (
            a in _Q_association and _Q_association[a] not in _Q_excluded)
    else:
        def goodness(a): return a not in _Q_excluded
    if not _Q_which_slice_is_None:
        s = _Q_which_slice[ab]
        lcid = {a: (_Q_apt[a][0][s] if ab[0] != 0 and s
                    in _Q_apt[a][0] else _Q_apt[a][1]) for a in _Q_apt}
        cids = np.zeros(sum(lcid.values()), 'u2')
        tcid = 0
        if _Q_association:
            def mycid(
                a): return a if a not in _Q_association else _Q_association[a]
        else:
            def mycid(a): return a
        for a in lcid:
            cids[tcid:tcid+lcid[a]] = mycid(a)
            tcid += lcid[a]

        if _Q_excluded:
            good = np.zeros(cids.size, dtype='b')
            tcid = 0
            for a in lcid:
                good[tcid:tcid+lcid[a]] = goodness(a)
                tcid += lcid[a]
        else:
            good = np.ones_like(cids, dtype='b')
    else:
        cids = np.array(list(_Q_apt.keys()), 'u2')
        if _Q_excluded:
            good = np.array([goodness(a) for a in _Q_apt], 'b')
        else:
            good = np.ones_like(cids, dtype='b')
    ncids = len(cids)
    if 'fn' in _Q_fra.columns:
        fn = _Q_fra.loc[ab, "fn"]
        lengths = _Q_fra.loc[ab, "f_siz"]
    else:
        fn = 1
        lengths = np.array([_Q_fra.loc[ab, "f_siz"]])

    return ncids, cids, good, fn, lengths


def attributor(bookList):
    global _Q_association, _Q_proFile, _Q_missingDelta, _Q_topRes, _Q_F
    # print(globals())
    ab = bookList[0]
    parts = bookList[1]
    offs, N, buck = parts
    ncids, cids, good, fn, lengths = preproc_pro(ab)
    truecid = support.associator(_Q_association)
    with support.resultsExtractor(_Q_proFile, buck, offs, N, ncids, cids, good, fn, lengths, _Q_F) as resGen:
        ret = {}
        for nowDelta in _Q_missingDelta:
            try:
                scores = resGen.attribute(nowDelta)
                scores['cid'] = truecid[scores['cid'].astype(int)]
                tmp = {}
                for typ in ["FNN", "TOP", "WP"]:
                    tarr = scores[scores[typ].argsort()][["cid", typ]]
                    tmp[typ] = [[[truecid[int(q["cid"])]], q[typ]]
                                for q in tarr[slice(None, _Q_topRes, -1)]]
                for typ in ["MR", "TMR", "WMR"]:
                    tmp[typ] = top_performers(
                        scores[scores[typ].argsort()][["cid", typ]], truecid)

                if ab[0] > 0:
                    tarr = scores[np.where(scores["cid"] == ab[0])]
                    tmp["FRA"] = [np.sum(tarr[typ])
                                  for typ in ["MR", "TMR", "WMR"]]
                    if not tmp["FRA"]:
                        tmp["FRA"] = [0, 0, 0]
                else:
                    tmp["FRA"] = [0, 0, 0]
            except Exception as e:
                print(f"{ab}: something strange. . .", e, file=sys.stderr)
                dumb = [([-10, ], -np.inf), ([-10, ], -np.inf)]
                dumbi = [([-10, ], -1), ([-10, ], -1)]
                tmp = {"FNN": dumb, "TOP": dumb, "WP": dumb, "MR": dumbi,
                       "TMR": dumbi, "WMR": dumbi, "FRA": [-1, -1, -1]}
            ret[nowDelta] = tmp
    return ab, ret


def assign(attributions: tp.Mapping[tp.Sequence[int], tp.Mapping[float, tp.Mapping[str, tp.Sequence]]]) -> pd.DataFrame:
    """From the proposed attributions compute base statistics.

    Args:
        attributions (tp.Mapping[tp.Sequence[int,int],tp.Mapping[float,tp.Mapping[str,tp.Sequence]]]): Attribution object, a map with (author,book) tuples as keys and the result of cp2dExperiment.attributor as values.

    Returns:
        pd.DataFrame: The res object (a map of maps indexed by delta and then by kind of attribution), and the optional per slice res object (not implemented).
    """
    global _Q_which_slice
    ret = {}
    if _Q_which_slice is not None:
        slices = _Q_which_slice
    else:
        slices = support.fakeGet()

    for tab in attributions:
        break   
    for nowDelta in attributions[tab]:
        ret[nowDelta]=[]
        for ab in attributions:
            toppers = {"A":ab[0],"B":ab[1], "S":slices[ab]}
            for attributionType in attributions[tab][nowDelta]:
                if attributionType == "FRA":
                    toppers.update({"FRA_"+k:attributions[ab][nowDelta]["FRA"][i] for i, k in enumerate(["MR", "TMR", "WMR"])})
                else:
                    toppers[attributionType]=attributions[ab][nowDelta][attributionType][0][0]
            ret[nowDelta].append(toppers)
    ret_df = {delta:pd.DataFrame(ret[delta]).set_index(["A","B","S"]).sort_index() for delta in ret}
    return pd.concat(ret_df, axis=1)

def return_dict(attributions: tp.Mapping[tp.Sequence[int], tp.Mapping[float, tp.Mapping[str, tp.Sequence]]], nonAttri: tp.Sequence[int], allowPartial: bool, authors: tp.Sequence, sliceSeparated: tp.Optional[int] = None) -> tp.Tuple[tp.Mapping[float, tp.Mapping[str, tp.Mapping[str, float]]], tp.Optional[tp.List[tp.Mapping[float, tp.Mapping[str, tp.Mapping[str, float]]]]]]:
    """From the proposed attributions compute base statistics.

    Args:
        attributions (tp.Mapping[tp.Sequence[int,int],tp.Mapping[float,tp.Mapping[str,tp.Sequence]]]): Attribution object, a map with (author,book) tuples as keys and the result of cp2dExperiment.attributor as values.
        nonAttri (tp.Sequence[int]): The authors whose book are not to be attributed.
        totalFragments (int): Total number of fragments.
        allowPartial (bool): Allow partial attributions i MR (half attribution if two authors with the same number of fragments and so on)
        nAfn (int): Total number of fragments after excluding the due authors.
        nAbn (int): Total number of books after excluding the due authors.
        sliceSeparated (tp.Optional[int]): If evaluates true must contain the number of slices.
        slices (tp.Optional[tp.Mapping[tp.Sequence[int,int],int]], optional): Association between (author,book) tuples and the corresponding slice number if per-slice results are desired. Defaults to None.

    Returns:
        tp.Tuple[tp.Mapping[float,tp.Mapping[str,tp.Mapping[str,float]]], tp.Optional[tp.List[tp.Mapping[float,tp.Mapping[str,tp.Mapping[str,float]]]]]]: The res object (a map of maps indexed by delta and then by kind of attribution), and the optional per slice res object (not implemented).
    """
    global _Q_fra, _Q_which_slice
    ret = {}
    if sliceSeparated:
        slSep = {}
        slices = _Q_which_slice
    else:
        slSep = None
        slices = support.fakeGet()
        sliceSeparated = 1
    if 'fn' in _Q_fra.columns:
        fra = _Q_fra
    else:
        class fakefra:
            def __init__(self) -> None:
                pass
            loc = support.fakeGet(1)
        fra = fakefra()
    if allowPartial:
        fractyp = float
    else:
        fractyp = int
    if len(attributions)==0:
        raise RuntimeError("Empty attributions in return_dict")
    for tab in attributions:
        break
    if nonAttri:
        notExcludedAuthor = np.array([(a not in nonAttri) and (a > 0)
                      for a in authors], dtype=bool)  # not attributed
    knownAuthor = np.array([(a > 0)
                  for a in authors], dtype=bool)  # known author
    if knownAuthor.all():
        print("Useless array!!!!")
        
    authorOrder = {a: i for i, a in enumerate(authors)}  # author order
    authorBooksInSlice = np.ma.zeros((len(authors), sliceSeparated))  # author books in slice
    authorFragmentsInSlice = np.ma.zeros((len(authors), sliceSeparated)
                     )  # author fragments in slice
    for ab in attributions:
        if ab[0]:
            authorBooksInSlice[authorOrder[ab[0]], slices[ab]] += 1
            authorFragmentsInSlice[authorOrder[ab[0]], slices[ab]] += fra.loc[ab, 'fn']
    authorBooksTotal = np.sum(authorBooksInSlice, 1)
    authorFragmentsTotal = np.sum(authorFragmentsInSlice, 1)
    for nowDelta in attributions[tab]:
        authorCorrectForSliceType = {typ: np.ma.zeros((len(authors), sliceSeparated))
               for typ in attributions[tab][nowDelta]}
        for ab in attributions:
            if not ab[0]:
                continue
            for attributionType in attributions[tab][nowDelta]:
                if attributionType == "FRA":
                    authorCorrectForSliceType["FRA"][authorOrder[ab[0]], slices[ab]
                               ] += attributions[ab][nowDelta]["FRA"][0]
                else:
                    if ab[0] in attributions[ab][nowDelta][attributionType][0][0]:
                        authorCorrectForSliceType[attributionType][authorOrder[ab[0]], slices[ab]
                                 ] += fractyp(1/len(attributions[ab][nowDelta][attributionType][0][0]))
        ret[nowDelta] = {"all": {}}
        if nonAttri:
            ret[nowDelta]["exc"] = {}
        for attributionType in authorCorrectForSliceType:
            tmpTyp = np.sum(authorCorrectForSliceType[attributionType], 1)
            if attributionType == 'FRA':
                ret[nowDelta]['all']['FRA'] = np.sum(
                    tmpTyp[knownAuthor])/np.sum(authorFragmentsTotal[knownAuthor])
            else:
                ret[nowDelta]['all'][attributionType] = np.sum(tmpTyp[knownAuthor])/np.sum(authorBooksTotal[knownAuthor])
            if nonAttri:
                if attributionType == 'FRA':
                    ret[nowDelta]["exc"]['FRA'] = np.sum(
                        tmpTyp[notExcludedAuthor])/np.sum(authorFragmentsTotal[notExcludedAuthor])
                else:
                    ret[nowDelta]['exc'][attributionType] = np.sum(
                        tmpTyp[notExcludedAuthor])/np.sum(authorBooksTotal[notExcludedAuthor])
        if sliceSeparated > 1:
            slSep[nowDelta] = [{} for s in range(sliceSeparated)]
            for s in range(sliceSeparated):
                slSep[nowDelta][s] = {"all": {}}
                if nonAttri:
                    slSep[nowDelta][s]["exc"] = {}
                for attributionType in authorCorrectForSliceType:
                    authorCorrectForSliceType[attributionType][:, s] = np.ma.masked
                    authorBooksInSlice[:, s] = np.ma.masked
                    authorFragmentsInSlice[:, s] = np.ma.masked
                    tmpTyp = np.sum(authorCorrectForSliceType[attributionType], 1)
                    if attributionType == 'FRA':
                        slSep[nowDelta][s]['all']['FRA'] = np.sum(
                            tmpTyp[knownAuthor])/np.sum(authorFragmentsInSlice[knownAuthor, :])
                    else:
                        slSep[nowDelta][s]['all'][attributionType] = np.sum(
                            tmpTyp[knownAuthor])/np.sum(authorBooksInSlice[knownAuthor, :])
                    if nonAttri:
                        if attributionType == 'FRA':
                            slSep[nowDelta][s]["exc"]['FRA'] = np.sum(
                                tmpTyp[notExcludedAuthor])/np.sum(authorFragmentsInSlice[notExcludedAuthor, :])
                        else:
                            slSep[nowDelta][s]['exc'][attributionType] = np.sum(
                                tmpTyp[notExcludedAuthor])/np.sum(authorBooksInSlice[notExcludedAuthor, :])
                    authorCorrectForSliceType[attributionType].mask = np.ma.nomask
                    authorBooksInSlice.mask = np.ma.nomask
                    authorFragmentsInSlice.mask = np.ma.nomask

    return ret, slSep


def assign_unknown(attributions):
    unk = {}
    for tab in attributions:
        break
    for nowDelta in attributions[tab]:
        unk[nowDelta] = {}
        for ab in attributions:
            if ab[0] <= 0:
                unk[nowDelta][str(ab)] = {k: v for k,
                                          v in attributions[ab][nowDelta].items() if k != "FRA"}
    return unk


def convert(o):
    if isinstance(o, np.generic):
        return o.item()
    raise TypeError


def add_results_line(dataDir, ret, nonAttri, allowPartial, logDelta, unk, saveUnk=True, F=None, A=None):
    global _Q_excluded, _Q_association
    res_dic = {"ret": ret, "exc": list(set(_Q_excluded)), "noa": list(
        set(nonAttri)), "asso": _Q_association, "part": allowPartial, "delta": logDelta, "unknowns": unk if saveUnk else {}}
    if F is not None:
        res_dic["F"] = F
    if A is not None:
        res_dic["A"] = A
    with open(os.path.join(dataDir, "results.json"), "a") as fp:
        print(json.dumps(res_dic, default=convert), file=fp)


def previous_results(dataDir, nonAttri, logDelta, allowPartial, allDeltas=False, DBDir=None):
    res = {}
    unk = {}
    if os.path.isfile(os.path.join(dataDir, "results.json")):
        with open(os.path.join(dataDir, "results.json")) as fp:
            previous_results = fp.readlines()
        for line in previous_results:
            try:
                res_dic = json.loads(line)
            except:
                continue
            if {"part", "exc", "noa", "asso", "ret", "delta", "unknowns"} == set(res_dic.keys()):
                if set(res_dic["exc"]) == set(_Q_excluded) and (set(res_dic["noa"]) == set(nonAttri) and res_dic["part"] == allowPartial and res_dic["delta"] in logDelta or allDeltas):
                    if loadAssociation(res_dic["asso"], DBDir) == _Q_association:
                        res[res_dic["delta"]] = res_dic["ret"]
                        if res_dic["unknowns"]:
                            unk[res_dic["delta"]]={k:v for k, v in res_dic["unknowns"].items()}#tuple(map(int,k[1:-1].split(", ")))
            elif len(previous_results) == 1:
                with open(os.path.join(dataDir, "results.json"), "a") as fp:
                    print("", file=fp)
    return res, unk

def loadAssociation(asso, DBDir):
    """Returns the association from dict or file or an easily recognisable invalid value with some info."""
    if len(asso):
        if type(asso) == dict:
            tmp_association = asso
        elif type(asso) == str and DBDir is not None:
            try:
                tmp_association = get_association_from_file(DBDir, asso)
            except FileNotFoundError as e:
                return {"_association_": e}
        else:
            return {"_association_": TypeError("association must be dict or str")}
        return {int(a): b for a, b in tmp_association.items()}
    else:
        return {}

def availableResults (dataDir, DBDir=None):
    if os.path.isfile(os.path.join(dataDir, "results.json")):
        previous_settings = set()
        with open(os.path.join(dataDir, "results.json")) as fp:
            previous_results = fp.readlines()
        for line in previous_results:
            try:
                res_dic = json.loads(line)
            except:
                continue
            if {"part", "exc", "noa", "asso", "ret", "delta", "unknowns"} == set(res_dic.keys()):
                previous_settings.add((frozenset(res_dic["exc"]),tuple(sorted(loadAssociation(res_dic["asso"], DBDir).items()))))
        return previous_settings
    else:
        raise FileNotFoundError("No previous results to juice!")

def top_performers(_tarr, asso):
    tarr = structured_to_unstructured(_tarr)
    ret = [[[asso[int(tarr[-1][0])]], int(tarr[-1][1])], ]
    for a in tarr[-2::-1, :]:
        if a[1] < ret[0][1]:
            ret.append([[asso[int(a[0])]], int(a[1])])
            break
        ret[0][0].append(asso[int(a[0])])

    return ret


def print_margout(dataDir, attributions, nonAttri, logDelta, allowPartial, returnDict, nowDelta):
    global _Q_fra, _Q_excluded, _Q_association
    for tab in attributions:
        break
    if 'fn' in _Q_fra.columns:
        def fragnum(ab):
            return _Q_fra.loc[ab]['fn']
    else:
        def fragnum(ab): return 1
    baseMR = sequential_file_name(dataDir, "MR_Attribution")
    outpar = open(baseMR, "w")
    print_margout_head(set(_Q_excluded), set(nonAttri),
                       _Q_association, logDelta, file=outpar)
    if nonAttri:
        outnop = open(os.path.join(
            dataDir, f"{baseMR}_nonAttri.txt"), "w")
        print_margout_head(set(_Q_excluded), set(nonAttri),
                           _Q_association, logDelta, file=outnop)
    titles = {}
    for ab in sorted(attributions, key=lambda x: x[0]*1e6+x[1]):
        au, bo = ab
        if au in _Q_excluded:
            continue
        if not au in titles:
            titles[au] = []
        titles[au].append(bo)

        print_margout_book_line(
            ab, attributions[ab][nowDelta], fragnum(ab), outpar)
        if nonAttri and ab[0] not in nonAttri:
            print_margout_book_line(
                ab, attributions[ab][nowDelta], fragnum(ab), outnop)
    print_margout_summary_line(returnDict[nowDelta], outpar)
    if nonAttri:
        print_margout_summary_line(returnDict[nowDelta], outnop)
    if _Q_excluded:
        print_margout_summary_line(returnDict[nowDelta], outpar, "exc")
        if nonAttri:
            print_margout_summary_line(returnDict[nowDelta], outnop, "exc")
    for author in sorted(titles.keys(), key=lambda x: len(titles[x]), reverse=True):
        tmp = {}
        if author <= 0 or author in _Q_excluded:
            continue
        for typ in attributions[tab][nowDelta]:
            if typ == "FRA":
                tmp["FRA"] = sum([attributions[ab][nowDelta]["FRA"][0]
                                  for ab in attributions if ab[0] == author])
                continue
            tmtyre = 0
            for ab in attributions:
                if ab[0] == author and (typ in attributions[ab][nowDelta]) and (ab[0] in attributions[ab][nowDelta][typ][0][0]):
                    tmtyre += round(1 /
                                    len(attributions[ab][nowDelta][typ][0][0]), 5*allowPartial)
            tmp[typ] = tmtyre
        aF = sum([fragnum(ab) for ab in attributions if ab[0] == author])
        print_margout_author_line(author, scores=tmp, numBooks=len(
            titles[author]), numFragments=aF, file=outpar)
        if nonAttri and author not in nonAttri:
            print_margout_author_line(author, scores=tmp, numBooks=len(
                titles[author]), numFragments=aF, file=outnop)

    if nonAttri:
        outnop.close()
    outpar.close()


def sequential_file_name(dataDir, baseName, ext="txt"):
    nextName = baseName
    if os.path.isfile(os.path.join(dataDir, f"{nextName}.{ext}")):
        i = 0
        nextName = baseName+f"{i:02}"
        while os.path.isfile(os.path.join(dataDir, f"{nextName}.{ext}")):
            i += 1
            nextName = baseName+f"{i:02}"
    return os.path.join(dataDir, f"{nextName}.{ext}")


def print_margout_head(excluded, nonAttri, association, logDelta, file):
    print("excluded: ", set(excluded), file=file)
    print("nonAttri: ", set(nonAttri), file=file)
    print("associat: ", association, file=file)
    print("delta: ", logDelta, file=file)


def print_margout_book_line(authorBook, attribution, fragmentNum, file):
    print(f"A{authorBook[0]}B{authorBook[1]}", end="\t", file=file)
    for att in ["FNN", "TOP", "WP"]:
        res = attribution[att]
        try:
            print(
                f"{att}: {res[0][0]}: {res[0][1]:.4} {f'{res[1][0][0]}: {res[1][1]:.4}' if len(res[1][0]) else '--'}", end="\t", file=file)
        except IndexError:
            print(f"A{authorBook[0]}B{authorBook[1]} {att} {res}")
            raise
    for att in ["MR", "TMR", "WMR"]:
        res = attribution[att]
        try:
            print(
                f"{att}: {res[0][0]}: {res[0][1]} {f'{res[1][0][0]}: {res[1][1]}' if len(res[1][0]) else '--'}", end="\t", file=file)
        except IndexError:
            print(f"A{authorBook[0]}B{authorBook[1]} {att} {res}")
            raise
    print(f"FRA: {attribution['FRA'][0]/fragmentNum:.4}", file=file)


def print_margout_summary_line(returnDict, file, which='all'):
    print(f'\n\nSuccess Rate: ', end="\t", file=file)
    for r in returnDict[which].items():
        print(f"{r[0]}: {r[1]*100:.4}%", end="\t", file=file)
    print("\n\n", file=file)


def print_margout_author_line(author, scores, numBooks, numFragments, file):
    print(f"A{author} ({numBooks})", end="\t", file=file)
    for typ in scores:
        if typ != "FRA":
            print(
                f"{typ}: {round(scores[typ],1)} {scores[typ]/numBooks*100:.4} %", end="\t", file=file)
    print(
        f"FRA: {scores['FRA']} {scores['FRA']/numFragments*100:.4} %", file=file)


def check_sources(el_path):
    warnFlag = False
    source_info = {}
    with open(os.path.join(el_path, "src/bookprob.cpp")) as fp:
        for line in fp:
            if line.startswith("#define _P0_NORMALIZATION_"):
                a, b = line.split()[1:3]
                break
        else:
            a = "undefined in source"
            b = ""
    P0norm = " ".join([a, b])
    source_info["_P0_NORMALIZATION_"] = b
    logger.info(f"P0 version in use: {P0norm}")
    print(f"P0 version in use: {P0norm}")
    if os.path.getmtime(os.path.join(el_path, "bin/bookprob.o")) < os.path.getmtime(os.path.join(el_path, "src/bookprob.cpp")):
        logger.warning(
            "Source newer than executable, actual normalization may be different")
        warnFlag = True

    with open(os.path.join(el_path, "src/base_experiment.cpp")) as fp:
        for line in fp:
            if line.startswith("#define _USE_JOINT_PARAMS_"):
                source_info["_USE_JOINT_PARAMS_"] = True
                jointPar = line.split()[1]
                break
        else:
            jointPar = "NOT using joint params"
            source_info["_USE_JOINT_PARAMS_"] = False
    logger.info(f"Active: {jointPar}")
    print(f"Active: {jointPar}")
    if os.path.getmtime(os.path.join(el_path, "bin/fragprob.o")) < os.path.getmtime(os.path.join(el_path, "src/fragprob.cpp")):
        logger.warning(
            "Source newer than executable, actual normalization may be different")
        warnFlag = True

    with open(os.path.join(el_path, "src/authorSplitter.cpp")) as fp:
        for line in fp:
            if line.startswith("#define _FLEXIBLE_AUTHOR_SIZE_"):
                source_info["_FLEXIBLE_AUTHOR_SIZE_"] = True
                flexSize = line.split()[1]
                break
        else:
            flexSize = "NOT using flexible size"
            source_info["_FLEXIBLE_AUTHOR_SIZE_"] = False
    logger.info(f"Active: {flexSize}")
    print(f"Active: {flexSize}")
    if os.path.getmtime(os.path.join(el_path, "bin/authorSplitter.o")) < os.path.getmtime(os.path.join(el_path, "src/authorSplitter.cpp")):
        logger.warning(
            "Source newer than executable, actual normalization may be different")
        warnFlag = True
    if warnFlag:
        source_info["warnSourceNewer"] = True
    else:
        source_info["warnSourceNewer"] = False
    return source_info


def check_parameters(param):
    param.register_checks({o: lambda x: isinstance(
        x, (int, np.integer)) and x >= -1 for o in ['window', 'LZ77baseFragments']})
    param.register_checks(
        {o: lambda x: isinstance(x, (int, np.integer)) and x >= 0 for o in ['fragment', 'authorLength']})
    param.register_checks(
        ngramSize=lambda x: isinstance(x, (int, np.integer)) and 30 > x >= 0)
    param.register_checks(
        leaveNout=lambda x: isinstance(x, (int, np.integer)) and x >= 0)
    param.register_checks(
        delta=lambda x: (np.isscalar(x) and np.isreal(x) and np.isfinite(x)) or (not np.isscalar(x) and np.isreal(x).all() and np.isfinite(x).all()))
    param.register_checks(database=lambda x: os.path.isdir(x))
    param.register_checks(
        baseCorpusDir=lambda x: not x or os.path.isdir(x))
    param.register_checks(
        {o:lambda x: not x or os.path.isfile(x)} for o in ['configFile', 'P0file'])
    param.register_checks({o: lambda x: type(x) == bool for o in [
        'LZ77', 'retrieve', 'overwriteOutputs', 'keepTemporary', 'keepNonAlfabetical', 'mantainSlicing', 'allowPartial', 'dumpP0']})
    param.register_checks(suffix=lambda x: not x or (
        type(x) == str and "seq" not in x and "wnt" not in x))

    if param['LZ77'] + (param['ngramSize'] > 0) + param['keepNonAlfabetical'] > 1:
        raise ValueError("Only one of z or n or U options can be passed.")
    try:
        param.run_checks()
    except Exception as e:
        print(param)
        raise e


def seq_to_wnt(fragment: str, compression: bool, window: int, Encoding: str, translation):
    if compression:
        fragment = LZ77dict(
            fragment, window, Encoding)

    clean = ("" if fragment[0] == " " else " ")+fragment.translate(
        translation).replace(u"\u2026", "...")+("" if fragment[-1] == " " else " ")
    retWords = clean.split()
    return retWords


def create_tmpfile_dir(dataDir, tmpKind):
    outputDir = os.path.join(
        dataDir, tmpKind)

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
        logger.info(f"Created {outputDir}")

    return outputDir


def read_config(configPath=None, workers=4, Encoding="latin1"):
    if not configPath:
        if os.path.isfile("config.ini"):
            configPath = "config.ini"
        # looks for the config file in this sort of "PATH"
        elif os.path.isfile(os.path.join(el_path,"config.ini")):
            configPath = os.path.join(el_path,"config.ini")
        else:
            logger.info("Configuration file not found")
            print(
                f"Can't find config file {configPath}", file=sys.stderr)
            if support.get_input(f"Create a config file?"):
                logger.info("Creating a configuration file")
                configPath = support.build_config(el_path)
            else:
                logger.info("Using default values everywhere")
                print(f"Using default values everywhere.", file=sys.stderr)
                # If the user decides to go on with default values everywhere.
                configPath = None

    config = configparser.ConfigParser()
    config.read(configPath)
    logger.info("Loading configuration")
    if "elaboration" in config:
        if config["elaboration"].getint("workers") == None:
            logger.warning(
                "Missing 'workers' value in config file. Trying to deduce from numcomp config")
            print(
                "Missing 'workers' value in config file. Trying to deduce from numcomp config.", file=sys.stderr)
            if "numcomp" in config and "numCpus" in config["numcomp"] and "buffer" in config["numcomp"]:
                workers = config["numcomp"].getint(
                    "numCpus")-config["numcomp"].getint("buffer")
            else:
                logger.warnng(
                    f"Failed to deduce 'workers' number. Going on with default value {workers}.")
                print(
                    f"Failed to deduce 'workers' number. Going on with default value {workers}.", file=sys.stderr)
        else:
            workers = config["elaboration"].getint('workers')
        if config["elaboration"].get("encoding") == None:
            logger.warning(
                f"Missing 'encoding' value in config file. Going on with default value {Encoding}.")
            print(
                f"Missing 'encoding' value in config file. Going on with default value {Encoding}.", file=sys.stderr)
        else:
            Encoding = config["elaboration"]['encoding']
    else:
        logger.warning(
            "Missing 'elaboration' section in config file. Going on with default values")
        print(
            "Missing 'elaboration' section in config file. Going on with default values.", file=sys.stderr)

    return workers, Encoding, configPath


def foldsize(database, folds, authorFiles):
    fileList = [f for f in os.listdir(database) if f.endswith("txt")]
    if authorFiles:
        def buf_count_newlines_gen(fname):
            """https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python/68385697#68385697"""
            def _make_gen(reader):
                while True:
                    b = reader(2 ** 16)
                    if not b:
                        break
                    yield b

            with open(fname, "rb") as f:
                count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))
            return count
        N = 0
        for f in fileList:
            N += buf_count_newlines_gen(os.path.join(database, f))
        N /= 2
    else:
        N = len(fileList)

    def true_folds(N, l):
        f = N//l
        if N > l * np.sqrt(f * (f + 1)):
            return f+1
        else:
            return f
    l = int(N//folds)
    if l > 0:
        tf = true_folds(N, l)
        if tf == folds:
            return l
        else:
            tf2 = true_folds(N, l+1)
            if tf2 == folds:
                return l+1
            if abs(tf2-folds) < abs(tf-folds):
                tf = int(tf2)
                l += 1
            print(f"Too few documents for a nice {folds}-fold cross-validation.\n"
                  f"The closest thing is {tf}-fold cross-validation with {l}\n"
                  f"documents per fold.")
            return l
    else:
        print(f"Too few documents for {folds}-fold cross-validation.\n"
              f"Falling back to a leave-one-out", file=sys.stderr)
        return 1


def dirnames(param):
    if param['ngramSize']:
        methodPart = f"OSFG"
        midpart = f"{param['ngramSize']:02}N"
    elif param['keepNonAlfabetical']:
        methodPart = "KEEP"
        midpart = "DW"
    elif param['LZ77']:
        methodPart = "LZ77"
        midpart = "infBw" if param['window'] <= 0 else f"{param['window']:06}Bw"
    else:
        methodPart = "DICT"
        midpart = "WOR"
    leaveNoutPart = f"L{param['leaveNout']}O"
    suffixPart = param['suffix'] if param['suffix'] else "".join(
        map(chr, np.random.randint(65, 91, 5)))
    basename = "-".join([methodPart, midpart, leaveNoutPart, suffixPart])

    authorPart = "fullA" if param[
        'authorLength'] == 0 else f"{param['authorLength']:06}BA"

    fragmentPart = "infBF" if not param[
        'fragment'] else f"{param['fragment']:06}BF"

    return basename, authorPart, fragmentPart


def derive_translation(LZ77: bool = False, ngramSize: int = 0, keepUnder: bool = False):
    if LZ77:
        # for LZ77 the only spaces I need are those from newlines
        equivalen = support.equivalen[:-5] + "_"*5
        control = support.control
    else:
        control = support.control + ''.join([chr(i)
                                             for i in range(int('a', 16), int('e', 16))])
        equivalen = support.equivalen  # for ngram (and dict) I keep spaces
    if LZ77:
        translation = str.maketrans(
            support.unicpunt+control, equivalen+"_"*len(control))
    elif ngramSize or keepUnder:
        translation = str.maketrans(
            support.unicpunt+control, equivalen+" "*len(control))
    else:
        translation = str.maketrans(support.nonAlfa+control+support.unicpunt, ' ' *
                                    len(support.nonAlfa+control+equivalen))

    return translation
