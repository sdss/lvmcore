# !/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This is code to generate the "lvm_fidicual_fibermap.yaml" file.
# To use: run "python create_fibermap -f [path_to_excel_fibermap]"
#
import argparse
import datetime
import itertools
import os
import pathlib

import numpy as np
import pandas as pd
import yaml


# finblick - column 2 (index=2)
# fiberids - columns 3-21
# spec fibers - rows 4 - 40
# blocks
# spec1 blockid row 3
# spec2 blockid row 44
# spec3 blockid row 85


today = datetime.datetime.now().date().isoformat()


def total_fiber_num(ring: int) -> int:
    """ compute the total fiber number within a ring id """
    return 3 * ((ring - 1) ** 2) + 3 * (ring - 1) + 1


# basic definitions
n_ring = 25
ringnum = np.arange(n_ring) + 1
totals = total_fiber_num(ringnum)
segment_count = list(reversed(np.arange(1, n_ring)))  # fiber number in each hex segment per ring
rstart = totals[:-1] + 1  # starting fiber number at vertical x=0 in each ring (segment 1 of 6)


def read_excel_map(filepath: str) -> tuple:
    """ Read the excel file

    Read the excel file into a pandas DataFrame. Returns
    a tuple of dataframes, one for each sheet of the excel
    file.  Sheet 1, "Fibre Mapping", is the spectrograph
    fiber assignments. Sheet 2, "Science IFU Coordinates", is
    the science fiber xy position coordinates in mm.

    Parameters
    ----------
    filepath : str
        the input excel filepath

    Returns
    -------
    tuple
        two dataframes
    """
    df = pd.read_excel(filepath, 1)  # Fiber Mapping
    df.name = filepath.stem
    df.version = filepath.stem.rsplit('_v', 1)[1]
    df.created = pd.Timestamp.fromtimestamp(os.stat(filepath).st_ctime).date().isoformat()
    coords = pd.read_excel(filepath, 2)  # Science IFU Coordinates
    return df, coords


def build_ifu(coords: pd.DataFrame) -> pd.DataFrame:
    """ Build a df to represent the full science IFU

    Constructs a dataframe container for the full 1801
    science fiber IFU.  The central fiber is id 1.  Fiber id
    increases along the vertical x=0 line, then clockwise around
    the ring.  E.g. ring 1 contains fibers 2-7. ring 2 fibers 8-19, etc.
    Also assigns ring ids, a fiber number per ring, science x/y mm coords,
    and maps spectrographs fiber ids 1-600/601 to relative IFU
    fiber positions.

    Parameters
    ----------
    coords : pd.DataFrame
        the xy science fiber coords df (sheet 2 excel file)

    Returns
    -------
    pd.DataFrame
        the IFU df
    """
    fibnum_per_ring = totals - np.roll(totals, 1)
    ring = pd.DataFrame({"ringid": ringnum, "total": totals, "fibpr": fibnum_per_ring})
    ring.iloc[0, 2] = 1

    n_fibers = totals[-1]  # 1801 for 25 rings
    ringexp = [1] + ring["ringid"].iloc[1:].repeat(ring["fibpr"].iloc[1:]).tolist()
    fibid = np.arange(n_fibers) + 1
    ringfibnum = [1] + sum(
        (pd.RangeIndex(1, i + 1).tolist() for i in ring.iloc[1:, 2]), []
    )
    spec = pd.DataFrame({"fibid": fibid, "ringid": ringexp, 'ringfibnum': ringfibnum})

    tt = assign_fibers()
    ww = pd.DataFrame(tt, columns=['fid', 'specfib', 'specid'])
    ww = ww.sort_values('fid')
    spec = spec.merge(ww, left_on="fibid", right_on="fid")
    spec = spec.drop(columns="fid")
    spec['label'] = spec.apply(lambda x: f"S{x.specid}-{x.specfib}", axis=1)

    # add alt standards
    pos = build_standards_df()
    spec.loc[pos.fibid - 1, "altP"] = pos.label.tolist()

    # add alt sky
    sky = build_sky_df()
    sa = sky[sky.ifu == "A"]
    spec.loc[sa.fibid - 1, "altA"] = sa.label.tolist()
    sb = sky[sky.ifu == "B"]
    spec.loc[sb.fibid - 1, "altB"] = sb.label.tolist()

    # add x, y coords to ifu
    tmp = coords.set_index(["Sector", "Fiber #"])
    tt = spec.set_index(['specid', 'specfib'], drop=False)
    tt.loc[tmp.index, "xpmm"] = tmp["X [mm]"]
    tt.loc[tmp.index, "ypmm"] = tmp["Y [mm]"]
    spec = tt.reset_index(drop=True)
    return spec


def compute_standards_pos() -> list:
    """ Compute the IFU fiber id for each standard P1/2 fiber """

    # define some P ring number, index, and fiber spacing
    pr_num = np.array([24, 19])
    pr_idx = n_ring - pr_num
    pr_spacing = [[6, 5], [17, 13]]

    # for each IFU segment identify the P1/2 fiber positions
    total = []
    for seg in range(1, 7):
        # get the new segment fiber starting positions and find the standards pos
        start = np.array(list(reversed(rstart))) + (np.array(segment_count)) * (seg - 1)
        s = start[pr_idx] + pr_spacing

        # create a repeating spec id array
        specid = 1 if seg in [1, 2] else 2 if seg in [3, 4] else 3
        ss = [itertools.repeat(specid)] * 2

        # create a label array
        z = np.array([
            [f"P2-{seg+(seg-1)}", f"P1-{seg+(seg-1)}"],
            [f"P2-{seg+1+(seg-1)}", f"P1-{seg+1+(seg-1)}"]
            ]
                    )

        # zip the arrays together
        tmp = sum((list(i) for i in map(zip, z, s, ss)), [])
        total.extend(tmp)
    return total


def build_standards_df() -> pd.DataFrame:
    """ Builds a standards P fiber df

    Builds dataframe for the standard fibers,
    P1 and P2, with their label, fiber id within the
    IFU, ring number, spec id, etc.

    Returns
    -------
    pd.DataFrame
        the standards df
    """
    pr_num = np.array([24, 19])

    fiberpos = compute_standards_pos()
    df = pd.DataFrame(fiberpos, columns=['label', 'fibid', 'specid'])
    df['ring'] = df.label.apply(lambda x: pr_num[0] if 'P2' in x else pr_num[1])
    df[['name', 'num']] = df['label'].str.split('-', expand=True)
    return df


def chunk_spectro(specid: int = 1, seg: int = 1) -> list:
    """ Chunk the spectrograph fiber ranges by hex segment

    Maps the spectrograph fiber id to the IFU fiber id.
    Chunks the spectrograph fiber ranges 1-600/601 into 24,
    i.e. "n_ring - 1" sub-lists, each containing the fiber ids
    for the given spectograph hex segment in their respective rings.

    Parameters
    ----------
    specid : int, optional
        The spectrograph id, by default 1
    seg : int, optional
        IFU hex segment id, by default 1

    Returns
    -------
    list
        a list of spectrograph fiber ids chunked by ring and hex segment
    """
    # set range and identify segment spectrograph
    spec = np.arange(1, 601)
    seg1 = seg % 2 == 1
    schunk = spec[:300] if seg1 else spec[300:]

    if specid == 3 and not seg1:
        schunk += 1
    x = 0
    aa = []

    siter = range(24) if seg1 else range(23, -1, -1)

    # chunk the fibers
    for i in siter:
        odd = (i % 2) == 1
        tt = schunk[x : x + segment_count[i]]
        x += segment_count[i]
        aa.append(list(reversed(tt)) if odd else list(tt))

    return aa if seg % 2 == 1 else list(reversed(aa))


def chunk_fibers(seg: int = 1) -> list:
    """ Chunk the IFU fiber range by hex segment

    Chunks the IFU fiber range 1-1801 into 24, i.e. "n_ring - 1"
    sub-lists, each containing the fiber ids for
    the given IFU hex segment in their respective rings.

    Parameters
    ----------
    seg : int, optional
        IFU hex segment id, by default 1

    Returns
    -------
    list
        a list of IFU fiber ids chunked by ring and hex segment
    """
    start = np.array(list(reversed(rstart))) + (np.array(segment_count)) * (seg - 1)
    return [list(start[i] + np.arange(segment_count[i])) for i in range(24)]


def chunk_sky(specid: int = 1, label: str = 'A') -> list:
    """ Chunk the sky fiber ranges by specid

    Maps the sky spectro. fiber id to the IFU fiber id.
    Chunks the sky fiber ranges 1-19/20/21 into tuples,
    each containing the fiber ids, the sky fiber id, the
    spectrograph id, and the sky IFU label.

    Parameters
    ----------
    specid : int, optional
        The spectrograph id, by default 1
    label : str, optional
        The sky IFU label A or B, by default A

    Returns
    -------
    list
        a list of sky fiber ids for each sky IFU
    """
    sky_ring = 5
    skystart = rstart[:sky_ring - 1]
    skyseg = segment_count[-(sky_ring - 1):]
    seg = 1 if specid == 1 else 3 if specid == 2 else 5

    nfib = 20 if specid == 1 or (specid == 2 and label == 'B') else 21

    srange = np.arange(1, nfib)
    if specid == 3:
        srange += 1

    start = np.array(list(reversed(skystart))) + (np.array(skyseg)) * (seg - 1)

    chunked = np.array_split(srange, 5)  # number of rows of skies

    # manual fiber mapping
    start = np.array(list(reversed(skystart))) + (np.array(skyseg)) * (seg - 1)
    a=[list(start[i] + np.arange(skyseg[i])) for i in range(4)]
    start = np.array(list(reversed(skystart))) + (np.array(skyseg)) * (seg)
    b=[list(start[i] + np.arange(skyseg[i])) for i in range(4)]
    aa=sum(a, [])
    bb=sum(b, [])
    fibers = [(aa[0], aa[4], aa[7], aa[9]), (aa[1], aa[5], aa[8], bb[9]),
              (aa[2], aa[6], bb[7], bb[8]), (aa[3], bb[4], bb[5], bb[6]),
              (bb[0], bb[1], bb[2], bb[3])]
    for i in range(0, 5, 2):
        fibers[i] = tuple(reversed(fibers[i]))

    ss = [itertools.repeat(specid)] * len(fibers)
    ll = [itertools.repeat(label)] * len(fibers)

    return sum((list(i) for i in map(zip, fibers, chunked, ss, ll)), [])


def build_sky_df() -> pd.DataFrame:
    """ Builds a sky A/B fiber df

    Builds dataframe for the sky ifus,
    A and B, with their label, fiber id within the
    IFU, spec id, etc.

    Returns
    -------
    pd.DataFrame
        the sky df
    """

    # create A sky ifu
    tmp = [(1, 1, 3, 'A')]
    for i in range(1, 4):
        tmp.extend(chunk_sky(specid=i, label='A'))

    # create B sky ifu
    tmp.extend([(1, 1, 3, 'B')])
    for i in range(1, 4):
        tmp.extend(chunk_sky(specid=i, label='B'))

    # combine together
    df = pd.DataFrame(tmp, columns=['fibid', 'skyfib', 'specid', 'ifu'])
    df['label'] = df.apply(lambda x: f'{x.ifu}{x.specid}-{x.skyfib}', axis=1)
    return df


def assign_fibers() -> list:
    """ Assign the spectrograpgh fibers to the IFU fibers """

    total = [(1, 301, 3)]
    for seg in range(1, 7):
        specid = 1 if seg in [1,2] else 2 if seg in [3,4] else 3
        fibers = chunk_fibers(seg=seg)
        spec = chunk_spectro(specid=specid, seg=seg)
        ss = [itertools.repeat(specid)] * len(spec)
        tmp = sum((list(i) for i in map(zip, fibers, spec, ss)), [])
        total.extend(tmp)
    return total


def create_spectro_df(df: pd.DataFrame, specid: int = 1, coords: pd.DataFrame = None) -> pd.DataFrame:
    """ Create a df for the specified spectrograph

    Parameters
    ----------
    df : pd.DataFrame
        the input fiber df (sheet 1 excel df)
    specid : int, optional
        the spectrograph id, by default 1
    coords : pd.DataFrame, optional
        the input xy coords df (sheet 2 excel df), by default None

    Returns
    -------
    pd.DataFrame
        output fiber df for a spectrograph
    """

    # set start column to search for slit title
    start_col = 3 if df.version == '0.3' else 2

    # get the starting spectrograph rows
    idxes = df.index[df.iloc[:, start_col].str.contains('SPECTROGRAPH SLIT').fillna(False)]
    if len(idxes) < 3:
        # insert a starting index if there is less than 3
        idxes = idxes.insert(0, -1)

    # get the starting block and and fiber rows for the input spectrograph id
    block_idxes = (idxes + 1)
    start_idx = block_idxes[specid - 1]
    fib_idx = start_idx + 1

    # extract the blocks in columns 3-20 (v0.3) or 2-19 (v0.5)
    n_blocks = 18
    blocks = df.iloc[start_idx][start_col:start_col + n_blocks].tolist()

    # extract the fiber num within the block
    n_finblock = 36
    finblock = df.iloc[fib_idx:fib_idx + n_finblock, start_col - 1].tolist() * n_blocks

    # extract fiber block
    sub = df.iloc[fib_idx:fib_idx + n_finblock, start_col:start_col + n_blocks]
    fibers = sub.unstack()
    n_fibs = len(fibers)

    # expand the blocks
    blocks = df.iloc[start_idx, start_col:start_col + n_blocks].apply(lambda x: x.split(f'S{specid}')[1]).repeat(n_finblock)

    # expand the spectrograph id
    spectroid = [specid] * n_fibs

    # create a new df
    new = pd.DataFrame(sub.unstack(), columns=['label'])
    new['spectrographid'] = spectroid
    new['blockid'] = blocks.tolist()
    new['finblock'] = finblock
    new['targettype'] = new['label'].apply(lambda x: 'science' if x.startswith('S') else 'standard' if x.startswith('P') else 'SKY')
    new[['ifulabel', 'finifu']] = new['label'].str.split('-', expand=True)
    new['finifu'] = new.finifu.astype(int)

    # add telescope column
    # TODO - check which sky IFU A, B is SkyE and SkyW
    tmap = {'S': 'Sci', 'P': 'Spec', 'A': 'SkyW', 'B': 'SkyE'}
    new['telescope'] = new['ifulabel'].apply(lambda x: tmap[x[0]])

    # return nothing if no coords provided
    if coords is None:
        return new

    # add in the x and y pos science coords
    sc = coords['Sector'] == specid
    new['xpmm'] = new['ypmm'] = ''
    new.loc[new['targettype'] == 'science', 'xpmm'] = coords[sc].iloc[:, 0].tolist()
    new.loc[new['targettype'] == 'science', 'ypmm'] = coords[sc].iloc[:, 1].tolist()

    return new


def set_xypos(col: str, main: pd.DataFrame, ifu: pd.DataFrame) -> pd.DataFrame:
    """ Sets the xypmm columns

    Updates the x/ypmm columns for the sky and standard fibers
    using the altP, altA, and altB columns.  It uses their relative
    positions equivalent to the same science fiber at the same id
    location.

    Parameters
    ----------
    col : str
        the name of the column to match
    main : pd.DataFrame
        the input fiber df
    ifu : pd.DataFrame
        the input IFU df

    Returns
    -------
    pd.DataFrame
        updated fiber df
    """

    # identify the locations of the sky and standard with
    # respect to their science fiber locations
    sub = ifu[~ifu[col].isna()]
    sub = sub.set_index(col, drop=False)
    sub = sub.reindex(index=main[main.label.isin(ifu[col])].label)

    # update the xy coord columns with the relevant ones from ifu
    main = main.set_index("label", drop=False)
    main.loc[sub.index, 'xpmm'] = sub["xpmm"].tolist()
    main.loc[sub.index, 'ypmm'] = sub["ypmm"].tolist()
    main = main.reset_index(drop=True)
    return main


def set_ring(col: str, main: pd.DataFrame, ifu: pd.DataFrame) -> pd.DataFrame:
    """ Sets the ringnum column

    Updates the ringnum column for the sky and standard fibers
    using the altP, altA, and altB columns.  It uses their relative
    positions equivalent to the same science fiber at the same id
    location.

    Parameters
    ----------
    col : str
        the name of the column to match
    main : pd.DataFrame
        the input fiber df
    ifu : pd.DataFrame
        the input IFU df

    Returns
    -------
    pd.DataFrame
        updated fiber df
    """

    # identify the locations of the sky and standard with
    # respect to their science fiber locations
    sub = ifu[~ifu[col].isna()]
    sub = sub.set_index(col, drop=False)
    sub = sub.reindex(index=main[main.label.isin(ifu[col])].label)

    # update the ringnum column with the relevant ring id
    main = main.set_index("label", drop=False)
    main.loc[sub.index, 'ringnum'] = sub["ringid"].tolist()
    main = main.reset_index(drop=True)
    return main


def build_fiber_data(main: pd.DataFrame, coords: pd.DataFrame, ifu: pd.DataFrame) -> pd.DataFrame:
    """ Build the output fiber map dataframe

    Reformats the input excel dataframe into a single fibermap
    dataframe.  Combines the 3 spectrograph science fibermaps.
    Adds a global fiberid column. Updates the xy positions for
    the sky and standard fibers. Assigns ringnum and fiber
    status columns for all the fibers.

    Parameters
    ----------
    main : pd.DataFrame
        the input fibers df (excel sheet 1)
    coords : pd.DataFrame
        the science xy coords df (excel sheet 2)
    ifu : pd.DataFrame
        the input IFU dataframe

    Returns
    -------
    pd.DataFrame
        the final fibermap df
    """

    # generate science fiber dfs for each spectrograph
    spec1 = create_spectro_df(main, 1, coords)
    spec2 = create_spectro_df(main, 2, coords)
    spec3 = create_spectro_df(main, 3, coords)

    # create a combined science fiber dataframe
    final = pd.concat([spec1, spec2, spec3])
    final = final.reset_index(drop=True)

    # update sky and standard xy positions
    final = set_xypos('altP', final, ifu)
    final = set_xypos('altA', final, ifu)
    final = set_xypos('altB', final, ifu)

    # add ringnum column
    ifu = ifu.set_index(["label"], drop=False)
    tmp = final.set_index(['label'], drop=False)
    tmp.loc[ifu.index, "ringnum"] = ifu.loc[ifu.index].ringid.tolist()
    tmp = set_ring('altP', tmp, ifu)
    tmp = set_ring('altA', tmp, ifu)
    tmp = set_ring('altB', tmp, ifu)
    final = tmp.reset_index(drop=True)

    # deal with broken, dead, misrouted fibers
    broke = pd.read_table('metrology/fiber_status.dat', sep=',', header=17,
                          names=['specid', 'label', 'finifu', 'status', 'misroute', 'blockid', 'finblock'])
    broke['label'] = broke['label'].str.strip()
    broke['statnum'] = pd.Categorical(broke['status'].str.strip(),
                                      categories=['ok', 'dead', 'low', 'repair']).codes
    broke = broke.set_index(['specid', 'label', 'finifu'])

    # add fiber status column
    tmp = final.set_index(['spectrographid', 'ifulabel', 'finifu'], drop=False)
    tmp['fibstatus'] = 0
    tmp.loc[broke.index, 'fibstatus'] = broke.loc[broke.index, 'statnum']

    # check misrouted fibers
    misb = broke[broke['misroute'] == 1]
    bb = tmp.loc[misb.index]['finblock'] == misb['finblock']
    if not bb.all():
        tmp.loc[misb.index, 'finblock'] = misb['finblock']
        bb = tmp.loc[misb.index]['finblock'] == misb['finblock']
        if not bb.all():
            raise ValueError('Mismatch between misrouted fibers in broken status table and main fiber table.  Double check correct fiber numbers and blocks.')

    # reset index
    final = tmp.reset_index(drop=True)

    # sort by specid, blockid, finblock to reorder misrouted fibers
    final['block'] = final['blockid'].str[1:].astype(int)
    final = final.sort_values(['spectrographid', 'block', 'finblock'])

    # add fiberid column
    final.insert(0, 'fiberid', pd.RangeIndex(len(final)) + 1)

    # drop unneeded columns
    final = final.drop(columns=["label", "block"])

    return final


class Dumper(yaml.Dumper):
    """ Custom YAML dumper class """
    def increase_indent(self, flow=False, *args, **kwargs):
        """ increase the sequence indent """
        return super().increase_indent(flow=flow, indentless=False)

    def represent_mapping(self, tag, mapping, **kwargs):
        """ switch the mapping flow style to False """
        return super().represent_mapping(tag, mapping, flow_style=False)


def write_yaml(output: str = "lvm_fiducial_fibermap.yaml",
               final: pd.DataFrame = None, comments: str = None):
    """ Write out a YAML file

    Write out the fiber mapping into a machine-readable
    YAML file.

    Parameters
    ----------
    output : str, optional
        the name of the output file, by default "lvm_fiducial_fibermap.yaml"
    final : pd.DataFrame, optional
        the fiber map dataframe, by default None
    comments : str, optional
        the yaml file comments, by default None
    """
    # write the output yaml file
    with open(output, "w") as f:
        yaml.dump(
            {"schema": schema, "fibers": final.values.tolist()},
            f,
            default_flow_style=None, Dumper=Dumper, sort_keys=False
        )

    # do nothing if no comments
    if not comments:
        return

    # write comments
    with open(output, "r+") as f:
        data = f.read()
        f.seek(0)
        f.write(comments + "\n" + data)


# schema description of columns in YAML file

schema = [
 {'name': 'fiberid',
  'dtype': 'int',
  'description': 'number of the fiber along the slithead',
  'unit': None},
 {'name': 'spectrographid',
  'dtype': 'int',
  'description': 'the spectrograph id number, either 1, 2 or 3',
  'unit': None},
 {'name': 'blockid',
  'dtype': 'str',
  'description': 'the ID label of the block along the slit',
  'unit': None},
 {'name': 'finblock',
  'dtype': 'int',
  'description': 'the fiber number within the v-groove block',
  'unit': None},
 {'name': 'targettype',
  'dtype': 'str',
  'description': 'the type of fiber, either science, standard or sky',
  'unit': None},
 {'name': 'ifulabel',
  'dtype': 'str',
  'description': 'the ID label for the IFU + spectrograph',
  'unit': None},
 {'name': 'finifu',
  'dtype': 'int',
  'description': 'the fiber number within the IFU',
  'unit': None},
 {'name': 'telescope',
  'dtype': 'str',
  'description': 'the name of the telescope; Sci, Spec, SkyE/W for science, standards, or skies',
  'unit': None},
 {'name': 'xpmm',
  'dtype': 'float',
  'description': 'the x coordinate in mm of the fiber relative to the centroid',
  'unit': 'mm'},
 {'name': 'ypmm',
  'dtype': 'float',
  'description': 'the x coordinate in mm of the fiber relative to the centroid',
  'unit': 'mm'},
 {'name': 'ringnum',
  'dtype': 'int',
  'description': 'the number of the IFU ring the fiber belongs to',
  'unit': None},
 {'name': 'fibstatus',
  'dtype': 'int',
  'description': 'the status of the fiber, 0=live, 1=dead, 2=low, 3=repair, 4=short',
  'unit': None}]


# comments for YAML file
def create_comments(filename, df):

    if 'Fibre mapping' in df.columns:
        version = df['Fibre mapping'][0]
        date = df['Fibre mapping'][1]
    else:
        version = f'Version {df.version}'
        date = f'Date: {df.created}'
    return rf"""
# LVM fiber mapping file: lvm_fiducial_fibermap
#
# This file is a machine-readble translation of the fiducial fiber mapping
# information specified in the given excel sheet.
#   Filename: {filename}
#   Revision: {version}
#   Source File: {date}
#   Yaml File Generated: {today}
#   Generated by: create_fibermap.py

# The data structure is as follows.
# 1) fiberid: Number of the fiber along the slithead (should be 1 --> 1944 where 1944
#    is the total number of fibers, 1--> 648 should be on
#    spectrograph 1, 649 --> 1296 should be on spectrograph 2, 1297 --> 1944 should
#    be on spectrograph 3
# 2) spectrographid: Either 1, 2, or 3
# 3) blockid: ID label of block along the slit
# 4) finblock: Fiber number within the v-groove block
# 5) targettype: science, standard, or sky
# 6) ifulabel: ID label for the IFU + spectrograph
# 7) finifu: Fiber number within the IFU
# 8) telescope: The name of the telescope; Sci, Spec, SkyE/W for science, standards, or skies
# 9) xpmm: The x coordinate in mm of the fiber relative to the centroid
# 10) ypmm: The y coordinate in mm of the fiber relative to the centroid
# 11) ringnum: Number of the IFU ring the fiber belongs to
# 12) fibstatus: 0=live, 1=dead, 2=low, 3=repair, 4=short
        """


def create_fibermap(filename: str, output: str):
    """ Create the fibermap YAML file

    Parameters
    ----------
    filename : str
        the path name of the fibermap excel file
    output : str
        the output YAML filename
    """
    # resolve the filepath
    filepath = pathlib.Path(filename).resolve()
    filename = filepath.name

    # read the excel file
    fmap, coords = read_excel_map(filepath)

    # create YAML comments
    comments = create_comments(filename, fmap)

    # build the input IFU and output fiber mapping dataframes
    ifu = build_ifu(coords)
    fibers = build_fiber_data(fmap, coords, ifu)

    # write the fibers out to a YAML file
    write_yaml(output, final=fibers, comments=comments)


def set_column(col: str, main: pd.DataFrame, ifu: pd.DataFrame,
               main_col: str, ifu_col: str) -> pd.DataFrame:
    """ Updates a column in main fiber df with values from ifu df

    Updates a column in main fiber df with values from the
    matching column in the ifu df.  ``main_col`` specifies the
    main column to update.  ``ifu_col`` specifies the ifu column
    to pull values from.  ``col`` specifies the name of the column
    of fibers you want to match, either "label" for science fibers,
    "altP" for standards, "altA" for sky IFU A, or "altB" for sky
    IFU B.  It uses their relative positions equivalent to the
    same science fiber at the same id location.

    Parameters
    ----------
    col : str
        the name of the column to match
    main : pd.DataFrame
        the input fiber df
    ifu : pd.DataFrame
        the input IFU df
    main_col : str
        the name of the column in the main fibers df
    ifu_col : str
        the name of the column in the main IFU df

    Returns
    -------
    pd.DataFrame
        updated fiber df
    """
    # identify the locations of the sky and standard with
    # respect to their science fiber locations
    sub = ifu[~ifu[col].isna()]
    sub = sub.set_index(col, drop=False)
    sub = sub.reindex(index=main[main.label.isin(ifu[col])].label)

    # update the column with the relevant IFU values
    main = main.set_index("label", drop=False)
    main.loc[sub.index, main_col] = sub[ifu_col].tolist()
    main = main.reset_index(drop=True)
    return main


def insert_into_db(filename: str = None, ifu: pd.DataFrame = None, fibers: pd.DataFrame = None):
    """ Insert static fiber data into the LVM db

    Inserts the static fiber data from the IFU and fibers dataframes
    into the lvmdb drp schema, tables "ifu", and "fibers".

    Parameters
    ----------
    filename : str, optional
        the name of the fibermap excel file, by default None
    ifu : pd.DataFrame, optional
        the IFU dataframe, by default None
    fibers : pd.DataFrame, optional
        the fibermap dataframe, by default None
    """
    from sdssdb.sqlalchemy.lvmdb import database

    if filename and (ifu is not None or fibers is not None):
        fmap, coords = read_excel_map(filename)
        ifu = build_ifu(coords)
        ifu = add_radec_offs(ifu=ifu)
        fibers = build_fiber_data(fmap, coords, ifu)

    # insert ifu table
    sub = ifu.rename(columns={'fibid': 'fiberid', 'specfib': 'specfibid',
                              'altP': 'standard', 'altA': 'skya', 'altB': 'skyb'})
    sub.to_sql('ifu', database.engine, schema='drp', index=False, if_exists='append')

    # fill out ifu_pk column
    fibers['label'] = fibers.apply(lambda x: f'{x.ifulabel}-{x.finifu}', axis=1)
    sub = fibers.set_index('label', drop=False)
    tmp = ifu.set_index('label', drop=False)
    sub = set_column('altP', sub, tmp, main_col='ifu_pk', ifu_col='fibid')
    sub = set_column('altA', sub, tmp, main_col='ifu_pk', ifu_col='fibid')
    sub = set_column('altB', sub, tmp, main_col='ifu_pk', ifu_col='fibid')
    sub = set_column('label', sub, tmp, main_col='ifu_pk', ifu_col='fibid')
    fibers = sub.reset_index(drop=True).drop(columns="label")

    # insert fibers table
    sub = fibers.rename(columns={'spectrographid': 'specid', 'fibstatus': 'status'})
    sub.to_sql('fibers', database.engine, schema='drp', index=False, if_exists='append')


def load_offsets_to_db():
    """ Loads the RA/Dec offsets to IFU table

    This can be used to load the RA/Dec offsets from the
    simbmap file to the drp.IFU table independently from the rest
    of the IFU columns.  This is useful if you don't want to
    have the LVM datasimulator package and the db package
    dependencies in the same env.
    """
    from sdssdb.sqlalchemy.lvmdb import database, drp
    from sqlalchemy import update
    df = pd.read_table('lvm_simbmap_1801.dat', sep=',', header=21)
    sub = df[['fiberid', 'raoff', 'decoff']]
    out = sub.rename(columns={'fiberid': 'pk'}).to_dict('records')
    with database.Session() as session, session.begin():
        session.execute(update(drp.IFU), out)


def add_radec_offs(ifu: pd.DataFrame = None) -> pd.DataFrame:
    """ Add RA/Dec offsets to the IFU dataframe

    Adds the RA/Dec relative fiber offsets in units
    of arcsec to the IFU bundle dataframe.  Offsets are
    taken from the input science_array.dat file from the
    LVM data simulator bundle, and mapped onto the existing
    IFU dataframe fiber ordering. Note: running this
    method requires the LVM datasimulator package to be installed.

    Parameters
    ----------
    ifu : pd.DataFrame, optional
        the input IFU df, by default None

    Returns
    -------
    pd.DataFrame
        the updated IFU df
    """
    # create a fiber bundle, rotated to match the excel IFU spec
    from lvmdatasimulator.fibers import FiberBundle
    bundle = FiberBundle(bundle_name='full', nrings=25, angle=-90)

    # convert to a df
    df = bundle.fibers_table_science.to_pandas()

    # round the x and y arcsec offsets so pandas can properly sort values
    sub = df.copy()
    sub = sub.round({'x': 4, 'y': 4})

    # sort the bundle df in descending x (ra) and ascending y (dec)
    bsort = sub.sort_values(['x', 'y'], ascending=[False, True])

    # sort the IFU df in ascending x, y
    ifusort = ifu.sort_values(['xpmm', 'ypmm'], ascending=[True, True])

    # map the columns, default matches by df index
    bb = bsort.set_index('ring_id', drop=False)
    ii = ifusort.set_index('ringid', drop=False)
    ii.loc[bb.index, 'raoff'] = bb.loc[bb.index, 'x'].tolist()
    ii.loc[bb.index, 'decoff'] = bb.loc[bb.index, 'y'].tolist()

    # reset back to the original IFU df sort
    return ii.sort_values('fibid').reset_index(drop=True)


def create_simbmap(filename: str, output: str = 'lvm_simbmap_1801.dat'):
    """
    """
    # resolve the filepath
    filepath = pathlib.Path(filename).resolve()
    filename = filepath.name

    # read the excel file
    fmap, coords = read_excel_map(filepath)

    # build the IFU and add the RA/Dec offsets
    ifu = build_ifu(coords)
    ifu = add_radec_offs(ifu=ifu)

    # create a subset to write out
    sub = ifu[['fibid', 'ringid', 'ringfibnum', 'raoff', 'decoff']]
    sub = sub.rename(columns={'fibid': 'fiberid', 'ringfibnum': 'ringfnum'})

    comments = r"""
# LVM fiducial metrology file
#
# The purpose of this file is to give fiducial positions relative to the center of an IFU
# for fibers in a given bundle.  The fiber numbering scheme here is clockwise starting
# from the center, incrementing along the vertical x=0 axis, using the Excel IFU fiber map
# schematic as guide. The data structure is as follows.
#
# Data table entries:
# 1) fiberid: Fiber number within the bundle.
# 2) ringid: Ring id number of the fiber
# 3) ringfnum: Fiber number within a ring, starting from 1, at the vertical x=0 position
# 4) raoff: The ra offset of the fiber in arcsec relative to the ifu center
# 5) decoff: The dec offset of the fiber in arcsec relative to the ifu center
#
# raoff and decoff are both given in arcseconds, so they should be added to the central pointing to get
# the effective on-sky coordinates of the fiber
# E.g., if RA=213.3423D and DEC=52.3234D in decimal degrees (ALWAYS use double precision here!)
# then RAfiber=RA+raoff/3600.D/cos(DEC*!DPI/180.)    DECfiber=DEC+decoff/3600.D
#
# Assumes that fiberdiameter=35.3 arcsec
#
    """

    # write the fibers out to a dat file
    with open(output, 'w') as f:
        f.write(comments + '\n\n')
        sub.to_csv(f, sep=',', index=False)


parser = argparse.ArgumentParser(
    prog='create_fibermap',
    description='Creates the fiducial fibermap machine-readable YAML file')

parser.add_argument('-f', '--filename', help='the input fibermap excel file', default="SDSS-V_QQQQ-LVMI-fiber-mapping_v0.5.xlsx")
parser.add_argument('-o', '--output', help='the output YAML filename', default="lvm_fiducial_fibermap.yaml")

if __name__ == '__main__':
    # parse arguments
    opts = parser.parse_args()
    create_fibermap(opts.filename, opts.output)
