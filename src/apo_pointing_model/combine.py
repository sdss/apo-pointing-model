#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Conor Sayres (csayres@uw.edu)
# @Date: 2024-08-29
# @Filename: combine.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import datetime

import polars


def getTpointHeader():
    """Returns list of lines to include in TPOINT header"""
    currDateStr = datetime.datetime.now().isoformat()
    header = [
        "! Caption record:",
        "SDSS APO 2.5m Pointing Data " + currDateStr,
        "! Option record: ALTAZ followed by rotator code, if applicable",
        ": ALTAZ",
        "! Run parameters: telescope latitude deg min sec",
        "32 46 49.3",
        "! Pointing data (TPOINT format #4):",
        "!",
        "! All angles are in degrees",
        "! Azimuth uses the convention: N=0, E=90 (unlike the TCC, which uses S=0, E=90)",
        "!",
        "!  Desired Phys            Actual Mount          Rot Phys",
        "! Az          Alt         Az          Alt",
    ]
    return header


def processFiles(fileInputList, fileOutput):
    """Get a TPOINT formatted file from one or more pointing data csv files

    Inputs:
    - file_input_list: list of paths to csv files to be used
    - file_output: path to output file

    Returns None
    """
    dfList = []
    for fin in fileInputList:
        dfList.append(polars.read_csv(fin))
    df = polars.concat(dfList)
    nMeas = len(df)
    df = df.filter(~polars.col("failed"))
    nGood = len(df)
    print("using %i of %i pointings" % (nGood, nMeas))

    # convert from tcc az convention to TPOINT az convention (N=0, E=90)
    desPhysAz = 180 - df["ptdata_azphys"].to_numpy()
    desPhysAlt = df["ptdata_altphys"].to_numpy()
    mountPosAz = 180 - df["ptdata_azmount"].to_numpy()
    mountPosAlt = df["ptdata_altmount"].to_numpy()
    rotPhys = df["ptdata_rotphys"].to_numpy()

    fileLines = getTpointHeader()

    for ii in range(nGood):
        fileLines.append(
            "%11.6f %11.6f %11.6f %11.6f %10.5f"
            % (
                desPhysAz[ii],
                desPhysAlt[ii],
                mountPosAz[ii],
                mountPosAlt[ii],
                rotPhys[ii],
            )
        )

    print("writing", fileOutput)
    with open(fileOutput, "w") as f:
        for line in fileLines:
            f.write(line + "\n")
