###########################################################
## input: APL (Andromeda peak list) file                 ##
##                                                       ##
## very similar to mgf format, but using different field ##
## names and containing spectra from multiple raw files  ##
###########################################################


import sys
import os

from .spectrum import Spectrum
from .precursor import Precursor
from . import helpers


def parseAplHeaders(inputFile):
    with open(inputFile, 'r') as f:
        for line in f:
            if not line.startswith('peaklist start'):
                yield line
            else:
                break

def parse_apl(inputFile, storeFragmentIons = False):
    if not os.path.isfile(inputFile):
        sys.exit("ERROR: Could not open apl file: " + inputFile)
    
    if not inputFile.lower().endswith(".apl"):
        sys.exit("ERROR: File does not have the .apl file extension: " + inputFile)
    
    scannr = -1
    rtime = 0.0
    fragmentIons = list() if storeFragmentIons else None
    with open(inputFile, 'r') as f:
        for line in f:
            if line.startswith('peaklist start'):
                fragmentIons = list() if storeFragmentIons else None
            elif line.startswith('header'):
                header = line.split('=')[1]
                scannr = int(header.split("Index: ")[1].split()[0]) # e.g. header=RawFile: 01308_D03_P013387_B00_N20_R1 Index: 55569 Precursor: 0 _multi_
                raw_file = header.split("RawFile: ")[1].split()[0]
            elif line.startswith('mz'):
                fields = line.split('=')[1].split()
                pepmass = float(fields[0])
                intensity = None if len(fields) <= 1 else float(fields[1])
            elif line.startswith('charge'):
                charge = int(line.split('=')[1][0])
            elif line[0].isdigit():
                fragmentIons.append(map(float, line.split()))
            elif line.startswith('peaklist end'):
                if scannr != -1:
                    pepmz = helpers.precMzFromPrecMass(pepmass, charge)
                    yield Spectrum(scannr, [{'mz': pepmz, 'charge': charge, 'intensity': intensity}], fragmentIons, raw_file)

