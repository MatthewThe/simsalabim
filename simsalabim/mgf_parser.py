###########################################################
## input: MGF file                                       ##
##                                                       ##
## http://www.matrixscience.com/help/data_file_help.html ##
###########################################################


import sys
import os

from .spectrum import Spectrum
from .precursor import Precursor
from . import helpers


def parseMgfHeaders(inputFile):
  with open(inputFile, 'r') as f:
    for line in f:
      if not line.startswith('BEGIN IONS'):
        yield line
      else:
        break

def parse_mgf(inputFile, storeFragmentIons = False):
  if not os.path.isfile(inputFile):
    sys.exit("ERROR: Could not open mgf file: " + inputFile)
  
  if not inputFile.lower().endswith(".mgf"):
    sys.exit("ERROR: File does not have the .mgf file extension: " + inputFile)
  
  scannr = -1
  rtime = 0.0
  fragmentIons = list() if storeFragmentIons else None
  with open(inputFile, 'r') as f:
    for line in f:
      if line.startswith('TITLE'):
        fields = line.split('=')
        if len(fields) > 2:
          scannr = int(line.split('=')[2].split("\"")[0])
        else:
          #scannr = int(line.split('=')[1].split("[")[0])
          scannr = int(line.split('=')[1].split(".")[1])
      elif line.startswith('SCANS'):
        scannr = int(line.split('=')[1])
      elif line.startswith('RTINSECONDS'):
        rtime = float(line.split('=')[1])
      elif line.startswith('PEPMASS'):
        fields = line.split('=')[1].split()
        pepmz = float(fields[0])
        intensity = None if len(fields) <= 1 else float(fields[1])
      elif line.startswith('CHARGE'):
        charge = int(line.split('=')[1][0])
      elif line[0].isdigit():
        fragmentIons.append(map(float, line.split()))
      elif line.startswith('END IONS'):
        if scannr != -1:
          yield Spectrum(scannr, [{'mz': pepmz, 'charge': charge, 'intensity': intensity}], fragmentIons)
        fragmentIons = list() if storeFragmentIons else None

