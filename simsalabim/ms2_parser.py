###############################################################
## input: MS2 file                                           ##
##                                                           ##
## http://crux.ms/file-formats/ms2-format.html (base format) ##
## http://crux.ms/commands/bullseye.html (with EZ lines)     ##
###############################################################


import sys
import os

from .spectrum import Spectrum
from .precursor import Precursor
from . import helpers
    
def parseMs2Headers(inputFile):
  with open(inputFile, 'r') as f:
    for line in f:
      if line.startswith('H'):
        yield line
      else:
        break


def parse_ms2(inputFile, storeFragmentIons = False):  
  if not os.path.isfile(inputFile):
    sys.exit("ERROR: Could not open ms2 file: " + inputFile)
  
  if not inputFile.lower().endswith(".ms2"):
    sys.exit("ERROR: File does not have the .ms2 file extension: " + inputFile)
  
  scannr = -1
  ezLines = list()
  fragmentIons = list() if storeFragmentIons else None
  hasEZ = False
  with open(inputFile, 'r') as f:
    for line in f:
      if line.startswith('S'):
        if scannr != -1:
          yield Spectrum(scannr, ezLines, fragmentIons)
        fragmentIons = list() if storeFragmentIons else None
        scannr = int(line.split('\t')[1])
        ezLines = list()
      elif line.startswith('Z'):
        if not hasEZ:
          charge = int(line.split('\t')[1])
          pepmass = float(line.split('\t')[2])
          pepmz = precMzFromPrecMass(pepmass, charge)
          intensity = None
          ezLines.append({'mz': pepmz, 'charge': charge, 'intensity': intensity})
      elif line.startswith('I\tRTime'):
        rtime = float(line.split('\t')[2])
      elif line.startswith('I\tEZ'):
        fields = line.split('\t')
        charge = int(fields[2])
        pepmass = float(fields[3])
        rtimeTmp = float(fields[4])
        intensity = float(fields[5])
        if rtimeTmp != 0.0:
          rtime = rtimeTmp
        pepmz = precMzFromPrecMass(pepmass, charge)
        if not hasEZ:
          ezLines = list()
        hasEZ = True
        ezLines.append(Precursor(rtime, pepmz, charge, intensity))
      elif line.startswith('I') or line.startswith('H'):
        continue
      elif storeFragmentIons:
        fragmentIons.append(map(float, line.split()))
    
    if scannr != -1:
      yield Spectrum(scannr, ezLines, fragmentIons)
