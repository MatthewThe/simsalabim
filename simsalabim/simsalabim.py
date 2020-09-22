from __future__ import print_function

__version__ = "0.1.0"
__copyright__ = '''Copyright (c) 2020 Matthew The. All rights reserved.
Written by Matthew The (matthew.the@tum.de) in the
Chair of Proteomics and Bioanalytics at the 
Technical University of Munich.'''

import sys
import os

from collections import namedtuple

from pyteomics.mzml import MzML as mzMLParser

#######################
## input: MS/MS file ##
#######################

PrecursorBase = namedtuple("PrecursorBase", "rtime precmz charge intensity")
SpectrumBase = namedtuple("SpectrumBase", "scannr precursors fragmentIons")

class Precursor(PrecursorBase):
  def precMass(self):
    return precMassFromPrecMz(self.precmz, self.charge)
  
class Spectrum(SpectrumBase):  
  def toString(self, outputFormat = 'ms2'):
    if outputFormat == 'ms2':
      specString = 'S\t%d\t%d\t%f\n' % (self.scannr, self.scannr, self.precursors[0].precmz)
      for prec in self.precursors:
        specString += 'Z\t%d\t%f\n' % (prec.charge, prec.precMass())
      for prec in self.precursors:
        specString += 'I\tEZ\t%d\t%f\t%f\t%f\n' % (prec.charge, prec.precMass(), prec.rtime, prec.intensity)
      for mz, fInt in self.fragmentIons:
        specString += '%g %g\n' % (mz, fInt)
    elif outputFormat == 'mgf':
      for precursorIdx, prec in enumerate(self.precursors):
        specString = 'BEGIN IONS\n'
        specString += 'TITLE=scan=%d precursorIdx=%d\n' % (self.scannr, precursorIdx)
        specString += 'PEPMASS=%f %f\n' % (prec.precMass(), prec.intensity)
        specString += 'RTINSECONDS=%f\n' % (prec.rtime)
        specString += 'CHARGE=%d+\n' % (prec.charge)
        specString += 'SCANS=%d\n' % (self.scannr*100 + precursorIdx)
        for mz, fInt in self.fragmentIons:
          specString += '%g %g\n' % (mz, fInt)
        specString += 'END IONS\n'
    return specString
    
###########################################################
## input: MGF file                                       ##
##                                                       ##
## http://www.matrixscience.com/help/data_file_help.html ##
###########################################################

def parseMgfHeaders(inputFile):
  with open(inputFile, 'rb') as f:
    for line in f:
      if not line.startswith('BEGIN IONS'):
        yield line
      else:
        break

def parseMgf(inputFile, storeFragmentIons = False):
  if not os.path.isfile(inputFile):
    sys.exit("ERROR: Could not open mgf file: " + inputFile)
  
  if not inputFile.lower().endswith(".mgf"):
    sys.exit("ERROR: File does not have the .mgf file extension: " + inputFile)
  
  scannr = -1
  fragmentIons = list() if storeFragmentIons else None
  with open(inputFile, 'rb') as f:
    for line in f:
      if line.startswith('TITLE'):
        fields = line.split('=')
        if len(fields) > 2:
          scannr = int(line.split('=')[2].split("\"")[0])
        else:
          scannr = int(line.split('=')[1].split("[")[0])
      elif line.startswith('SCANS'):
        scannr = int(line.split('=')[1])
      elif line.startswith('RTINSECONDS'):
        rtime = float(line.split('=')[1])
      elif line.startswith('PEPMASS'):
        fields = line.split('=')[1].split()
        pepmass = float(fields[0])
        pepmz = precMzFromPrecMass(pepmass, charge)
        intensity = 0.0 if len(fields) > 1 else float(fields[1])
      elif line.startswith('CHARGE'):
        charge = int(line.split('=')[1][0])
      elif line[0].isdigit():
        fragmentIons.append(line.split())
      elif line.startswith('END IONS'):
        if scannr != -1:
          yield Spectrum(scannr, [Precursor(rtime, pepmz, charge, intensity)], fragmentIons)
        fragmentIons = list() if storeFragmentIons else None

###############################################################
## input: MS2 file                                           ##
##                                                           ##
## http://crux.ms/file-formats/ms2-format.html (base format) ##
## http://crux.ms/commands/bullseye.html (with EZ lines)     ##
###############################################################
    
def parseMs2Headers(inputFile):
  with open(inputFile, 'rb') as f:
    for line in f:
      if line.startswith('H'):
        yield line
      else:
        break

def generateMs2Headers(scriptName, authorInfo):
  from datetime import datetime
  
  headerLines = "H\tCreationDate\t%s\n" % str(datetime.now())
  headerLines += "H\tExtractor\t%s\n" % scriptName
  headerLines += "H\tExtractorVersion\t1.0\n"
  headerLines += "H\tComments\t%s written by %s\n" % (scriptName, authorInfo)
  headerLines += "H\tExtractorOptions\tN/A\n"
  return headerLines

def parseMs2(inputFile, storeFragmentIons = False):  
  if not os.path.isfile(inputFile):
    sys.exit("ERROR: Could not open ms2 file: " + inputFile)
  
  if not inputFile.lower().endswith(".ms2"):
    sys.exit("ERROR: File does not have the .ms2 file extension: " + inputFile)
  
  scannr = -1
  ezLines = list()
  fragmentIons = list() if storeFragmentIons else None
  hasEZ = False
  with open(inputFile, 'rb') as f:
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
          intensity = 0.0
          ezLines.append(Precursor(rtime, pepmz, charge, intensity))
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

####################################################
## input: mzML file                               ##
##                                                ##
## http://www.peptideatlas.org/tmp/mzML1.1.0.html ##
####################################################

# this function was loosely adapted from MzMLTransformer:
# https://github.com/mobiusklein/psims/blob/master/psims/transform/mzml.py
def parseMzML(inputFile, storeFragmentIons = False):
  parser = mzMLParser(inputFile, iterative = True)
  for i, spectrum in enumerate(parser.iterfind("spectrum")):
    if spectrum["ms level"] != 2:
      continue
    
    if storeFragmentIons:
      fragmentIons = zip(spectrum["m/z array"], spectrum["intensity array"])
    else:
      fragmentIons = None
    
    if "scan=" in spectrum["id"]:
      scannr = int(spectrum["id"].split("scan=")[1])
    else:
      scannr = i + 1
    
    scans = spectrum.get("scanList", {}).get('scan', [{}])[0]
    for key, value in list(scans.items()):
      if not hasattr(key, 'accession'):
        continue
      accession = key.accession
      if accession == "MS:1000016":
        rtime = value
    
    ezLines = list()
    precursors = spectrum.get("precursorList", {}).get("precursor")
    if precursors:
      for prec in precursors:
        ion = prec['selectedIonList'].get("selectedIon")[0]
        for key, value in list(ion.items()):
          if key == "selected ion m/z":
            pepmz = value
          elif key == "peak intensity":
            intensity = value
          elif key in ("charge state", "possible charge state"):
            charge = value
        
        #for key, value in prec.get("isolationWindow", {}).items():
        #  if key == "isolation window target m/z":
        #    isolationWindowTarget = value
        #  elif key == "isolation window lower offset":
        #    isolationWindowLowerOffset = value
        #  elif key == "isolation window upper offset":
        #    isolationWindowUpperOffset = value
        
        ezLines.append(Precursor(rtime, pepmz, charge, intensity))
    
    yield Spectrum(scannr, ezLines, fragmentIons)

######################
## Helper functions ##
######################

def precMzFromPrecMass(pmass, z):
  return (float(pmass) + 1.00727646677 * (int(z) - 1)) / int(z)

def precMassFromPrecMz(pmz, z):
  return pmz * z - 1.00727646677 * (z - 1)

