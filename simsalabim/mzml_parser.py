from pyteomics.mzml import MzML as mzMLParser

from .spectrum import Spectrum
from .precursor import Precursor
from . import helpers

####################################################
## input: mzML file                               ##
##                                                ##
## http://www.peptideatlas.org/tmp/mzML1.1.0.html ##
####################################################

# this function was loosely adapted from MzMLTransformer:
# https://github.com/mobiusklein/psims/blob/master/psims/transform/mzml.py
def parse_mzml(inputFile, storeFragmentIons = False):
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
        
        ezLines.append({'mz': pepmz, 'charge': charge, 'intensity': intensity})
    
    yield Spectrum(scannr, ezLines, fragmentIons)
