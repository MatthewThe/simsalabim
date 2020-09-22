from __future__ import print_function

import sys
import os
import csv
import copy
import re

import numpy as np
from psims.transform.mzml import MzMLTransformer

from .simsalabim import __version__, __copyright__
from . import helpers

def main(argv):
  print('add-quant-info version %s\n%s' % (__version__, __copyright__))
  print('Issued command:', os.path.basename(__file__) + " " + " ".join(map(str, sys.argv[1:])))
  
  args, params = parseArgs()
  
  add_accurate_precursors(args.feature_fn, args.mzml_fn, args.output_fn, params)

def parseArgs():
  import argparse
  apars = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  requiredNamed = apars.add_argument_group('required arguments')
  
  apars.add_argument('mzml_fn', default=None, metavar = "IN_FILE",
                     help='''mzML file.
                          ''')
  
  requiredNamed.add_argument('--feature_fn', metavar='F',
                     help='feature csv table from Dinosaur or OpenMS',
                     required = True)
  
  apars.add_argument('--output_fn', default = "spectra.mzML", metavar='OUT', 
                     help='''output file in ms2, mgf or mzML format (optional).
                          ''')
  
  apars.add_argument('--split_precursors',
                     help='''for .mzML or .ms2 output this creates a new spectrum for each precursor, e.g. 
                             if spectrum with scan number 132 matches two precursors, we generate two spectra 
                             with scan numbers 13201 and 13202. This can be useful if your downstream 
                             analysis includes tools that do not support multiple precursors per spectrum,
                             such as MSGF+. For MGF output this flag is always set, as it does not support 
                             multiple precursors per spectrum.
                          ''',
                     action='store_true')
                     
  apars.add_argument('--map_fn', default = "", metavar='M', 
                     help='''tab separated file containing a mapping from spectra to precursor information.
                          ''')
  
  # ------------------------------------------------
  args = apars.parse_args()
  
  params = dict()
  params['specPrecMapFile'] = args.map_fn
  params['splitPrecursors'] = args.split_precursors or args.mzml_fn.lower().endswith(".mgf")
  
  return args, params
  
def load_feature_table(fn):
  table = []
  with helpers.openVersionSafe(fn, 'r') as fh:
    rd = csv.reader(fh, delimiter=',')
    for row in rd:
      if row[0] == 'FEATURE':
        _, rt, mz, _, chg, _, _, _, _, rtl, rtr = row
        table.append([float(mz), int(chg), float(rtl), float(rtr), float(rt)])
  table.sort(key=lambda x: x[3])
  return table

def load_dinosaur_table(fn, minCharge=1):
  table = []
  with helpers.openVersionSafe(fn, 'r') as fh:
    rd = csv.reader(fh, delimiter='\t')
    for row in rd:
      if row[0] != 'mz':
        mz, _, chg, rtl, rt, rtr, _, _, _, _, _, _, int_apex, int_sum = row
        z = int(chg)
        if z >= minCharge:
          table.append([float(mz), z, float(rtl)*60, float(rtr)*60, float(rt)*60, float(int_sum)])
  table.sort(key=lambda x: x[3])
  return table

def add_accurate_precursors(feature_fn, mzml_fn, ms2_outpath, params):
  if feature_fn[-13:] == ".features.tsv" or feature_fn[-13:] == ".features.csv":
    features = load_dinosaur_table(feature_fn)
  else:
    features = load_feature_table(feature_fn)
  features.sort(key=lambda x: x[0])
  fmz_all = np.array([f[0] for f in features])

  if ms2_outpath and not ms2_outpath.endswith(".dummy.txt") and not (ms2_outpath.lower().endswith(".ms2") or ms2_outpath.lower().endswith(".mzml") or ms2_outpath.lower().endswith(".mgf")):
    sys.exit("ERROR: Could not detect output format from filename. Please use the extension .mzML, .ms2 or .mgf for your output file")  
  
  transformSpectra(features, fmz_all, mzml_fn, ms2_outpath, params)
      
def transformSpectra(features, fmz_all, mzml_fn, ms2_outpath, params):
  transform = replace_precursors_multi_out if params['splitPrecursors'] else replace_precursors
  if len(params['specPrecMapFile']) > 0:
    specPrecMapWriter = helpers.openVersionSafe(params['specPrecMapFile'], 'w')
    specPrecMapWriter.write('FileName\tScanNr\tPrecMz\tCharge\tRTime\tIntensity\n')
  else:
    specPrecMapWriter = None
  
  mzml_fn_base = os.path.basename(mzml_fn)
  
  writeMode = 'w'
  if ms2_outpath.lower().endswith(".mzml"):
    writeMode = 'wb'
  
  with helpers.openVersionSafe(ms2_outpath, writeMode) as out:
    trans = MzMLTransformerMultiOut(mzml_fn, out, transform = lambda x : transform(x, features, fmz_all, specPrecMapWriter, mzml_fn_base), transform_description = "Assigned accurate precursor information from Dinosaur")
    if ms2_outpath.lower().endswith(".mzml"):
      trans.write()
    elif ms2_outpath.lower().endswith(".ms2"):
      trans.write_ms2()
    elif ms2_outpath.lower().endswith(".mgf"):
      trans.write_mgf()
    else:
      trans.transform_only()
  
  print("Done")

  
def replace_precursors(spectrum, features, fmz_all, specPrecMapWriter, mzml_fn_base):
  if spectrum["ms level"] == 2:
    rt = get_rtime(spectrum)
    precursors = spectrum.get("precursorList", {}).get("precursor")
    returnSpectrum = False
    if precursors:
      new_precursors = list()
      for prec in precursors:
        pmz, iso_width_lower, iso_width_upper = get_isolation_window(prec)
        candidates = get_candidate_precursors(features, fmz_all, pmz, iso_width_lower, iso_width_upper, rt)
        if len(candidates) > 0:
          returnSpectrum = True
          for f in candidates:
            prec_copy = update_precursor(prec, f)
            new_precursors.append(prec_copy)
            
            writeSpecPrecRow(specPrecMapWriter, spectrum["id"], f, mzml_fn_base)
          
          spectrum["precursorList"]["precursor"] = new_precursors
    
    if returnSpectrum:
      return [spectrum]
    else:
      return []
  else:
    return [spectrum]

def replace_precursors_multi_out(spectrum, features, fmz_all, specPrecMapWriter, mzml_fn_base):
  if spectrum["ms level"] == 2:
    rt = get_rtime(spectrum)
    precursors = spectrum.get("precursorList", {}).get("precursor")
    if precursors:
      spectra = list()
      for prec in precursors:
        pmz, iso_width_lower, iso_width_upper = get_isolation_window(prec)
        candidates = get_candidate_precursors(features, fmz_all, pmz, iso_width_lower, iso_width_upper, rt)
        if len(candidates) > 0:
          for precIdx, f in enumerate(candidates[:99]):
            fmz, fz, frt_left, frt_right, frt, fint = f
            spectrum_copy = spectrum.copy()
            spectrum_copy["precursorList"] = spectrum["precursorList"].copy()
            
            originalScanNr = getScanNr(spectrum_copy["id"])
            newScanNr = originalScanNr*100+precIdx+1
            spectrum_copy["id"] = replaceScanNr(spectrum_copy["id"], str(newScanNr) + " originalScan=" + str(originalScanNr))
            #spectrum_copy["id"] += " precursorIdx=" + str(precIdx+1)
            prec_copy = update_precursor(prec, f)
            spectrum_copy["precursorList"]["precursor"] = [prec_copy]
            spectra.append(spectrum_copy)
            
            writeSpecPrecRow(specPrecMapWriter, spectrum_copy["id"], f, mzml_fn_base)
      return spectra
    else:
      return []
  else:
    return [spectrum]

def update_precursor(prec, new_prec_info):
  fmz, fz, _, _, _, fint = new_prec_info
  prec_copy = copy.deepcopy(prec)
  ion = prec_copy['selectedIonList'].get("selectedIon")[0]
  for key, value in list(ion.items()):
    if key == "selected ion m/z":
      ion[key] = fmz
    elif key == "peak intensity":
      ion[key] = fint
    elif key in ("charge state", "possible charge state"):
      ion[key] = fz
  return prec_copy
  
def get_isolation_window(prec):
  ion = prec['selectedIonList'].get("selectedIon")[0]
  for key, value in prec.get("isolationWindow", {}).items():
    if key == "isolation window target m/z":
      isolationWindowTarget = float(value)
    elif key == "isolation window lower offset":
      isolationWindowLowerOffset = float(value)
    elif key == "isolation window upper offset":
      isolationWindowUpperOffset = float(value)
  return isolationWindowTarget, isolationWindowLowerOffset, isolationWindowUpperOffset

def get_rtime(spectrum):
  scans = spectrum.get("scanList", {}).get('scan', [{}])[0]
  rtime = 0.0
  timescale = 1
  for key, value in list(scans.items()):
    if not hasattr(key, 'accession'):
      continue
    accession = key.accession
    if accession == "MS:1000016": # scan_start_time
      rtime = value
      if key.unit_accession == 'UO:0000031': # UO:0000031 = minute; UO:0000010 = second
        timescale = 60
  
  return rtime * timescale

def get_candidate_precursors(features, fmz_all, pmz, iso_width_lower, iso_width_upper, rt):
  l_idx = fmz_all.searchsorted(pmz - iso_width_lower, side='left')
  r_idx = fmz_all.searchsorted(pmz + iso_width_upper, side='right')
  candidates = list()
  for f in features[l_idx:r_idx]:
    fmz, fz, frt_left, frt_right, frt, fint = f
    if frt_left < rt < frt_right:
      candidates.append(f)
  return candidates

def getScanNr(specId):
  return int(specId.split('scan=')[-1].split()[0])

def replaceScanNr(specId, newScanNr):
  return re.sub(r'scan=([0-9]*)', 'scan=%s' % newScanNr, specId)

def writeSpecPrecRow(specPrecMapWriter, specId, f, mzml_fn_base):
  if specPrecMapWriter:
    fmz, fz, _, _, frt, fint = f
    specPrecMapWriter.write("%s\t%s\t%f\t%d\t%f\t%f\n" % (mzml_fn_base, specId, fmz, fz, frt, fint))
  
class MzMLTransformerMultiOut(MzMLTransformer):
  def write(self):
    '''Write out the the transformed mzML file
    '''
    writer = self.writer
    with writer:
      writer.controlled_vocabularies()
      self.copy_metadata()
      with writer.run(id="transformation_run"):
        with writer.spectrum_list(len(self.reader._offset_index)):
          self.reader.reset()
          for i, spectrum in enumerate(self.iterspectrum()):
            spectra = self.transform(spectrum)
            for spectrum in spectra:
              self.writer.write_spectrum(**self.format_spectrum(spectrum))
            if i % 1000 == 0:
              self.log("Handled %d spectra" % (i, ))
  
  def transform_only(self):
    for i, spectrum in enumerate(self.iterspectrum()):
      spectra = self.transform(spectrum)
      if i % 1000 == 0:
        self.log("Handled %d spectra" % (i, ))
        
  def write_ms2(self):
    writer = Ms2Writer(self.output_stream)
    writer.write_ms2_headers()
    for i, spectrum in enumerate(self.iterspectrum()):
      spectra = self.transform(spectrum)
      for spectrum in spectra:
        if spectrum['ms level'] == 2:
          writer.write_ms2_spectrum(**self.format_spectrum(spectrum))
      if i % 1000 == 0:
        self.log("Handled %d spectra" % (i, ))
  
  def write_mgf(self):
    writer = MgfWriter(self.output_stream)
    for i, spectrum in enumerate(self.iterspectrum()):
      spectra = self.transform(spectrum)
      for spectrum in spectra:
        if spectrum['ms level'] == 2:
          writer.write_mgf_spectrum(**self.format_spectrum(spectrum))
      if i % 1000 == 0:
        self.log("Handled %d spectra" % (i, ))

class Ms2Writer:
  def __init__(self, output_stream):
    self.ms2_writer = output_stream
      
  def write_ms2_spectrum(self, mz_array=None, intensity_array=None, charge_array=None, id=None,
       polarity='positive scan', centroided=True, precursor_information=None,
       scan_start_time=None, params=None, compression=None,
       encoding=None, other_arrays=None, scan_params=None, scan_window_list=None,
       instrument_configuration_id=None, intensity_unit=None):
    if precursor_information:
      # in the MS2 format retention times are in minutes
      timescale = 1
      if scan_start_time['unit_name'] == "second":
        timescale = 60
      
      written_scan_number = False
      for p in precursor_information:
        if not written_scan_number:
          idx = getScanNr(id)
          self.ms2_writer.write("S\t%d\t%d\t%f\n" % (idx, idx, p['mz']))
          self.ms2_writer.write("I\tRTime\t%f\n" % (scan_start_time['value'] / timescale))
          written_scan_number = True
        
        fmass = helpers.precMassFromPrecMz(p['mz'], p['charge'])
        if p['intensity']:
          self.ms2_writer.write("I\tEZ\t%d\t%f\t%f\t%f\n" % (p['charge'], fmass, scan_start_time['value'] / timescale, p['intensity']))
        self.ms2_writer.write("Z\t%d\t%f\n" % (p['charge'], fmass))
      
      if written_scan_number:
        for a, b in zip(mz_array, intensity_array):
          self.ms2_writer.write("%f %f\n" % (a, b))
    
  def write_ms2_headers(self):
    from datetime import datetime
    
    self.ms2_writer.write("H\tCreationDate\t%s\n" % str(datetime.now()))
    self.ms2_writer.write("H\tExtractor\t%s\n" % os.path.basename(__file__))
    self.ms2_writer.write("H\tExtractorVersion\t1.0\n")
    self.ms2_writer.write("H\tComments\t%s written by Matthew The, 2019\n" % os.path.basename(__file__))
    self.ms2_writer.write("H\tExtractorOptions\tN/A\n")

class MgfWriter:
  def __init__(self, output_stream):
    self.mgf_writer = output_stream
      
  def write_mgf_spectrum(self, mz_array=None, intensity_array=None, charge_array=None, id=None,
       polarity='positive scan', centroided=True, precursor_information=None,
       scan_start_time=None, params=None, compression=None,
       encoding=None, other_arrays=None, scan_params=None, scan_window_list=None,
       instrument_configuration_id=None, intensity_unit=None):
    if precursor_information:
      if len(precursor_information) != 1:
        sys.exit("ERROR: cannot have 0 or multiple precursors for one spectrum in MGF format.")
      
      timescale = 1
      if scan_start_time['unit_name'] == "minute":
        timescale = 60
      
      p = precursor_information[0]
      fmass = helpers.precMassFromPrecMz(p['mz'], p['charge'])
      idx = getScanNr(id)
      self.mgf_writer.write("BEGIN IONS\n")
      self.mgf_writer.write("TITLE=%s\n" % id)
      self.mgf_writer.write("PEPMASS=%f %f\n" % (fmass, p['intensity']))
      self.mgf_writer.write("RTINSECONDS=%f\n" % (scan_start_time['value']*timescale))
      self.mgf_writer.write("CHARGE=%d+\n" % (p['charge']))
      self.mgf_writer.write("SCANS=%d\n" % (idx))
      
      for a, b in zip(mz_array, intensity_array):
        self.mgf_writer.write("%f %f\n" % (a, b))
      
      self.mgf_writer.write("END IONS\n")
  
if __name__ == '__main__':
  main(sys.argv[1:])
