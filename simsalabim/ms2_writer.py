import os

from . import helpers


class Ms2Writer:
  def __init__(self, output_stream):
      self.ms2_writer = output_stream
  
  def write(self, spectra):
      self.write_ms2_headers()
      for i, spectrum in enumerate(spectra):
          if spectrum.get_ms_level() == 2:
              self.write_ms2_spectrum(**spectrum.spectrum)
          if i % 1000 == 0:
              print("Handled %d spectra" % (i, ))

  def write_ms2_spectrum(self, mz_array=None, intensity_array=None, charge_array=None, id=None,
       polarity='positive scan', centroided=True, precursor_information=None,
       scan_start_time=None, params=None, compression=None,
       encoding=None, other_arrays=None, scan_params=None, scan_window_list=None,
       instrument_configuration_id=None, intensity_unit=None):
    if precursor_information:
      # in the MS2 format retention times are in minutes
      timescale = 1
      if scan_start_time and scan_start_time['unit_name'] == "second":
        timescale = 60
      
      written_scan_number = False
      for p in precursor_information:
        if not written_scan_number:
          idx = helpers.getScanNr(id)
          self.ms2_writer.write("S\t%d\t%d\t%f\n" % (idx, idx, p['mz']))
          if scan_start_time:
              self.ms2_writer.write("I\tRTime\t%f\n" % (scan_start_time['value'] / timescale))
          written_scan_number = True
        
        fmass = helpers.precMassFromPrecMz(p['mz'], p['charge'])
        if p['intensity'] and scan_start_time:
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
