from psims.transform.mzml import MzMLTransformer

from .mgf_writer import MgfWriter
from .ms2_writer import Ms2Writer


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

