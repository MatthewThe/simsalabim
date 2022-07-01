from psims.transform.mzml import MzMLWriter

from .mgf_writer import MgfWriter
from .ms2_writer import Ms2Writer


class MzMLWriter2(MzMLWriter):
    def write(self, spectra):
        '''Write out the the transformed mzML file
        '''
        with self:
            self.controlled_vocabularies()
            #self.copy_metadata()
            with self.run(id="id"):
                for i, spectrum in enumerate(spectra):
                    self.write_spectrum(**spectrum.spectrum)
                    if i % 1000 == 0:
                        print("Handled %d spectra" % (i, ))

