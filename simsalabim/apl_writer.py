import sys
from pathlib import Path

from . import helpers


class AplWriter:
    def __init__(self, output_stream):
        self.apl_writer = output_stream
    
    def write(self, spectra):
        for i, spectrum in enumerate(spectra):
            if spectrum.get_ms_level() == 2:
                self.write_apl_spectrum(**spectrum.spectrum, raw_file=spectrum.raw_file)
            if i % 1000 == 0:
                print("Handled %d spectra" % (i, ))

    def write_apl_spectrum(self, mz_array=None, intensity_array=None, charge_array=None, id=None,
             polarity='positive scan', centroided=True, precursor_information=None,
             scan_start_time=None, params=None, compression=None,
             encoding=None, other_arrays=None, scan_params=None, scan_window_list=None,
             instrument_configuration_id=None, intensity_unit=None, raw_file=None):
        if precursor_information:
            if len(precursor_information) != 1:
                sys.exit("ERROR: cannot have 0 or multiple precursors for one spectrum in MGF format.")
            
            timescale = 1
            if scan_start_time and scan_start_time['unit_name'] == "minute":
                timescale = 60
            
            p = precursor_information[0]
            fmass = helpers.precMassFromPrecMz(p['mz'], p['charge'])
            idx = helpers.getScanNr(id)
            self.apl_writer.write("peaklist start\n")
            self.apl_writer.write(f"header=RawFile: {raw_file} Index: {idx}\n")  # e.g. header=RawFile: 01308_D03_P013387_B00_N20_R1 Index: 55569 Precursor: 0 _multi_
            self.apl_writer.write("mz=%f\n" % (p['mz']))
            self.apl_writer.write("charge=%d+\n" % (p['charge']))
            
            for a, b in zip(mz_array, intensity_array):
                self.apl_writer.write("%f %f\n" % (a, b))
            
            self.apl_writer.write("peaklist end\n\n")


class AplWriterMultiFile:
    def __init__(self, output_dir):
        self.output_dir = Path(output_dir)
        self.apl_writers = dict()
    
    def write(self, spectra):
        for i, spectrum in enumerate(spectra):
            if spectrum.get_ms_level() == 2:
                if spectrum.raw_file not in self.apl_writers:
                    output_file = self.output_dir / Path(spectrum.raw_file + ".apl")
                    output_stream = open(output_file, 'w')
                    self.apl_writers[spectrum.raw_file] = AplWriter(output_stream)
                self.apl_writers[spectrum.raw_file].write_apl_spectrum(**spectrum.spectrum)
            if i % 1000 == 0:
                print("Handled %d spectra" % (i, ))
