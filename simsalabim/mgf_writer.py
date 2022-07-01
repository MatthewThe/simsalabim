import sys

from . import helpers


class MgfWriter:
    def __init__(self, output_stream):
        self.mgf_writer = output_stream
    
    def write(self, spectra):
        for i, spectrum in enumerate(spectra):
            if spectrum.get_ms_level() == 2:
                self.write_mgf_spectrum(**spectrum.spectrum)
            if i % 1000 == 0:
                print("Handled %d spectra" % (i, ))

    def write_mgf_spectrum(self, mz_array=None, intensity_array=None, charge_array=None, id=None,
             polarity='positive scan', centroided=True, precursor_information=None,
             scan_start_time=None, params=None, compression=None,
             encoding=None, other_arrays=None, scan_params=None, scan_window_list=None,
             instrument_configuration_id=None, intensity_unit=None):
        if precursor_information:
            if len(precursor_information) != 1:
                sys.exit("ERROR: cannot have 0 or multiple precursors for one spectrum in MGF format.")
            
            timescale = 1
            if scan_start_time and scan_start_time['unit_name'] == "minute":
                timescale = 60
            
            p = precursor_information[0]
            fmass = helpers.precMassFromPrecMz(p['mz'], p['charge'])
            idx = helpers.getScanNr(id)
            self.mgf_writer.write("BEGIN IONS\n")
            self.mgf_writer.write("TITLE=%s\n" % id)
            if p['intensity']:
                self.mgf_writer.write("PEPMASS=%f %f\n" % (fmass, p['intensity']))
            else:
                self.mgf_writer.write("PEPMASS=%f\n" % (fmass))
            if scan_start_time:
                self.mgf_writer.write("RTINSECONDS=%f\n" % (scan_start_time['value']*timescale))
            self.mgf_writer.write("CHARGE=%d+\n" % (p['charge']))
            self.mgf_writer.write("SCANS=%d\n" % (idx))
            
            for a, b in zip(mz_array, intensity_array):
                self.mgf_writer.write("%f %f\n" % (a, b))
            
            self.mgf_writer.write("END IONS\n")
