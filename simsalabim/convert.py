import sys
import os

from .simsalabim import __version__, __copyright__

from .ms2_parser import parse_ms2
from .mgf_parser import parse_mgf
from .apl_parser import parse_apl
from .mzml_parser import parse_mzml

from .mzml_writer import MzMLWriter2
from .ms2_writer import Ms2Writer
from .mgf_writer import MgfWriter

from . import helpers

def main(argv):
    print('simsalabim-convert version %s\n%s' % (__version__, __copyright__))
    print('Issued command:', os.path.basename(__file__) + " " + " ".join(map(str, sys.argv[1:])))
    
    args, params = parseArgs()
    
    convert(args.input_file, args.output_fn, params)


def parseArgs():
    import argparse
    apars = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    requiredNamed = apars.add_argument_group('required arguments')
    
    apars.add_argument('input_file', default=None, metavar = "IN_FILE",
                                         help='''input MS/MS file, either ms2, mgf, apl or mzML format.
                                                    ''')
    
    apars.add_argument('--output_fn', default = "spectra.mzML", metavar='OUT', 
                                         help='''output file in ms2, mgf or mzML format.
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
    params['splitPrecursors'] = args.split_precursors or args.output_fn.lower().endswith(".mgf")
    
    return args, params


def convert(input_file, output_file, params):
    if input_file.lower().endswith(".mzml"):
        parser = parse_mzml
    elif input_file.lower().endswith(".ms2"):
        parser = parse_ms2
    elif input_file.lower().endswith(".mgf"):
        parser = parse_mgf
    elif input_file.lower().endswith(".apl"):
        parser = parse_apl
        
    spectra = parser(input_file, storeFragmentIons=True)
    
    writeMode = 'w'
    if output_file.lower().endswith(".mzml"):
        writeMode = 'wb'

    with helpers.openVersionSafe(output_file, writeMode) as out:
        if output_file.lower().endswith(".mzml"):
            writer = MzMLWriter2(out)
        elif output_file.lower().endswith(".ms2"):
            writer = Ms2Writer(out)
        elif output_file.lower().endswith(".mgf"):
            writer = MgfWriter(out)
    
        writer.write(spectra)

    print("Done")


if __name__ == '__main__':
    main(sys.argv[1:])
