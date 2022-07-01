import sys
import os
from pathlib import Path

from .simsalabim import __version__, __copyright__

from .apl_parser import parse_apl
from .apl_writer import AplWriterMultiFile

from . import helpers

def main(argv):
    print('simsalabim-split-apl version %s\n%s' % (__version__, __copyright__))
    print('Issued command:', os.path.basename(__file__) + " " + " ".join(map(str, sys.argv[1:])))
    
    args, params = parseArgs()
    
    split_apl_by_raw_file(args.input_dir, args.output_dir)


def parseArgs():
    import argparse
    apars = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    requiredNamed = apars.add_argument_group('required arguments')
    
    apars.add_argument('input_dir', default=None, metavar = "IN",
                                         help='''directory with files in APL format.
                                                    ''')
    
    apars.add_argument('--output_dir', default=None, metavar='OUT', 
                                         help='''output directory.
                                                    ''')
    
    # ------------------------------------------------
    args = apars.parse_args()
    
    params = dict()
    
    return args, params


def split_apl_by_raw_file(input_dir, output_dir):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    spectra = get_all_spectra(input_dir)
    
    writer = AplWriterMultiFile(output_dir)
    writer.write(spectra)

    print("Done")


def get_all_spectra(input_dir):
    for input_file in Path(input_dir).glob('*.apl'):
        yield from parse_apl(str(input_file), storeFragmentIons=True)


if __name__ == '__main__':
    main(sys.argv[1:])
