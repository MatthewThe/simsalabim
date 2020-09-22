from __future__ import print_function

import sys
import os
import subprocess

from .simsalabim import __version__, __copyright__
from . import add_quant_info as quant
from . import helpers

def main(argv):
  print('dinosaur-adapter version %s\n%s' % (__version__, __copyright__))
  print('Issued command:', os.path.basename(__file__) + " " + " ".join(map(str, sys.argv[1:])))
  
  args, params = parseArgs()
  
  run_dinosaur(args.dinosaur_jar_path, args.mzml_fns, args.output_folder, args.spectrum_output_format, params)

def parseArgs():
  import argparse
  apars = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  requiredNamed = apars.add_argument_group('required arguments')
  
  requiredNamed.add_argument('--dinosaur_jar_path', metavar = "JAR", required = True,
                     help='''Path to the Dinosaur .jar file.
                          ''')
                          
  apars.add_argument('--mzml_fns', default=None, metavar = "M", nargs='*',
                     help='''mzML file(s). To easily specify multiple files one can use wildcards, e.g. my_spectrum_files/*.mzML
                          ''')
  
  apars.add_argument('--file_list_file', default=None, metavar = "L",
                     help='''Text file with paths to mzML files, one per line.
                          ''')
  
  apars.add_argument('--output_folder', default="./dinosaur/", metavar='O', 
                     help='''Output folder.
                          ''')
  
  apars.add_argument('--dinosaur_mem', default=8.0, metavar='M', type=float,
                     help='''Memory for allocated for Dinosaur in GB.
                          ''')
  
  apars.add_argument('--dinosaur_flags', default="", metavar='O', 
                     help='''Extra command line flags to pass to Dinosaur, as indicated in Dinosaur's help text.
                          ''')                        
                                                  
  apars.add_argument('--spectrum_output_format', default=None, metavar='F', 
                     help='''If you want updated spectrum files with the new MS1 features assigned to the MS2 spectra, set this to the desired output format (ms2, mgf or mzML).
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
  
  # ------------------------------------------------
  args = apars.parse_args()
  
  if not args.mzml_fns:
    if args.file_list_file and len(args.file_list_file) > 0:
      with open(args.file_list_file, 'r') as f:
        args.mzml_fns = list(filter(lambda x : len(x) > 0, map(lambda x : re.sub(r"[\n\t\s]*", "", x), f.read().splitlines())))
    else:
      sys.exit("No input mzML files specified. Use either --mzml_fns or --file_list_file.")
  elif args.file_list_file and len(args.file_list_file) > 0:
    sys.exit("Ambiguous mzML input. Use either --mzml_fns or --file_list_file, not both.")
  
  params = dict()
  params['splitPrecursors'] = args.split_precursors
  params['dinosaurMemory'] = args.dinosaur_mem
  params['dinosaurFlags'] = args.dinosaur_flags
  
  return args, params

def run_dinosaur(dinosaur_jar_path, mzml_fns, output_folder, spectrum_output_format, params):
  dinosaur_binary = "java -Xmx%dM -jar %s --seed=1" % (int(params['dinosaurMemory']*1000), dinosaur_jar_path)
  helpers.createDir(output_folder)
  for mzml_fn in mzml_fns:
    baseFN = helpers.getBase(helpers.getFileName(mzml_fn))
    dinosaur_output_file = os.path.join(output_folder, baseFN + ".features.tsv")
    if not os.path.isfile(dinosaur_output_file):
      cmd_dinosaur = "%s --force --outDir=%s %s %s;" % (dinosaur_binary, output_folder, params['dinosaurFlags'], mzml_fn)
      helpers.executeCmd(cmd_dinosaur)
    else:
      print("Found dinosaur output file at %s, remove this file to re-run Dinosaur on this file" % (dinosaur_output_file))
    
    output_fn = os.path.join(output_folder, baseFN + ".dummy.txt")
    if spectrum_output_format:
      output_fn = os.path.join(output_folder, baseFN + ".recalibrated." + spectrum_output_format)
    
    params['specPrecMapFile'] = os.path.join(output_folder, baseFN + ".feature_map.tsv")
    if not os.path.isfile(params['specPrecMapFile']):
      quant.add_accurate_precursors(dinosaur_output_file, mzml_fn, output_fn, params)
      if output_fn.endswith(".dummy.txt"):
        os.remove(output_fn)
    else:
      print("Found dinosaur mapping file at %s, remove this file to re-run Dinosaur on this file" % (params['specPrecMapFile']))
 
    
if __name__ == '__main__':
  main(sys.argv[1:])
