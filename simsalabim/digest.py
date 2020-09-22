from __future__ import print_function

import sys
import csv
import itertools
import collections

def main(argv):
  args = parseArgs()
  
  writePeptideToProteinMap(args)

def parseArgs():
  import argparse
  apars = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  addArguments(apars)
  
  args = apars.parse_args()
  
  return args
  
def addArguments(apars):  
  apars.add_argument('-F', '--fasta_file', default = None, metavar='F', required = True,
                     help='''Fasta database used in the search 
                             against the spectra file.
                          ''')
  
  apars.add_argument('-o', '--output_file', default = None, metavar='pin.tab', required = True,
                     help='''Save peptide to protein map to this tab delimited file.
                          ''')
  
  apars.add_argument('-P', '--pattern', default = "REV__", metavar='P',
                     help='''Pattern used to identify the decoy PSMs.
                          ''')
                            
  apars.add_argument('-e', '--enzyme', default = "trypsin", metavar='E',
                     help='''Type of enzyme "no_enzyme", "elastase","pepsin",
                             "proteinasek","thermolysin","chymotrypsin",
                             "lys-n","lys-c","arg-c","asp-n","glu-c","trypsin",
                             "trypsinp".
                          ''')    
                            
  apars.add_argument('-c', '--cleavages', default = 2, metavar='C', type=int,
                     help='''Number of allowed miss cleavages used in the search 
                             engine.
                          ''')
  
  apars.add_argument('-l', '--min-length', default = 7, metavar='L', type=int,
                     help='''Minimum peptide length allowed used in the search 
                             engine.
                          ''')
  
  apars.add_argument('-t', '--max-length', default = 60, metavar='L', type=int,
                     help='''Maximum peptide length allowed used in the search 
                             engine.
                          ''')                   
  
  apars.add_argument('--special-aas', default = 'KR', metavar='S', 
                     help='''Special AAs that MaxQuant uses for decoy generation.
                          ''')
  
def writePeptideToProteinMap(args):  
  writer = getTsvWriter(args.output_file)
  pre, not_post = getCleavageSites(args.enzyme)
  for peptide, proteins in getPeptideToProteinMap(args.fasta_file, min_len = args.min_length, max_len = args.max_length, pre = pre, not_post = not_post, miscleavages = args.cleavages, decoyPattern = args.pattern).items():
    writer.writerow([peptide, ";".join(proteins)])

def readFasta(filePath, db = "target", parseId = lambda x : x.split(" ")[0], specialAAs = ['K', 'R'], decoyPattern = "REV__"):
  if db not in ["target", "decoy", "concat"]:
    sys.exit("unknown db mode: %s" % db)
  
  hasSpecialAAs = len(specialAAs) > 0
  name, seq = None, []
  with open(filePath, 'r') as fp:
    for line in itertools.chain(fp, [">"]):
      line = line.rstrip()
      if line.startswith(">"):
        if name: 
          if db in ["target", "concat"]:
            yield (name, "".join(seq))
          
          if db in ["decoy", "concat"]:
            seq = list("".join(seq)[::-1])
            if hasSpecialAAs:
              for i in range(1, len(seq)):
                if seq[i] in specialAAs:
                  swapPositions(seq, i, i-1)
            yield (decoyPattern + name, "".join(seq))
          
        if len(line) > 1:
          name, seq = parseId(line[1:]), []
      else: seq.append(line)

def swapPositions(seq, pos1, pos2):     
  seq[pos1], seq[pos2] = seq[pos2], seq[pos1] 

# method for generating the list of peptides
def getDigestedPeptides(seq, min_len = 6, max_len = 50, pre = ['K', 'R'], not_post = ['P'], digestion = 'full', miscleavages = 0, methionineCleavage = True):
  lenS, starts = len(seq), [0]
  length_accepted = lambda x : x >= min_len and x <= max_len
  methionineCleavage = methionineCleavage and seq[0] == "M"
  
  if digestion == 'none':
    for i in range(lenS + 1):
      for j in range(i + min_len, i + max_len + 1):
        if j <= lenS:
          yield seq[i:j]
  elif digestion == 'semi':
    for i in range(lenS + 1):
      isCleavageSite = (seq[min([lenS-1,i])] in pre and seq[min([lenS-1,i+1])] not in not_post)
      isMethionineCleavageSite = (i == 0 and methionineCleavage)
      if i == lenS or isCleavageSite or isMethionineCleavageSite:
        # peptides with enzymatic C-terminal (both enzymatic and non-enzymatic N-terminal)
        start = starts[0]
        for j in range(start, min([i+1, lenS])):
          lenP = min([i, lenS - 1]) - j + 1
          if length_accepted(lenP):
            yield (seq[j : i + 1])
        starts.append(i + 1)
        methionineCleaved = int(methionineCleavage and starts[0] == 0)
        if len(starts) > miscleavages + 1 + methionineCleaved or i == lenS:    
          starts = starts[1 + methionineCleaved:]
      else: # peptides with non enzymatic C-terminal
        for start in starts:
          lenP = i - start + 1
          if length_accepted(lenP) and i + 1 not in starts:
            yield (seq[start : i + 1])
  else:
    hasMethionineCleavage = False
    cleavageSites = []
    if methionineCleavage and seq[0] == 'M':
      cleavageSites.append(0)
      hasMethionineCleavage = True
    cleavageSites.extend([i for i in range(lenS) if seq[i] in pre and seq[min([lenS-1,i+1])] not in not_post])
    cleavageSites.append(lenS)
    for i in cleavageSites:
      for start in starts:
        lenP = i - start + 1
        if length_accepted(lenP):
          yield (seq[start : i + 1])
      starts.append(i + 1)
      methionineCleaved = int(starts[0] == 0 and hasMethionineCleavage)
      if len(starts) > miscleavages + 1 + methionineCleaved:    
        starts = starts[1 + methionineCleaved:]

def getPeptideToProteinMap(fastaFile, db = "concat", min_len = 6, max_len = 52, pre = ['K', 'R'], not_post = ['P'], digestion = 'full', miscleavages = 2, methionineCleavage = True, specialAAs = ['K', 'R'], decoyPattern = "REV__"):
  peptideToProteinMap = collections.defaultdict(list)
  proteinToSeqMap = dict()
  #parseId = lambda x : x.replace(" ", "") # default MaxQuant 
  #parseId = lambda x : x.split("|")[1] # extract Uniprot ID, for PrDB runs (mouse proteome, Dongxue original)
  parseId = lambda x : x.split(" ")[0] # for Dongxue mimic
  for proteinIdx, (protein, seq) in enumerate(readFasta(fastaFile, db, parseId, specialAAs = specialAAs, decoyPattern = decoyPattern)):
    if (proteinIdx+1) % 10000 == 0:
      print("Digesting protein", proteinIdx+1)
    seenPeptides = set()
    proteinToSeqMap[protein] = seq
    for peptide in getDigestedPeptides(seq, min_len, max_len, pre, not_post, digestion, miscleavages, methionineCleavage):
      if peptide not in seenPeptides:
        seenPeptides.add(peptide)
        peptideToProteinMap[peptide].append(protein)
  
  return peptideToProteinMap

def getCleavageSites(enzyme):
  if enzyme == "trypsinp":
    pre = ['K', 'R']
    not_post = []
  elif enzyme == "trypsin":
    pre = ['K', 'R']
    not_post = ['P']
  elif enzyme == "no_enzyme":
    pre = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    not_post = []
  elif enzyme == "chymotrypsin":
    pre = ['F', 'W', 'Y', 'L']
    not_post = ['P']
  elif enzyme == "proteinasek":
    pre = ['A', 'E', 'F', 'I', 'L', 'T', 'V', 'W', 'Y']
    not_post = []
  elif enzyme == "elastase":
    pre = ['L', 'V', 'A', 'G']
    not_post = ['P']
  elif enzyme == "lys-c":
    pre = ['K']
    not_post = ['P']
  elif enzyme == "arg-c":
    pre = ['R']
    not_post = ['P']
  elif enzyme == "glu-c":
    pre = ['E']
    not_post = ['P']
  elif enzyme == 'v8-de':
    pre = ['N', 'D', 'E', 'Q']
    not_post = ['P']
  else:
    sys.exit("Enzyme", enzyme, "not implemented yet")
  
  return pre, not_post

def getTsvWriter(filename, delimiter = '\t'):
  # Python 3
  if sys.version_info[0] >= 3:
    return csv.writer(open(filename, 'w', newline = ''), delimiter = delimiter)
  # Python 2
  else:
    return csv.writer(open(filename, 'wb'), delimiter = delimiter)
    
if __name__ == "__main__":
  main(sys.argv[1:])
