import os
import subprocess

### command line helpers ###

def executeCmd(cmd):
  print(cmd)
  print("")
  sys.stdout.flush()
  rc = subprocess.call(cmd, shell=True)
  #rc = 0
  if rc == 1:
    print("Error while processing " + cmd)
    return 1
  else:
    return 0

### IO helpers ###

def openVersionSafe(filename, flag):
  # Python 3
  if sys.version_info[0] >= 3:
    if 'b' in flag:
      return open(filename, flag)
    else:
      return open(filename, flag, newline = '')
  # Python 2
  else:
    return open(filename, flag + 'b')
    
def createDir(directory):
  if not os.path.isdir(directory):
    os.makedirs(directory)

def getBase(filename):
  return os.path.splitext(filename)[0]

def getExt(filename):
  return os.path.splitext(filename)[1]

def getFileFolder(filepath):
  return os.path.dirname(filepath)

def getFileName(filepath):
  return os.path.basename(filepath)
  
### MS helpers ###

def precMzFromPrecMass(pmass, z):
  return (float(pmass) + 1.00727646677 * (int(z) - 1)) / int(z)

def precMassFromPrecMz(pmz, z):
  return pmz * z - 1.00727646677 * (z - 1)
