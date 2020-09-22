import unittest
from .. import digest

class DigestTest(unittest.TestCase):
  def test_digest_getDigestedPeptides_noMiscleavages(self):
    self.assertSetEqual(set(digest.getDigestedPeptides("APEPTIDERBPEPTIDEKCPEPTIDEK")), set(["APEPTIDER", "BPEPTIDEK", "CPEPTIDEK"]))
  
  def test_digest_getDigestedPeptides_noMiscleavagesNoTrypticEnd(self):
    self.assertSetEqual(set(digest.getDigestedPeptides("APEPTIDERBPEPTIDEKCPEPTIDE")), set(["APEPTIDER", "BPEPTIDEK", "CPEPTIDE"]))
  
  def test_digest_getDigestedPeptides_noMiscleavagesWithProlineBeforeTryptic(self):
    self.assertSetEqual(set(digest.getDigestedPeptides("APEPTIDERPBPEPTIDEKCPEPTIDEK")), set(["APEPTIDERPBPEPTIDEK", "CPEPTIDEK"]))
    
  def test_digest_getDigestedPeptides_oneMiscleavage(self):
    self.assertSetEqual(set(digest.getDigestedPeptides("APEPTIDERBPEPTIDEKCPEPTIDEK", miscleavages = 1)), set(["APEPTIDER", "BPEPTIDEK", "CPEPTIDEK", "APEPTIDERBPEPTIDEK", "BPEPTIDEKCPEPTIDEK"]))
    
  def test_digest_getDigestedPeptides_twoMiscleavages(self):
    self.assertSetEqual(set(digest.getDigestedPeptides("APEPTIDERBPEPTIDEKCPEPTIDEK", miscleavages = 2)), set(["APEPTIDER", "BPEPTIDEK", "CPEPTIDEK", "APEPTIDERBPEPTIDEK", "BPEPTIDEKCPEPTIDEK", "APEPTIDERBPEPTIDEKCPEPTIDEK"]))
  
  def test_digest_getDigestedPeptides_tooShort(self):
    self.assertSetEqual(set(digest.getDigestedPeptides("APEPRBPEPTIDEKCPEPTIDEK")), set(["BPEPTIDEK", "CPEPTIDEK"]))
  
  def test_digest_getDigestedPeptides_tooLong(self):
    self.assertSetEqual(set(digest.getDigestedPeptides("APEPTIDERBPEPTIDEKCPEPTIDEK", max_len = 7)), set())
  
  def test_digest_getDigestedPeptides_partialDigestion(self):
    self.assertSetEqual(set(digest.getDigestedPeptides("APEPTIDERB", digestion = 'semi')), set(["APEPTIDER", "APEPTIDE", "APEPTID", "APEPTI", "PEPTIDER", "EPTIDER", "PTIDER"]))
  
  def test_digest_getDigestedPeptides_partialDigestionTrypticEnd(self):
    self.assertSetEqual(set(digest.getDigestedPeptides("APEPTIDER", digestion = 'semi')), set(["APEPTIDER", "APEPTIDE", "APEPTID", "APEPTI", "PEPTIDER", "EPTIDER", "PTIDER"]))
    
  def test_digest_getDigestedPeptides_unspecificDigestion(self):
    self.assertSetEqual(set(digest.getDigestedPeptides("APEPTIDERB", digestion = 'none')), set(["APEPTIDER", "APEPTIDE", "APEPTID", "APEPTI", "PEPTIDER", "EPTIDER", "PTIDER", "PEPTIDE", "PEPTID", "EPTIDE", "APEPTIDERB", "PEPTIDERB", "EPTIDERB", "PTIDERB", "TIDERB"]))
  
if __name__ == '__main__':
  unittest.main()
