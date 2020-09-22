import unittest
from .. import helpers

class HelpersIOTest(unittest.TestCase):
  def test_helpers_getBase_withoutPath(self):
    self.assertEqual(helpers.getBase("myfile.txt"), "myfile")
  
  def test_helpers_getBase_withPath(self):
    self.assertEqual(helpers.getBase("/tmp/tmp2/myfile.txt"), "/tmp/tmp2/myfile")
  
  def test_helpers_getExt_withoutPath(self):
    self.assertEqual(helpers.getExt("myfile.txt"), ".txt")
  
  def test_helpers_getExt_withPath(self):
    self.assertEqual(helpers.getExt("/tmp/tmp2/myfile.txt"), ".txt")
    
  def test_helpers_getFileFolder_simple(self):
    self.assertEqual(helpers.getFileFolder("/tmp/tmp2/myfile.txt"), "/tmp/tmp2")
  
  def test_helpers_getFileName_simple(self):
    self.assertEqual(helpers.getFileName("/tmp/tmp2/myfile.txt"), "myfile.txt")


class HelpersMSTest(unittest.TestCase):
  def test_helpers_precMzFromPrecMass_simple(self):
    self.assertAlmostEqual(helpers.precMzFromPrecMass(1000.0, 2), 500.5036, places = 4)
  
  def test_helpers_precMassFromPrecMz_simple(self):
    self.assertAlmostEqual(helpers.precMassFromPrecMz(1000.0, 2), 1998.9927, places = 4)
  
  def test_helpers_precMassMzConversion_simple(self):
    self.assertAlmostEqual(helpers.precMassFromPrecMz(helpers.precMzFromPrecMass(1000.0, 2), 2), 1000.0, places = 4)


if __name__ == '__main__':
  unittest.main()
