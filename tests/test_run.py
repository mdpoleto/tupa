import unittest
import subprocess
import os

from grep_run_values import grep_efield_value

###################

path = "tests"
runline = "TUPA.py -top " + path + "/unit_test.psf -traj " + path + "/unit_test.pdb -config " + path + "/config.conf -outdir " + path + "/tupa_test"

x = subprocess.call(runline, shell=True, stdout=open(os.devnull, 'wb'))

outpath = os.path.abspath(path + "/tupa_test")

magnitude = 14.399652
field_x   = 10.182092
field_y   = 10.182092
field_z   = 0.0

class CheckValues(unittest.TestCase):

    def test_elecfield_magnitude(self):
        v = grep_efield_value(1, outpath)
        self.assertEqual(v, magnitude)

    def test_elecfield_field_x(self):
        v = grep_efield_value(2, outpath)
        self.assertEqual(v, field_x)

    def test_elecfield_field_y(self):
        v = grep_efield_value(3, outpath)
        self.assertEqual(v, field_y)

    def test_elecfield_field_z(self):
        v = grep_efield_value(4, outpath)
        self.assertEqual(v, field_z)

unittest.main()




