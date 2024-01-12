"""
test how FF_function will be called and how the exception and error msg be handled/present
"""
import FF_functions.frag_gen_all as frag_gen_all
import unittest
import subprocess
from utilities.parse import read_whole_file

class test_FFfunctionCall(unittest.TestCase):
    def test_regular_shell_call(self):
        # for shell call the exception will be handeled and the error msg will be saved to log
        subprocess.call("python ../FF_functions/frag_gen_all.py \"1.xyz 2.xyz\" -q 1", shell=True)
        lines = read_whole_file('frag_gen_all.log')
        self.assertEqual(lines[-1], 'ERROR: Could not find file 1.xyz 2.xyz. Exiting...')
    def test_function_call(self):
        # exception will be raised if this is a function call
        frag_gen_all.run(xyz="1.xyz 2.xyz", q_tot=1)



