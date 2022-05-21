import os
import sys
import unittest

from ralf import get_rotation_limit


class Test_ralf(unittest.TestCase):
    """
    Test the various functionalities of ralf.
    """

    @classmethod
    def setUpClass(self):
        self._knipholone_pdb = os.path.join(
            os.getcwd(),
            "test",
            "data",
            "knipholone.pdb",
        )
        self._binap_pdb = os.path.join(
            os.getcwd(),
            "test",
            "data",
            "binap.pdb",
        )
        self._segphos_pdb = os.path.join(
            os.getcwd(),
            "test",
            "data",
            "segphos.pdb",
        )
        self._tolbinap_pdb = os.path.join(
            os.getcwd(),
            "test",
            "data",
            "tolbinap.pdb",
        )

    def test_knipholone(self):
        """knipholone natrual product"""
        get_rotation_limit(
            self._knipholone_pdb,
            14,
            7,
        )

    def test_binap(self):
        """binap ligand"""
        get_rotation_limit(
            self._binap_pdb,
            23,
            24,
        )

    def test_segphos(self):
        """segphos ligand"""
        get_rotation_limit(
            self._segphos_pdb,
            23,
            6,
        )

    def test_tolbinap(self):
        """tolbinap ligand"""
        get_rotation_limit(
            self._tolbinap_pdb,
            25,
            26,
        )

    def test_custom_cutoff(self):
        """run with custom cutoff distance"""
        get_rotation_limit(
            self._segphos_pdb,
            23,
            6,
            cutoff_distance=1.4,
        )

    def test_unrestricted_mol(self):
        """non-atropisomer should raise an error"""
        with self.assertRaises(RuntimeError):
            raise RuntimeError

    def test_nonexistent_atom_ids(self):
        pass

    def test_unbonded_atom_ids(self):
        pass


if __name__ == "__main__":
    unittest.main()
