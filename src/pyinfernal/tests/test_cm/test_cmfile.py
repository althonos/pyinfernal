import io
import itertools
import os
import shutil
import unittest
import tempfile

from pyinfernal.cm import CMFile

from .. import __name__ as __package__
from .utils import INFERNAL_FOLDER, resource_files

# --- Mixins -------------------------------------------------------------------

class _TestCMFile:

    ID = NotImplemented
    NAMES = NotImplemented

    @classmethod
    def setUpClass(cls):
        cls.cms_folder = resource_files(__package__).joinpath("data", "cms")

    def open_cm(self, path):
        raise NotImplementedError()

    def check_cmfile(self, cmfile):
        for cm, expected in itertools.zip_longest(cmfile, self.NAMES):
            self.assertIsNot(cm, None, "missing CM: {}".format(expected))
            self.assertIsNot(expected, None, "unexpected extra CM: {}".format(cm))
            self.assertEqual(cm.name, expected)
            # self.assertIsNot(hmm.cutoffs, None)
            # self.assertIsNot(hmm.evalue_parameters, None)

    def test_empty(self):
        with tempfile.NamedTemporaryFile() as empty:
            self.assertRaises(EOFError, self.open_cm, empty.name)

    def test_read(self):
        path = self.cms_folder.joinpath("{}.cm".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_cm(path) as f:
            self.check_cmfile(f)


class _TestCMFileFileobj:

    def open_cm(self, path):
        with open(path, "rb") as f:
            buffer = io.BytesIO(f.read())
        return CMFile(buffer)

    def test_name(self):
        path = self.cms_folder.joinpath("{}.cm".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_cm(path) as f:
            self.assertIs(f.name, None)


class _TestCMFilePath:

    def open_cm(self, path):
        return CMFile(path)

    def test_name(self):
        path = self.cms_folder.joinpath("{}.cm".format(self.ID))
        if not path.exists():
            self.skipTest("data files not available")
        with self.open_cm(path) as f:
            self.assertEqual(f.name, str(path))
    


class _TestRF5c(_TestCMFile):
    ID = "5.c"
    NAMES = [
        b"tRNA",
        b"Vault",
        b"snR75",
        b"Plant_SRP",
        b"tRNA-Sec",
    ]


# --- Test cases ---------------------------------------------------------------

@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestRF5cFileobj(_TestCMFileFileobj, _TestRF5c, unittest.TestCase):
    pass

@unittest.skipUnless(resource_files, "importlib.resources.files not available")
class TestRF5cPath(_TestCMFilePath, _TestRF5c, unittest.TestCase):
    pass
