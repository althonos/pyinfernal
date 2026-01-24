import abc
import collections
import math
import io
import itertools
import os
import platform
import unittest
import tempfile
import threading
import multiprocessing.resource_sharer

import pyhmmer
import pyinfernal
from pyhmmer.easel import Alphabet, DigitalMSA, MSAFile, SequenceFile, TextSequence
from pyinfernal.cm import CM, CMFile, TopHits, Hit, Alignment, Pipeline

from ..utils import resource_files

class _TestSearch(metaclass=abc.ABCMeta):

    def tearDown(self):
        multiprocessing.resource_sharer.stop()
        self.assertEqual(threading.active_count(), 1, threading.enumerate())

    @abc.abstractmethod
    def get_hits(self, cm, sequences, **options):
        return NotImplemented

    def get_hits_multi(self, cms, sequences, **options):
        return [self.get_hits(cm, sequences, **options) for cm in cms]

    def table(self, name):
        path = resource_files("pyinfernal.tests").joinpath("data", "tables", name)
        if not path.exists():
            self.skipTest("data files not available")
        return path.open()

    def cm_file(self, name):
        path = resource_files("pyinfernal.tests").joinpath("data", "cms", "{}.cm".format(name))
        if not path.exists():
            self.skipTest(f"data files not available: {str(path)!r}")
        return CMFile(path)

    def seqs_file(self, name, digital=False, alphabet=None):
        path = resource_files("pyinfernal.tests").joinpath("data", "seqs", "{}.fa".format(name))
        if not path.exists():
            self.skipTest(f"data files not available: {str(path)!r}")
        return SequenceFile(path, digital=digital, alphabet=alphabet)

    def assert_hits_match_table(self, hits, table):
        lines = filter(lambda line: not line.startswith("#"), table)
        for hit, line in itertools.zip_longest(hits, lines):
            self.assertIsNot(line, None)
            self.assertIsNot(hit, None)

            fields = line.split()

            self.assertEqual(hit.name, fields[0])
            self.assertEqual(hit.accession, None if fields[1] == "-" else fields[1])
            self.assertEqual(hit.hits.query.name, fields[2])
            self.assertEqual(hit.hits.query.accession, fields[3])
            self.assertEqual(hit.alignment.cm_from, int(fields[5]))
            self.assertEqual(hit.alignment.cm_to, int(fields[6]))
            self.assertEqual(hit.alignment.target_from, int(fields[7]))
            self.assertEqual(hit.alignment.target_to, int(fields[8]))
            self.assertEqual(hit.strand, fields[9])
            # self.assertEqual(hit.trunc, fields[10])
            # self.assertEqual(hit.pipeline_pass, fields[11])
            # self.assertEqual(hit.gc, fields[12])
            self.assertAlmostEqual(hit.bias, float(fields[13]), places=1)
            self.assertAlmostEqual(hit.score, float(fields[14]), places=1)
            self.assertAlmostEqual(hit.evalue, float(fields[15]), delta=hit.evalue / 10)
            
            if fields[16] == "!":
                self.assertTrue(hit.included)
            elif fields[16] == "?":
                self.assertTrue(hit.reported)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_pANT_5c(self):
        with self.cm_file("5.c") as cm_file:
            cms = list(cm_file)
        with self.seqs_file("pANT_R100", digital=True, alphabet=cms[0].alphabet) as seqs_file:
            seqs = seqs_file.read_block()

        all_hits = self.get_hits_multi(cms, seqs, Z=1e5)  # 0.1 Mbp
        all_hits.sort(key=lambda hits: ["tRNA", "Vault", "snR75", "Plant_SRP", "tRNA-Sec"].index(hits.query.name))
        
        nreported = sum(len(hits.reported) for hits in all_hits)
        nincluded = sum(len(hits.included) for hits in all_hits)
        self.assertEqual(nreported, 6)
        self.assertEqual(nincluded, 0)

        with self.table("pANT.5c.Z0_1.tbl") as tbl:
            hits_it = itertools.chain.from_iterable(all_hits)
            self.assert_hits_match_table(hits_it, tbl)

    def _test_pANT_RF(self, rfam_id, nincluded, nreported):
        with self.cm_file(rfam_id) as cm_file:
            cms = list(cm_file)
        with self.seqs_file("pANT_R100", digital=True, alphabet=cms[0].alphabet) as seqs_file:
            seqs = seqs_file.read_block()

        all_hits = self.get_hits_multi(cms, seqs, Z=1e5)  # 0.1 Mbp
        
        reported = sum(len(hits.reported) for hits in all_hits)
        included = sum(len(hits.included) for hits in all_hits)
        self.assertEqual(reported, nreported)
        self.assertEqual(included, nincluded)

        with self.table(f"pANT.{rfam_id}.Z0_1.tbl") as tbl:
            hits_it = itertools.chain.from_iterable(all_hits)
            self.assert_hits_match_table(hits_it, tbl)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_pANT_RF00029(self):
        self._test_pANT_RF("RF00029", 3, 3)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_pANT_RF00042(self):
        self._test_pANT_RF("RF00042", 2, 4)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_pANT_RF00107(self):
        self._test_pANT_RF("RF00107", 2, 5)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_pANT_RF00243(self):
        self._test_pANT_RF("RF00243", 1, 5)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_pANT_RF03523(self):
        self._test_pANT_RF("RF03523", 2, 2)
    
    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_pANT_RF_all(self):
        cms = []
        rfam_ids = ["RF00029", "RF00042", "RF00107", "RF00243", "RF03523"]
        for rfam_id in rfam_ids:
            with self.cm_file(rfam_id) as cm_file:
                cms.extend(cm_file)
        with self.seqs_file("pANT_R100", digital=True, alphabet=cms[0].alphabet) as seqs_file:
            seqs = seqs_file.read_block()

        all_hits = self.get_hits_multi(cms, seqs, Z=1e5)  # 0.1 Mbp
        all_hits.sort(key=lambda hits: rfam_ids.index(hits.query.accession))
        
        reported = sum(len(hits.reported) for hits in all_hits)
        included = sum(len(hits.included) for hits in all_hits)
        self.assertEqual(reported, 19)
        self.assertEqual(included, 10)

        with self.table(f"pANT.RFall.Z0_1.tbl") as tbl:
            hits_it = itertools.chain.from_iterable(all_hits)
            self.assert_hits_match_table(hits_it, tbl)


class TestCmsearch(_TestSearch, unittest.TestCase):
    parallel = "queries"

    @staticmethod
    def _random_sequences(n=20):
        rng = pyhmmer.easel.Randomness(42)
        alphabet = Alphabet.amino()
        return pyhmmer.easel.DigitalSequenceBlock(
            alphabet,
            [
                pyhmmer.easel.DigitalSequence.sample(alphabet, 200, rng)
                for _ in range(10)
            ]
        )

    def get_hits(self, cm, seqs, **options):
        return list(pyinfernal.cmsearch(cm, seqs, parallel=self.parallel, **options))[0]

    def get_hits_multi(self, cms, seqs, **options):
        return list(pyinfernal.cmsearch(cms, seqs, parallel=self.parallel, **options))

    # def test_callback_error_single_threaded(self):

    #     class MyException(Exception):
    #         pass

    #     def callback(cm, total):
    #         raise MyException("oopsie")

    #     rng = pyhmmer.easel.Randomness(42)
    #     alphabet = Alphabet.amino()
    #     hmm = HMM.sample(alphabet, 100, rng)
    #     seqs = self._random_sequences()

    #     hits = pyhmmer.hmmsearch(cm, seqs, cpus=1, callback=callback, parallel=self.parallel)
    #     with self.assertRaises(MyException):
    #         hit = next(hits)

    # def test_callback_error_multi_threaded(self):

    #     class MyException(Exception):
    #         pass

    #     def callback(cm, total):
    #         raise MyException("oopsie")

    #     rng = pyhmmer.easel.Randomness(42)
    #     alphabet = Alphabet.amino()
    #     hmm = HMM.sample(alphabet, 100, rng)
    #     seqs = self._random_sequences()

    #     hits = pyhmmer.hmmsearch(cm, seqs, cpus=2, callback=callback, parallel=self.parallel)
    #     with self.assertRaises(MyException):
    #         hit = next(hits)

    # def test_background_error(self):
    #     # check that errors occuring in worker threads are recovered and raised
    #     # in the main threads (a common error is mismatching the HMM and the
    #     # sequence alphabets).
    #     rng = pyhmmer.easel.Randomness(42)
    #     seqs = [TextSequence().digitize(Alphabet.dna())]
    #     hmm = HMM.sample(Alphabet.amino(), 100, rng)
    #     self.assertRaises(ValueError, self.get_hits, cm, seqs)

    def test_no_queries(self):
        seqs = self._random_sequences()
        hits = pyinfernal.cmsearch([], seqs, parallel=self.parallel)
        self.assertIs(None, next(hits, None))


class TestCmsearchSingle(TestCmsearch, unittest.TestCase):

    def get_hits(self, cm, seqs, **options):
        return list(pyinfernal.cmsearch(cm, seqs, cpus=1, parallel=self.parallel, **options))[0]

    def get_hits_multi(self, cms, seqs, **options):
        return list(pyinfernal.cmsearch(cms, seqs, cpus=1, parallel=self.parallel, **options))

    def test_no_queries(self):
        seqs = self._random_sequences()
        hits = pyinfernal.cmsearch([], seqs, cpus=1, parallel=self.parallel)
        self.assertIs(None, next(hits, None))

# @unittest.skipIf(platform.system() == "Windows", "may deadlock on Windows")
# @unittest.skipIf(platform.system() == "Darwin", "may deadlock on MacOS")
# @unittest.skipIf(platform.system() == "Emscripten", "no process support on Emscripten")
# class TestCmsearchProcess(TestCmsearch, unittest.TestCase):
#     def get_hits(self, cm, seqs, **options):
#         return list(pyinfernal.cmsearch(cm, seqs, cpus=2, backend="multiprocessing", parallel=self.parallel, **options))[0]

#     def get_hits_multi(self, cms, seqs, **options):
#         return list(pyinfernal.cmsearch(cms, seqs, cpus=2, backend="multiprocessing", parallel=self.parallel, **options))

#     def test_no_queries(self):
#         seqs = self._random_sequences()
#         hits = pyinfernal.cmsearch([], seqs, cpus=2, backend="multiprocessing", parallel=self.parallel)
#         self.assertIs(None, next(hits, None))


class TestCmsearchReverse(TestCmsearch):
    parallel = "targets"


class TestCmsearchReverseSingle(TestCmsearchSingle):
    parallel = "targets"


# @unittest.skipIf(platform.system() == "Windows", "may deadlock on Windows")
# @unittest.skipIf(platform.system() == "Darwin", "may deadlock on MacOS")
# @unittest.skipIf(platform.system() == "Emscripten", "no process support on Emscripten")
# class TestCmsearchReverseProcess(TestCmsearchProcess):
#     parallel = "targets"


class TestPipelinesearch(_TestSearch, unittest.TestCase):

    def get_hits(self, cm, seqs, **options):
        pipeline = Pipeline(alphabet=cm.alphabet, **options)
        hits = pipeline.search_cm(cm, seqs)
        return hits
