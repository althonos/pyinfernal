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

    def seqs_file(self, name, digital=False):
        path = resource_files("pyinfernal.tests").joinpath("data", "seqs", "{}.fa".format(name))
        if not path.exists():
            self.skipTest(f"data files not available: {str(path)!r}")
        return SequenceFile(path, digital=digital)

    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_trna(self):
        # $ cmsearch -Z1 tRNA.c.cm 100k-4.fa
        #
        # Query:       tRNA  [CLEN=71]
        # Accession:   RF00005
        # Description: tRNA
        # Hit scores:
        # rank     E-value  score  bias  sequence                       start    end   mdl trunc   gc  description
        # ----   --------- ------ -----  ----------------------------- ------ ------   --- ----- ----  -----------
        # (1) !   4.9e-15   65.4   0.0  tRNA-sample6/82250-82321       82250  82321 +  cm    no 0.51  -
        # (2) !   9.8e-14   60.8   0.0  tRNA-sample2/53622-53692       53622  53692 +  cm    no 0.58  -
        # (3) !   6.1e-12   54.4   0.0  tRNA-sample7/51572-51643       51572  51643 +  cm    no 0.38  -
        # (4) !   1.3e-11   53.2   0.0  tRNA-sample4/62081-62154       62081  62154 +  cm    no 0.53  -
        # (5) !   1.4e-10   49.6   0.0  tRNA-sample3/94923-94993       94923  94993 +  cm    no 0.56  -
        # (6) !   5.4e-09   43.8   0.0  tRNA-sample1/46421-46490       46421  46490 +  cm    no 0.40  -
        # (7) !   1.9e-08   41.9   0.0  tRNA-sample5/56066-56134       56066  56134 +  cm    no 0.54  -
        # (8) !   3.3e-08   41.1   0.0  tRNA-sample8/12834-12904       12834  12904 +  cm    no 0.38  -
        # (9) !   1.3e-07   38.9   0.0  tRNA-sample10/32353-32425      32353  32425 +  cm    no 0.40  -
        # (10) !   2.3e-05   30.9   0.0  tRNA-sample9/66515-66583       66516  66582 +  cm    no 0.57  -
        # ------ inclusion threshold ------
        # (11) ?     0.014   21.0   0.1  tRNA-sample10/32353-32425      59443  59368 -  cm    no 0.30  -
        # (12) ?     0.026   20.0   0.1  Plant_SRP-sample2/44686-44961  42608  42527 -  cm    no 0.34  -
        # (13) ?       0.4   15.8   0.0  snR75-sample7/20908-21006      86357  86283 -  cm    no 0.35  -
        # (14) ?      0.67   15.0   0.0  Vault-sample2/67800-67897      89308  89250 -  cm    no 0.34  -
        # (15) ?      0.83   14.6   3.9  tRNA-sample8/12834-12904       46204  46146 -  cm    no 0.19  -
        # (16) ?       2.8   12.7   3.1  Plant_SRP-sample4/64689-64987  97544  97488 -  cm    no 0.21  -
        # (17) ?       4.2   12.1   0.0  Plant_SRP-sample5/23875-24163  93227  93282 +  cm    no 0.27  -
        # (18) ?       6.5   11.4   0.0  Plant_SRP-sample2/44686-44961  95590  95544 -  cm    no 0.38  -
        # (19) ?       9.3   10.9   3.3  snR75-sample7/20908-21006      14188  14138 -  cm    no 0.18  -
        #
        # Internal CM pipeline statistics summary:
        # ----------------------------------------
        # Query model(s):                                                  1  (71 consensus positions)
        # Target sequences:                                               40  (8000000 residues searched)
        # Target sequences re-searched for truncated hits:                40  (34880 residues re-searched)
        # Windows   passing  local HMM SSV           filter:           10984  (0.342); expected (0.35)
        # Windows   passing  local HMM Viterbi       filter:                  (off)
        # Windows   passing  local HMM Viterbi  bias filter:                  (off)
        # Windows   passing  local HMM Forward       filter:             828  (0.02853); expected (0.02)
        # Windows   passing  local HMM Forward  bias filter:             694  (0.02375); expected (0.02)
        # Windows   passing glocal HMM Forward       filter:             381  (0.01398); expected (0.02)
        # Windows   passing glocal HMM Forward  bias filter:             371  (0.01362); expected (0.02)
        # Envelopes passing glocal HMM envelope defn filter:             363  (0.004153); expected (0.02)
        # Envelopes passing  local CM  CYK           filter:              26  (0.0002083); expected (0.0001)
        # Total CM hits reported:                                         19  (0.0001583); includes 0 truncated hit(s)

        with self.cm_file("tRNA.c") as cm_file:
            cm = cm_file.read()
        with self.seqs_file("100k-4", digital=True) as seqs_file:
            seqs = seqs_file.read_block()

        hits = self.get_hits(cm, seqs, Z=1e6)  # 1 Mbp
        self.assertEqual(len(hits.reported), 19)
        self.assertEqual(len(hits.included), 10)

        with self.table("tRNA.Z1.tbl") as tbl:
            lines = filter(lambda line: not line.startswith("#"), tbl)
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

    # @unittest.skipUnless(resource_files, "importlib.resources not available")
    # def test_5c(self):
    #     with self.cm_file("5.c") as cm_file:
    #         cms = list(cm_file)
    #     with self.seqs_file("100k-4", digital=True) as seqs_file:
    #         seqs = seqs_file.read_block()

    #     all_hits = self.get_hits_multi(cms, seqs, Z=1e6)  # 1 Mbp
    #     all_hits.sort(key=lambda hits: ["tRNA", "Vault", "snR75", "Plant_SRP", "tRNA-Sec"].index(hits.query.name))

    #     nreported = sum(len(hits.reported) for hits in all_hits)
    #     nincluded = sum(len(hits.included) for hits in all_hits)
    #     self.assertEqual(nreported, 71)
    #     self.assertEqual(nincluded, 39)

    #     with self.table("5c.Z1.tbl") as tbl:
    #         lines = filter(lambda line: not line.startswith("#"), tbl)
    #         hits_it = itertools.chain.from_iterable(all_hits)
    #         for hit, line in itertools.zip_longest(hits_it, lines):
    #             self.assertIsNot(line, None)
    #             self.assertIsNot(hit, None)

    #             fields = line.split()

    #             self.assertEqual(hit.name, fields[0])
    #             self.assertEqual(hit.accession, None if fields[1] == "-" else fields[1])
    #             self.assertEqual(hit.hits.query.name, fields[2])
    #             self.assertEqual(hit.hits.query.accession, fields[3])
    #             self.assertEqual(hit.alignment.cm_from, int(fields[5]))
    #             self.assertEqual(hit.alignment.cm_to, int(fields[6]))
    #             self.assertEqual(hit.alignment.target_from, int(fields[7]))
    #             self.assertEqual(hit.alignment.target_to, int(fields[8]))
    #             self.assertEqual(hit.strand, fields[9])
    #             # self.assertEqual(hit.trunc, fields[10])
    #             # self.assertEqual(hit.pipeline_pass, fields[11])
    #             # self.assertEqual(hit.gc, fields[12])
    #             self.assertAlmostEqual(hit.bias, float(fields[13]), places=1)
    #             self.assertAlmostEqual(hit.score, float(fields[14]), places=1)
    #             self.assertAlmostEqual(hit.evalue, float(fields[15]), delta=hit.evalue / 10)
                
    #             if fields[16] == "!":
    #                 self.assertTrue(hit.included)
    #             elif fields[16] == "?":
    #                 self.assertTrue(hit.reported)
    
    @unittest.skipUnless(resource_files, "importlib.resources not available")
    def test_snR75_T50(self):
        with self.cm_file("5.c") as cm_file:
            cm = next(cm for cm in cm_file if cm.name == "snR75")
        with self.seqs_file("100k-4", digital=True) as seqs_file:
            seqs = seqs_file.read_block()

        hits = self.get_hits(cm, seqs, T=50, Z=1e6)  # 1 Mbp

        self.assertEqual(len(hits.reported), 1)
        self.assertEqual(len(hits.included), 1)

        with self.table("snR75.Z1.T50.tbl") as tbl:
            lines = filter(lambda line: not line.startswith("#"), tbl)
            for hit, line in itertools.zip_longest(hits, lines):
                self.assertIsNot(line, None)
                self.assertIsNot(hit, None)

                fields = line.split()

                self.assertEqual(hit.name, fields[0])
                self.assertEqual(hit.accession, None if fields[1] == "-" else fields[1])
                # self.assertEqual(hits.query.name, fields[2])
                # self.assertEqual(hits.query.accession, fields[3])
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
    def test_5c_T50(self):
        with self.cm_file("5.c") as cm_file:
            cms = list(cm_file)
        with self.seqs_file("100k-4", digital=True) as seqs_file:
            seqs = seqs_file.read_block()

        all_hits = self.get_hits_multi(cms, seqs, T=50, Z=1e6)  # 1 Mbp
        all_hits.sort(key=lambda hits: ["tRNA", "Vault", "snR75", "Plant_SRP", "tRNA-Sec"].index(hits.query.name))
        
        nreported = sum(len(hits.reported) for hits in all_hits)
        nincluded = sum(len(hits.included) for hits in all_hits)
        self.assertEqual(nreported, 17)
        self.assertEqual(nincluded, 17)

        with self.table("5c.Z1.T50.tbl") as tbl:
            lines = filter(lambda line: not line.startswith("#"), tbl)
            hits_it = itertools.chain.from_iterable(all_hits)
            for hit, line in itertools.zip_longest(hits_it, lines):
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
