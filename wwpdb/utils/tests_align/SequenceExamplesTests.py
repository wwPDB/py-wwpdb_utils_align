"""
File:    SequenceExamplesTests.py

Test cases for the SequenceExamples example data and helper methods.

"""

import io
import logging
import runpy
import unittest
import warnings
from typing import List
from unittest import mock

from wwpdb.utils.align.SequenceExamples import SequenceExamples, formatCpp

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)

_RANDOM_TARGET = "wwpdb.utils.align.SequenceExamples.random"


class SequenceExamplesTests(unittest.TestCase):
    def setUp(self) -> None:
        self.sE = SequenceExamples()

    def testIdCode(self) -> None:
        self.assertEqual(self.sE.getIdCode(), "3IJE")
        self.assertEqual(SequenceExamples("1ABC").getIdCode(), "1ABC")

    def testSeqIdList(self) -> None:
        self.assertEqual(self.sE.getSeqIdList(), ["A", "B"])

    def testUnknownChain(self) -> None:
        """All accessors return an empty list for an unknown chain"""
        self.assertEqual(self.sE.getAuthSequenceList("Z"), [])
        self.assertEqual(self.sE.getAuthSequenceWithIndexList("Z"), [])
        self.assertEqual(self.sE.getXyzSequenceWithIndexList("Z"), [])
        self.assertEqual(self.sE.getRefSequenceList("Z"), [])
        self.assertEqual(self.sE.getRefSequenceWithIndexList("Z"), [])
        self.assertEqual(self.sE.getRefSequence3List("Z"), [])

    def testAuthSequenceList(self) -> None:
        for chainId in self.sE.getSeqIdList():
            sL = self.sE.getAuthSequenceList(chainId)
            self.assertGreater(len(sL), 0)
            for r in sL:
                self.assertIn(r, self.sE.aaList)
        self.assertEqual(self.sE.getAuthSequenceList("A")[:4], ["PHE", "ASN", "LEU", "ASP"])

    def testAuthSequenceWithIndexList(self) -> None:
        for chainId in self.sE.getSeqIdList():
            sL = self.sE.getAuthSequenceList(chainId)
            tL = self.sE.getAuthSequenceWithIndexList(chainId)
            self.assertEqual(len(tL), len(sL))
            for idx, (r, sIdx) in enumerate(tL):
                self.assertEqual(r, sL[idx])
                self.assertEqual(sIdx, str(idx))

    def testXyzSequenceWithIndexList(self) -> None:
        tLA = self.sE.getXyzSequenceWithIndexList("A")
        self.assertEqual(tLA[-1], ("ILE", "967"))
        tLB = self.sE.getXyzSequenceWithIndexList("B")
        self.assertEqual(tLB[0], ("GLY", "1"))
        self.assertEqual(tLB[-1], ("VAL", "695"))

    def testToList(self) -> None:
        self.assertEqual(self.sE.toList(" AB\nC\tD "), ["A", "B", "C", "D"])
        self.assertEqual(self.sE.toList(""), [])

    def testRefSequenceList(self) -> None:
        for chainId in self.sE.getSeqIdList():
            s1L = self.sE.getRefSequenceList(chainId)
            s3L = self.sE.getRefSequence3List(chainId)
            self.assertGreater(len(s1L), 0)
            self.assertEqual(len(s1L), len(s3L))
            for r in s1L:
                self.assertEqual(len(r), 1)
            for r in s3L:
                self.assertIn(r, self.sE.aaList)
        self.assertEqual(self.sE.getRefSequence3List("A")[:5], ["MET", "ALA", "PHE", "PRO", "PRO"])

    def testRefSequence3ListUnknownResidue(self) -> None:
        """One letter codes without a three letter mapping become UNK"""
        self.sE.refSequenceA = "MXA"
        self.assertEqual(self.sE.getRefSequence3List("A"), ["MET", "UNK", "ALA"])

    def testRefSequenceWithIndexList(self) -> None:
        for chainId in self.sE.getSeqIdList():
            s3L = self.sE.getRefSequence3List(chainId)
            tL = self.sE.getRefSequenceWithIndexList(chainId)
            self.assertEqual([r for r, _ in tL], s3L)
            self.assertEqual([sIdx for _, sIdx in tL], [str(i) for i in range(len(s3L))])

    def testMutatedSequenceNoChange(self) -> None:
        """Random counts of one produce no insertions or deletions"""
        with mock.patch(_RANDOM_TARGET) as mRandom:
            mRandom.randint.return_value = 1
            seqT = self.sE.getAuthSequenceListTest("A")
        self.assertEqual(seqT, self.sE.getAuthSequenceList("A"))

    def testMutatedSequenceInsert(self) -> None:
        with mock.patch(_RANDOM_TARGET) as mRandom:
            # one single residue insertion at 0, no run insertion, no deletion
            mRandom.randint.side_effect = [2, 0, 1, 1]
            mRandom.choice.return_value = "UNK"
            seqT = self.sE.getAuthSequenceListTest("A")
        self.assertEqual(seqT, ["UNK"] + self.sE.getAuthSequenceList("A"))

    def testMutatedSequenceInsertRun(self) -> None:
        with mock.patch(_RANDOM_TARGET) as mRandom:
            # no single insertion, one run of two GLY at 5, no deletion
            mRandom.randint.side_effect = [1, 2, 5, 3, 1]
            mRandom.choice.return_value = "GLY"
            seqT = self.sE.getAuthSequenceListTest("A")
        orig = self.sE.getAuthSequenceList("A")
        self.assertEqual(seqT, orig[:5] + ["GLY", "GLY"] + orig[5:])

    def testMutatedSequenceDelete(self) -> None:
        with mock.patch(_RANDOM_TARGET) as mRandom:
            # no insertions, one deletion at 0
            mRandom.randint.side_effect = [1, 1, 2, 0]
            seqT = self.sE.getAuthSequenceListTest("A")
        self.assertEqual(seqT, self.sE.getAuthSequenceList("A")[1:])

    def testMutatedSequenceDeleteOutOfRange(self) -> None:
        """A deletion past the end of the sequence is ignored"""
        orig = self.sE.getAuthSequenceList("A")
        with mock.patch(_RANDOM_TARGET) as mRandom:
            mRandom.randint.side_effect = [1, 1, 2, len(orig)]
            seqT = self.sE.getAuthSequenceListTest("A")
        self.assertEqual(seqT, orig)

    def testMutatedSequenceRandom(self) -> None:
        """Unpatched mutation leaves the source data unchanged"""
        orig = self.sE.getAuthSequenceList("A")
        for _ in range(5):
            seqT = self.sE.getAuthSequenceListTest("A")
            for r in seqT:
                self.assertIn(r, self.sE.aaList)
            # at most 4 single insertions, 9 run insertions, 4 deletions
            self.assertLessEqual(abs(len(seqT) - len(orig)), 13)
        self.assertEqual(self.sE.getAuthSequenceList("A"), orig)


class FormatCppTests(unittest.TestCase):
    def __format(self, sL: List[str], seqName: str) -> str:
        with mock.patch("sys.stdout", new_callable=io.StringIO) as mOut:
            formatCpp(sL, seqName)
            return str(mOut.getvalue())

    def testFormatCpp(self) -> None:
        sL = ["ALA"] * 12
        out = self.__format(sL, "test")
        lines = out.splitlines()
        self.assertEqual(lines[0], "// test Length = 12 ")
        self.assertEqual(lines[1], "static const char *_testSequence[] = {")
        self.assertEqual(lines[2], '"ALA",' * 10)
        self.assertEqual(lines[3], '"ALA",' * 2)
        self.assertEqual(lines[4], '""};')
        self.assertEqual(len(lines), 5)

    def testFormatCppEmpty(self) -> None:
        out = self.__format([], "empty")
        self.assertEqual(out.splitlines(), ["// empty Length = 0 ", "static const char *_emptySequence[] = {", '""};'])

    def testMain(self) -> None:
        """Module entry point writes C++ arrays for reference, author and 10 mutated sequences"""
        with mock.patch("sys.stdout", new_callable=io.StringIO) as mOut, warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            runpy.run_module("wwpdb.utils.align.SequenceExamples", run_name="__main__")
            out = str(mOut.getvalue())
        self.assertIn("_refASequence[]", out)
        self.assertIn("_authASequence[]", out)
        self.assertEqual(out.count("static const char"), 12)


if __name__ == "__main__":
    unittest.main()
