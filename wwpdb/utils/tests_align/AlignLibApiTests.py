"""
File:    AlignLibApiTests.py

Test cases covering the PairwiseAlign and PseudoMultiAlign wrapper API using
small sequences with known alignments.

"""

import logging
import unittest
from typing import List, Tuple

from wwpdb.utils.align.alignlib import PairwiseAlign, PseudoMultiAlign  # pylint: disable=no-name-in-module,import-error
from wwpdb.utils.align.SequenceExamples import SequenceExamples

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)

REF_SEQ = ["ALA", "GLY", "SER", "THR", "VAL", "LEU", "ILE", "MET", "PHE", "TRP"]
GAP = "."


class PairwiseAlignApiTests(unittest.TestCase):
    def setUp(self) -> None:
        self.pA = PairwiseAlign()
        self.pA.setVerbose(False)
        self.pA.setReferenceSequence(REF_SEQ, "REF")

    def testIdentical(self) -> None:
        self.pA.addTestSequence(list(REF_SEQ), "T")
        self.pA.doAlign()
        self.assertEqual(self.pA.getAlignment("T"), list(zip(REF_SEQ, REF_SEQ)))

    def testDeletion(self) -> None:
        """Residues missing from the test sequence align to gaps"""
        self.pA.addTestSequence(REF_SEQ[:4] + REF_SEQ[6:], "T")
        self.pA.doAlign()
        myAlign = self.pA.getAlignment("T")
        self.assertEqual(len(myAlign), len(REF_SEQ))
        self.assertEqual([r for r, _ in myAlign], REF_SEQ)
        testCol = [t for _, t in myAlign]
        self.assertEqual(testCol, REF_SEQ[:4] + [GAP, GAP] + REF_SEQ[6:])
        self.assertEqual(self.pA.countGaps(testCol), 2)

    def testMismatch(self) -> None:
        tS = list(REF_SEQ)
        tS[5] = "TYR"
        self.pA.addTestSequence(tS, "T")
        self.pA.doAlign()
        conflicts = [(i, p) for i, p in enumerate(self.pA.getAlignment("T")) if p[0] != p[1]]
        self.assertEqual(conflicts, [(5, ("LEU", "TYR"))])

    def testCountGaps(self) -> None:
        self.assertEqual(self.pA.countGaps([]), 0)
        self.assertEqual(self.pA.countGaps(REF_SEQ), 0)
        self.assertEqual(self.pA.countGaps([GAP, "ALA", GAP, GAP]), 3)

    def testUnknownSequence(self) -> None:
        self.pA.addTestSequence(list(REF_SEQ), "T")
        self.pA.doAlign()
        self.assertEqual(self.pA.getAlignment("NOT_THERE"), [])

    def testClear(self) -> None:
        self.pA.addTestSequence(list(REF_SEQ), "T")
        self.pA.doAlign()
        self.assertEqual(len(self.pA.getAlignment("T")), len(REF_SEQ))
        self.pA.clear()
        self.assertEqual(self.pA.getAlignment("T"), [])
        # Usable again after clear
        self.pA.setReferenceSequence(REF_SEQ, "REF")
        self.pA.addTestSequence(REF_SEQ[:5], "T2")
        self.pA.doAlign()
        self.assertEqual([t for _, t in self.pA.getAlignment("T2")], REF_SEQ[:5] + [GAP] * 5)

    def testAddWithLink(self) -> None:
        """Return value is the number of test sequences stored"""
        self.assertEqual(self.pA.addTestSequenceWithLink(list(REF_SEQ), "L1", [1] * len(REF_SEQ)), 1)
        # Linkage length mismatch is ignored but the sequence is still added
        self.assertEqual(self.pA.addTestSequenceWithLink(list(REF_SEQ), "L2", [1, 1, 1]), 2)
        self.pA.doAlign()
        self.assertEqual(self.pA.getAlignment("L1"), list(zip(REF_SEQ, REF_SEQ)))
        self.assertEqual(self.pA.getAlignment("L2"), list(zip(REF_SEQ, REF_SEQ)))

    def testAddWithLinkAndRange(self) -> None:
        """Test sequence restricted to reference positions 2-7 (one based, inclusive begin)"""
        tS = REF_SEQ[2:8]
        nSeq = self.pA.addTestSequenceWithLinkAndRange(tS, "R", [1] * len(tS), 2, 7)
        self.assertEqual(nSeq, 1)
        self.pA.doAlign()
        myAlign = self.pA.getAlignment("R")
        self.assertEqual([r for r, _ in myAlign], REF_SEQ)
        self.assertEqual([t for _, t in myAlign], [GAP, GAP] + tS + [GAP, GAP])

    def testAlignConsensus(self) -> None:
        self.pA.addTestSequence(list(REF_SEQ), "T1")
        self.pA.addTestSequence(REF_SEQ[:4] + REF_SEQ[6:], "T2")
        myFails = self.pA.doAlignConsensus()
        self.assertEqual(myFails, [])
        self.assertEqual(self.pA.getAlignment("T1"), list(zip(REF_SEQ, REF_SEQ)))
        self.assertEqual(self.pA.countGaps([t for _, t in self.pA.getAlignment("T2")]), 2)

    def testMultipleAlign(self) -> None:
        self.pA.addTestSequence(list(REF_SEQ), "T1")
        self.pA.addTestSequence(REF_SEQ[:4] + REF_SEQ[6:], "T2")
        myFails = self.pA.doMultipleAlign()
        self.assertIsInstance(myFails, list)
        self.assertEqual(len(self.pA.getAlignment("T1")), len(REF_SEQ))

    def testPrintMethods(self) -> None:
        """Report methods write via C++ streams and must not raise"""
        tS = list(REF_SEQ)
        tS[5] = "TYR"
        self.pA.addTestSequence(tS, "T")
        self.pA.doAlignConsensus()
        self.pA.prAlignment("T")
        self.pA.prAlignmentConflicts("T")
        self.pA.prAlignmentFull()

    def testSelfAlignExample(self) -> None:
        """Example reference sequence aligned against itself has no conflicts"""
        refL = SequenceExamples().getRefSequence3List("A")
        pA = PairwiseAlign()
        pA.setVerbose(False)
        pA.setReferenceSequence(refL, "REFA")
        pA.addTestSequence(list(refL), "A")
        pA.doAlign()
        myAlign = pA.getAlignment("A")
        self.assertEqual(len(myAlign), len(refL))
        self.assertTrue(all(r == t for r, t in myAlign))

    def testAuthIsSubsequence(self) -> None:
        """Every author residue appears, in order, in the test column of the alignment"""
        sE = SequenceExamples()
        authL = sE.getAuthSequenceList("A")
        pA = PairwiseAlign()
        pA.setVerbose(False)
        pA.setReferenceSequence(sE.getRefSequence3List("A"), "REFA")
        pA.addTestSequence(authL, "A")
        pA.doAlign()
        testCol = [t for _, t in pA.getAlignment("A") if t != GAP]
        self.assertEqual(testCol, authL)


class PseudoMultiAlignApiTests(unittest.TestCase):
    def setUp(self) -> None:
        self.pM = PseudoMultiAlign()
        self.authSeq: List[Tuple[str, str]] = [(r, str(i)) for i, r in enumerate(REF_SEQ)]

    @staticmethod
    def __coorSeq(seq: List[str]) -> List[Tuple[str, str, str, str]]:
        return [(r, str(i), "1", str(i + 1)) for i, r in enumerate(seq)]

    @staticmethod
    def __coorSeqList(seq: List[str]) -> List[List[str]]:
        return [[r, str(i), "1", str(i + 1)] for i, r in enumerate(seq)]

    def testEmpty(self) -> None:
        self.assertEqual(self.pM.getAlignIndices(), [])
        self.assertEqual(self.pM.getAlignSequences(), [])

    def testIdentical(self) -> None:
        self.pM.setAuthSequence(self.authSeq)
        self.pM.addAlignSequence(self.__coorSeq(REF_SEQ))
        self.assertEqual(self.pM.getAlignIndices(), [[i, i] for i in range(len(REF_SEQ))])
        self.assertEqual(self.pM.getAlignSequences(), [[r, r] for r in REF_SEQ])

    def testDeletion(self) -> None:
        """Coordinate sequence missing two residues - all residues are placed in order"""
        coorL = REF_SEQ[:3] + REF_SEQ[5:]
        self.pM.setAuthSequence(self.authSeq)
        self.pM.addAlignSequence(self.__coorSeq(coorL))
        indices = self.pM.getAlignIndices()
        seqs = self.pM.getAlignSequences()
        self.assertEqual(len(indices), len(REF_SEQ))
        self.assertEqual(len(seqs), len(REF_SEQ))
        self.assertEqual([a for a, _ in indices], list(range(len(REF_SEQ))))
        self.assertEqual([c for _, c in indices if c >= 0], list(range(len(coorL))))
        self.assertEqual([c for _, c in seqs if c != GAP], coorL)
        for (_, c), (_, cIdx) in zip(seqs, indices):
            self.assertEqual(c == GAP, cIdx == -1)

    def testClear(self) -> None:
        self.pM.setAuthSequence(self.authSeq)
        self.pM.addAlignSequence(self.__coorSeq(REF_SEQ))
        self.assertEqual(len(self.pM.getAlignIndices()), len(REF_SEQ))
        self.pM.clear()
        self.assertEqual(self.pM.getAlignIndices(), [])
        self.assertEqual(self.pM.getAlignSequences(), [])

    def testAddWithRange(self) -> None:
        coorL = REF_SEQ[2:6]
        self.pM.setAuthSequence(self.authSeq)
        self.pM.addAlignSequenceWithRange(self.__coorSeqList(coorL), 2, 5)
        expIdx = [[i, i - 2 if 2 <= i <= 5 else -1] for i in range(len(REF_SEQ))]
        self.assertEqual(self.pM.getAlignIndices(), expIdx)
        expSeq = [[r, r if 2 <= i <= 5 else GAP] for i, r in enumerate(REF_SEQ)]
        self.assertEqual(self.pM.getAlignSequences(), expSeq)

    def testAddWithLinkageAndRange(self) -> None:
        self.pM.setAuthSequence(self.authSeq)
        self.pM.addAlignSequenceWithLinkageAndRange(self.__coorSeqList(REF_SEQ), [1] * len(REF_SEQ), 0, len(REF_SEQ) - 1)
        self.assertEqual(self.pM.getAlignIndices(), [[i, i] for i in range(len(REF_SEQ))])
        self.assertEqual(self.pM.getAlignSequences(), [[r, r] for r in REF_SEQ])

    def testScoringOptions(self) -> None:
        """Alternate scoring settings still align identical sequences exactly"""
        for setScore in (self.pM.setRefScore, self.pM.setAuthScore):
            self.pM.clear()
            self.pM.setPenaltyFactor(2.0)
            setScore()
            self.pM.setAuthSequence(self.authSeq)
            self.pM.addAlignSequence(self.__coorSeq(REF_SEQ))
            self.assertEqual(self.pM.getAlignIndices(), [[i, i] for i in range(len(REF_SEQ))])

    def testExampleXyzSequence(self) -> None:
        """Example chain B coordinate sequence placed against its author sequence"""
        sE = SequenceExamples()
        authL = sE.getAuthSequenceWithIndexList("B")
        xyzL = sE.getXyzSequenceWithIndexList("B")
        coorL = [(r, str(i), "1", num) for i, (r, num) in enumerate(xyzL)]
        self.pM.setAuthSequence(authL)
        self.pM.addAlignSequence(coorL)
        seqs = self.pM.getAlignSequences()
        self.assertEqual([c for _, c in seqs if c != GAP], [r for r, _ in xyzL])
        self.assertEqual([a for a, _ in seqs if a != GAP], [r for r, _ in authL])


if __name__ == "__main__":
    unittest.main()
