"""
File:    PairwiseAlignTests.py
Author:  jdw
Date:    11-Jan-2009
Version: 0.001

Updates:

14-Oct-2018  jdw updated for Python packaging, Pybind11 wrapping, and Py3

Test cases for pairwise sequence alignment wrapper class.

"""

import logging
import unittest

from wwpdb.utils.align.alignlib import PairwiseAlign  # pylint: disable=no-name-in-module
from wwpdb.utils.align.SequenceExamples import SequenceExamples

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s')
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class PairwiseAlignTests(unittest.TestCase):

    def setUp(self):
        self.__verbose = True
        sE = SequenceExamples()
        self.seqRef3L = sE.getRefSequence3List('A')
        self.seqAuth3L = sE.getAuthSequenceList('A')
        #
        # Test sequence with random insertions and deletions
        #
        self.sTests = {}
        for tt in ['T1', 'T2', 'T3', 'T4', 'T5']:
            self.sTests[tt] = sE.getAuthSequenceListTest('A')
        #

    def tearDown(self):
        pass

    def testAlign1(self):
        """ Run internal alignment test embedded in the class -  This is a basic santity check.
        """
        logger.info(" -------------------------")
        try:
            pA = PairwiseAlign()
            pA.testExample()
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testAlign2(self):
        """ Run author vs reference sequence alignment returning a copy of the alignment
            via getAlignment() -
        """
        logger.info("------------------------")
        try:
            pA = PairwiseAlign()
            pA.setVerbose(self.__verbose)
            pA.setReferenceSequence(self.seqRef3L, "REFA")
            pA.addTestSequence(self.seqAuth3L, "A")
            pA.doAlign()
            pA.prAlignmentConflicts("A")
            logger.info("Length of reference sequence = %d", len(self.seqRef3L))
            logger.info("Length of    author sequence = %d", len(self.seqAuth3L))
            myAlign = pA.getAlignment("A")
            logger.info("Length   of alignment     = %d", len(myAlign))
            ii = 0
            jj = 0
            for myPr in myAlign:
                if myPr[0] != myPr[1]:
                    logger.debug("Py - conflict at alignment position %d  -  %s - %s", ii, myPr[0], myPr[1])
                    jj += 1
                ii += 1
            logger.info("Conflicts in alignment %d", jj)
            self.assertGreaterEqual(jj, 80)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testAlign3(self):
        """  Consensus alignment for author and reference sequences.
             Returning the name of any sequence that is not part of the consensus.
        """
        logger.info("------------------------ ")
        try:
            pA = PairwiseAlign()
            pA.setVerbose(self.__verbose)
            pA.setReferenceSequence(self.seqRef3L, "REFA")
            pA.addTestSequence(self.seqAuth3L, "A")
            #
            myFails = pA.doAlignConsensus()
            logger.info("Failed sequences = %d", len(myFails))
            pA.prAlignmentConflicts("A")
            logger.info("Length reference sequence = %d", len(self.seqRef3L))
            logger.info("Length    author sequence = %d", len(self.seqAuth3L))
            myAlign = pA.getAlignment("A")
            logger.info("Length   of alignment     = %d", len(myAlign))
            ii = 0
            jj = 0
            for myPr in myAlign:
                if myPr[0] != myPr[1]:
                    logger.debug("Py - conflict position %8d  -  %3s - %3s", ii, myPr[0], myPr[1])
                    jj += 1
                ii += 1
            logger.info("Conflicts in alignment %d", jj)
            self.assertGreaterEqual(jj, 80)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testAlign4(self):
        """  Consensus alignment for reference and 5 test sequences

        """
        logger.info("------------------------")
        try:
            pA = PairwiseAlign()
            pA.setVerbose(self.__verbose)
            pA.setReferenceSequence(self.seqRef3L, "REFA")
            for k, v in self.sTests.items():
                logger.debug("Added test sequence %10s len = %8d", k, len(v))
                pA.addTestSequence(v, k)
            #
            myFails = pA.doAlignConsensus()
            logger.info("Failed sequences = %d", len(myFails))
            logger.info("Length reference sequence = %d", len(self.seqRef3L))
            self.assertEqual(len(myFails), 0)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testAlign5(self):
        """  Consensus alignment for reference and 5 test sequences

        """
        logger.info("------------------------ ")
        try:
            pA = PairwiseAlign()
            pA.setVerbose(self.__verbose)
            pA.setReferenceSequence(self.seqRef3L, "REFA")
            logger.info("Length reference sequence = %d", len(self.seqRef3L))
            #
            for k, v in self.sTests.items():
                logger.info("Added test sequence %10s len = %8d", k, len(v))
                pA.addTestSequence(v, k)
            #
            myFails = pA.doAlignConsensus()
            logger.info("Failed sequences = %d", len(myFails))
            self.assertEqual(len(myFails), 0)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def pairAlignSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PairwiseAlignTests("testAlign1"))
    suiteSelect.addTest(PairwiseAlignTests("testAlign2"))
    suiteSelect.addTest(PairwiseAlignTests("testAlign3"))
    suiteSelect.addTest(PairwiseAlignTests("testAlign4"))
    suiteSelect.addTest(PairwiseAlignTests("testAlign5"))
    return suiteSelect


if __name__ == '__main__':
    #
    mySuite = pairAlignSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
