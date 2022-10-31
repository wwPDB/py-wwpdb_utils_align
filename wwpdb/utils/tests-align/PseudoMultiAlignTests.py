import string
import unittest
from wwpdb.utils.align.alignlib import PseudoMultiAlign  # pylint: disable=no-name-in-module


class PseudoMultiAlignTest(object):
    def __init__(self):
        self.__authSeq = [
            ( "A", "0" ),
            ( "A", "1" ),
            ( "A", "2" ),
            ( "G", "3" ),
            ( "C", "4" ),
            ( "A", "5" ),
            ( "A", "6" ),
            ( "G", "7" ),
            ( "C", "8" ),
            ( "G", "9" ),
            ( "G", "10"),
            ( "G", "11"),
            ( "C", "12"),
            ( "C", "13"),
            ( "G", "14"),
            ( "C", "15"),
            ( "A", "16"),
            ( "C", "17"),
            ( "G", "18"),
            ( "C", "19"),
            ( "G", "20"),
            ( "G", "21"),
            ( "C", "22"),
            ( "C", "23"),
            ( "C", "24"),
            ( "G", "25"),
            ( "C", "26"),
            ( "A", "27"),
            ( "A", "28"),
            ( "A", "29")
        ]
        self.__coorSeq = [
            ( "A", "0",  "1", "648"  ),
            ( "A", "1",  "1", "649"  ),
            ( "A", "2",  "1", "650"  ),
            ( "G", "3",  "1", "651"  ),
            ( "C", "4",  "1", "652"  ),
            ( "A", "5",  "0", "652A" ),
            ( "G", "6",  "1", "652C" ),
            ( "C", "7",  "1", "652D" ),
            ( "G", "8",  "1", "652E" ),
            ( "G", "9",  "0", "652F" ),
            ( "C", "10", "1", "652S" ),
            ( "C", "11", "1", "652T" ) ,
            ( "G", "12", "1", "652U" ),
            ( "C", "13", "1", "652V" ),
            ( "A", "14", "1", "653"  ),
            ( "A", "15", "1", "654"  ),
            ( "A", "16", "1", "655"  )
        ]
        self.__seqA = """GUCAAGAUGGUAAGGGCCCACGGUGGAUGCCUCGGCACCCGAGCCGAUGAAGGACGUGGCUACCUGCGAUAAGCCAGGGG
        GAGCCGGUAGCGGGCGUGGAUCCCUGGAUGUCCGAAUGGGGGAACCCGGCCGGCGGGAACGCCGGUCACCGCGCUUUUGC
        GCGGGGGGAACCUGGGGAACUGAAACAUCUCAGUACCCAGAGGAGAGGAAAGAGAAAUCGACUCCCUGAGUAGCGGCGAG
        CGAAAGGGGACCAGCCUAAACCGUCCGGCUUGUCCGGGCGGGGUCGUGGGGCCCUCGGACACCGAAUCCCCAGCCUAGCC
        GAAGCUGUUGGGAAGCAGCGCCAGAGAGGGUGAAAGCCCCGUAGGCGAAAGGUGGGGGGAUAGGUGAGGGUACCCGAGUA
        CCCCGUGGUUCGUGGAGCCAUGGGGGAAUCUGGGCGGACCACCGCCUAAGGCUAAGUACUCCGGGUGACCGAUAGCGCAC
        CAGUACCGUGAGGGAAAGGUGAAAAGAACCCCGGGAGGGGAGUGAAAUAGAGCCUGAAACCGUGGGCUUACAAGCAGUCA
        CGGCCCCGCAAGGGGUUGUGGCGUGCCUAUUGAAGCAUGAGCCGGCGACUCACGGUCGUGGGCGAGCUUAAGCCGUUGAG
        GCGGAGGCGUAGGGAAACCGAGUCCGAACAGGGCGCAAGCGGGCCGCACGCGGCCCGCAAAGUCCGCGGCCGUGGACCCG
        AAACCGGGCGAGCUAGCCCUGGCCAGGGUGAAGCUGGGGUGAGACCCAGUGGAGGCCCGAACCGGUGGGGGAUGCAAACC
        CCUCGGAUGAGCUGGGGCUAGGAGUGAAAAGCUAACCGAGCCCGGAGAUAGCUGGUUCUCCCCGAAAUGACUUUAGGGUC
        AGCCUCAGGCGCUGACUGGGGCCUGUAGAGCACUGAUAGGGCUAGGGGGCCCACCAGCCUACCAAACCCUGUCAAACUCC
        GAAGGGUCCCAGGUGGAGCCUGGGAGUGAGGGCGCGAGCGAUAACGUCCGCGUCCGAGCGCGGGAACAACCGAGACCGCC
        AGCUAAGGCCCCCAAGUCUGGGCUAAGUGGUAAAGGAUGUGGCGCCGCGAAGACAGCCAGGAGGUUGGCUUAGAAGCAGC
        CAUCCUUUAAAGAGUGCGUAAUAGCUCACUGGUCGAGUGGCGCCGCGCCGAAAAUGAUCGGGGCUUAAGCCCAGCGCCGA
        AGCUGCGGGUCUGGGGGAUGACCCCAGGCGGUAGGGGAGCGUUCCCGAUGCCGAUGAAGGCCGACCCGCGAGGGCGGCUG
        GAGGUAAGGGAAGUGCGAAUGCCGGCAUGAGUAACGAUAAAGAGGGUGAGAAUCCCUCUCGCCGUAAGCCCAAGGGUUCC
        UACGCAAUGGUCGUCAGCGUAGGGUUAGGCGGGACCUAAGGUGAAGCCGAAAGGCGUAGCCGAAGGGCAGCCGGUUAAUA
        UUCCGGCCCUUCCCGCAGGUGCGAUGGGGGGACGCUCUAGGCUAGGGGGACCGGAGCCAUGGACGAGCCCGGCCAGAAGC
        GCAGGGUGGGAGGUAGGCAAAUCCGCCUCCCAACAAGCUCUGCGUGGUGGGGAAGCCCGUACGGGUGACAACCCCCCGAA
        GCCAGGGAGCCAAGAAAAGCCUCUAAGCACAACCUGCGGGAACCCGUACCGCAAACCGACACAGGUGGGCGGGUGCAAGA
        GCACUCAGGCGCGCGGGAGAACCCUCGCCAAGGAACUCUGCAAGUUGGCCCCGUAACUUCGGGAGAAGGGGUGCUCCCUG
        GGGUGAUGAGCCCCGGGGAGCCGCAGUGAACAGGCUCUGGCGACUGUUUACCAAAAACACAGCUCUCUGCGAACUCGUAA
        GAGGAGGUAUAGGGAGCGACGCUUGCCCGGUGCCGGAAGGUCAAGGGGAGGGGUGCAAGCCCCGAACCGAAGCCCCGGUG
        AACGGCGGCCGUAACUAUAACGGUCCUAAGGUAGCGAAAUUCCUUGUCGGGUAAGUUCCGACCUGCACGAAAAGCGUAAC
        GACCGGAGCGCUGUCUCGGCGAGGGACCCGGUGAAAUUGAACUGGCCGUGAAGAUGCGGCCUACCCGUGGCAGGACGAAA
        AGACCCCGUGGAGCUUUACUGCAGCCUGGUGUUGGCUCUUGGUCGCGCCUGCGUAGGAUAGGUGGGAGCCUGUGAACCCC
        CGCCUCCGGGUGGGGGGGAGGCGCCGGUGAAAUACCACCCUGGCGCGGCUGGGGGCCUAACCCUCGGAUGGGGGGACAGC
        GCUUGGCGGGCAGUUUGACUGGGGCGGUCGCCUCCUAAAAGGUAACGGAGGCGCCCAAAGGUCCCCUCAGGCGGGACGGA
        AAUCCGCCGGAGAGCGCAAGGGUAGAAGGGGGCCUGACUGCGAGGCCUGCAAGCCGAGCAGGGGCGAAAGCCGGGCCUAG
        UGAACCGGUGGUCCCGUGUGGAAGGGCCAUCGAUCAACGGAUAAAAGUUACCCCGGGGAUAACAGGCUGAUCUCCCCCGA
        GCGUCCACAGCGGCGGGGAGGUUUGGCACCUCGAUGUCGGCUCGUCGCAUCCUGGGGCUGAAGAAGGUCCCAAGGGUUGG
        GCUGUUCGCCCAUUAAAGCGGCACGCGAGCUGGGUUCAGAACGUCGUGAGACAGUUCGGUCUCUAUCCGCCACGGGCGCA
        GGAGGCUUGAGGGGGGCUCUUCCUAGUACGAGAGGACCGGAAGGGACGCACCUCUGGUUUCCCAGCUGUCCCUCCAGGGG
        CAUAAGCUGGGUAGCCAUGUGCGGAAGGGAUAACCGCUGAAAGCAUCUAAGCGGGAAGCCCGCCCCAAGAUGAGGCCUCC
        CACGGCGUCAAGCCGGUAAGGACCCGGGAAGACCACCCGGUGGAUGGGCCGGGGGUGUAAGCGCCGCGAGGCGUUGAGCC
        GACCGGUCCCAAUCGUCCGAGGUCUUGACCCCUCC"""
        self.__seqB = "HHHHHHXXXGCGGGAGAACCCCGCCAAGGAACUCUGXXXXXXXXX"

    def runTest1(self):
        pA = PseudoMultiAlign()
        pA.setAuthSequence(self.__authSeq)
        pA.addAlignSequence(self.__coorSeq)
        alignIndexList = pA.getAlignIndices()
        alignSeqList = pA.getAlignSequences()
        #
        for idx, alignIdx in enumerate(alignIndexList):
            sym = "---"
            if (alignSeqList[idx][0] != ".") and (alignSeqList[idx][0] == alignSeqList[idx][1]):
                sym = "==="
            #
            print("%3d %3s (%3d) %s %3s (%3d)" % (idx, alignSeqList[idx][0], alignIdx[0], sym, alignSeqList[idx][1], alignIdx[1]))
        #

    def runTest2(self):
        pA = PseudoMultiAlign()
        pA.setRefScore()
        pA.setAuthSequence(self.__toList(self.__seqA))
        pA.addAlignSequence(self.__toList(self.__seqB))
        alignIndexList = pA.getAlignIndices()
        alignSeqList = pA.getAlignSequences()
        #
        print("")
        for idx, alignIdx in enumerate(alignIndexList):
            sym = "---"
            if (alignSeqList[idx][0] != ".") and (alignSeqList[idx][0] == alignSeqList[idx][1]):
                sym = "==="
            #
            print("%3d %3s (%3d) %s %3s (%3d)" % (idx, alignSeqList[idx][0], alignIdx[0], sym, alignSeqList[idx][1], alignIdx[1]))
        #

    def __toList(self, strIn):
        sL = []
        count = 0
        for ss in strIn:
            if ss in string.whitespace:
                continue
            #
            sL.append((ss, str(count)))
            count += 1
        #
        return sL


class PseudoMultiAlignUnitTest(unittest.TestCase):
    def testRun(self):
        pM = PseudoMultiAlignTest()
        pM.runTest1()
        pM.runTest2()


def __main():  # pragma: no cover
    pM = PseudoMultiAlignTest()
    pM.runTest1()
    pM.runTest2()


if __name__ == "__main__":  # pragma: no cover
    __main()
