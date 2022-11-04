import string
import unittest
from wwpdb.utils.align.alignlib import PseudoMultiAlign  # pylint: disable=no-name-in-module,import-error


class PseudoMultiAlignTest(object):
    def __init__(self):
        self.__authSeq = []
        self.__coorSeq = []
        # self.__authSeq = [
        #     ( "A", "0" ),
        #     ( "A", "1" ),
        #     ( "A", "2" ),
        #     ( "G", "3" ),
        #     ( "C", "4" ),
        #     ( "A", "5" ),
        #     ( "A", "6" ),
        #     ( "G", "7" ),
        #     ( "C", "8" ),
        #     ( "G", "9" ),
        #     ( "G", "10"),
        #     ( "G", "11"),
        #     ( "C", "12"),
        #     ( "C", "13"),
        #     ( "G", "14"),
        #     ( "C", "15"),
        #     ( "A", "16"),
        #     ( "C", "17"),
        #     ( "G", "18"),
        #     ( "C", "19"),
        #     ( "G", "20"),
        #     ( "G", "21"),
        #     ( "C", "22"),
        #     ( "C", "23"),
        #     ( "C", "24"),
        #     ( "G", "25"),
        #     ( "C", "26"),
        #     ( "A", "27"),
        #     ( "A", "28"),
        #     ( "A", "29")
        # ]
        # self.__coorSeq = [
        #     ( "A", "0",  "1", "648"  ),
        #     ( "A", "1",  "1", "649"  ),
        #     ( "A", "2",  "1", "650"  ),
        #     ( "G", "3",  "1", "651"  ),
        #     ( "C", "4",  "1", "652"  ),
        #     ( "A", "5",  "0", "652A" ),
        #     ( "G", "6",  "1", "652C" ),
        #     ( "C", "7",  "1", "652D" ),
        #     ( "G", "8",  "1", "652E" ),
        #     ( "G", "9",  "0", "652F" ),
        #     ( "C", "10", "1", "652S" ),
        #     ( "C", "11", "1", "652T" ) ,
        #     ( "G", "12", "1", "652U" ),
        #     ( "C", "13", "1", "652V" ),
        #     ( "A", "14", "1", "653"  ),
        #     ( "A", "15", "1", "654"  ),
        #     ( "A", "16", "1", "655"  )
        # ]
        # self.__seqA = """GUCAAGAUGGUAAGGGCCCACGGUGGAUGCCUCGGCACCCGAGCCGAUGAAGGACGUGGCUACCUGCGAUAAGCCAGGGG
        # GAGCCGGUAGCGGGCGUGGAUCCCUGGAUGUCCGAAUGGGGGAACCCGGCCGGCGGGAACGCCGGUCACCGCGCUUUUGC
        # GCGGGGGGAACCUGGGGAACUGAAACAUCUCAGUACCCAGAGGAGAGGAAAGAGAAAUCGACUCCCUGAGUAGCGGCGAG
        # CGAAAGGGGACCAGCCUAAACCGUCCGGCUUGUCCGGGCGGGGUCGUGGGGCCCUCGGACACCGAAUCCCCAGCCUAGCC
        # GAAGCUGUUGGGAAGCAGCGCCAGAGAGGGUGAAAGCCCCGUAGGCGAAAGGUGGGGGGAUAGGUGAGGGUACCCGAGUA
        # CCCCGUGGUUCGUGGAGCCAUGGGGGAAUCUGGGCGGACCACCGCCUAAGGCUAAGUACUCCGGGUGACCGAUAGCGCAC
        # CAGUACCGUGAGGGAAAGGUGAAAAGAACCCCGGGAGGGGAGUGAAAUAGAGCCUGAAACCGUGGGCUUACAAGCAGUCA
        # CGGCCCCGCAAGGGGUUGUGGCGUGCCUAUUGAAGCAUGAGCCGGCGACUCACGGUCGUGGGCGAGCUUAAGCCGUUGAG
        # GCGGAGGCGUAGGGAAACCGAGUCCGAACAGGGCGCAAGCGGGCCGCACGCGGCCCGCAAAGUCCGCGGCCGUGGACCCG
        # AAACCGGGCGAGCUAGCCCUGGCCAGGGUGAAGCUGGGGUGAGACCCAGUGGAGGCCCGAACCGGUGGGGGAUGCAAACC
        # CCUCGGAUGAGCUGGGGCUAGGAGUGAAAAGCUAACCGAGCCCGGAGAUAGCUGGUUCUCCCCGAAAUGACUUUAGGGUC
        # AGCCUCAGGCGCUGACUGGGGCCUGUAGAGCACUGAUAGGGCUAGGGGGCCCACCAGCCUACCAAACCCUGUCAAACUCC
        # GAAGGGUCCCAGGUGGAGCCUGGGAGUGAGGGCGCGAGCGAUAACGUCCGCGUCCGAGCGCGGGAACAACCGAGACCGCC
        # AGCUAAGGCCCCCAAGUCUGGGCUAAGUGGUAAAGGAUGUGGCGCCGCGAAGACAGCCAGGAGGUUGGCUUAGAAGCAGC
        # CAUCCUUUAAAGAGUGCGUAAUAGCUCACUGGUCGAGUGGCGCCGCGCCGAAAAUGAUCGGGGCUUAAGCCCAGCGCCGA
        # AGCUGCGGGUCUGGGGGAUGACCCCAGGCGGUAGGGGAGCGUUCCCGAUGCCGAUGAAGGCCGACCCGCGAGGGCGGCUG
        # GAGGUAAGGGAAGUGCGAAUGCCGGCAUGAGUAACGAUAAAGAGGGUGAGAAUCCCUCUCGCCGUAAGCCCAAGGGUUCC
        # UACGCAAUGGUCGUCAGCGUAGGGUUAGGCGGGACCUAAGGUGAAGCCGAAAGGCGUAGCCGAAGGGCAGCCGGUUAAUA
        # UUCCGGCCCUUCCCGCAGGUGCGAUGGGGGGACGCUCUAGGCUAGGGGGACCGGAGCCAUGGACGAGCCCGGCCAGAAGC
        # GCAGGGUGGGAGGUAGGCAAAUCCGCCUCCCAACAAGCUCUGCGUGGUGGGGAAGCCCGUACGGGUGACAACCCCCCGAA
        # GCCAGGGAGCCAAGAAAAGCCUCUAAGCACAACCUGCGGGAACCCGUACCGCAAACCGACACAGGUGGGCGGGUGCAAGA
        # GCACUCAGGCGCGCGGGAGAACCCUCGCCAAGGAACUCUGCAAGUUGGCCCCGUAACUUCGGGAGAAGGGGUGCUCCCUG
        # GGGUGAUGAGCCCCGGGGAGCCGCAGUGAACAGGCUCUGGCGACUGUUUACCAAAAACACAGCUCUCUGCGAACUCGUAA
        # GAGGAGGUAUAGGGAGCGACGCUUGCCCGGUGCCGGAAGGUCAAGGGGAGGGGUGCAAGCCCCGAACCGAAGCCCCGGUG
        # AACGGCGGCCGUAACUAUAACGGUCCUAAGGUAGCGAAAUUCCUUGUCGGGUAAGUUCCGACCUGCACGAAAAGCGUAAC
        # GACCGGAGCGCUGUCUCGGCGAGGGACCCGGUGAAAUUGAACUGGCCGUGAAGAUGCGGCCUACCCGUGGCAGGACGAAA
        # AGACCCCGUGGAGCUUUACUGCAGCCUGGUGUUGGCUCUUGGUCGCGCCUGCGUAGGAUAGGUGGGAGCCUGUGAACCCC
        # CGCCUCCGGGUGGGGGGGAGGCGCCGGUGAAAUACCACCCUGGCGCGGCUGGGGGCCUAACCCUCGGAUGGGGGGACAGC
        # GCUUGGCGGGCAGUUUGACUGGGGCGGUCGCCUCCUAAAAGGUAACGGAGGCGCCCAAAGGUCCCCUCAGGCGGGACGGA
        # AAUCCGCCGGAGAGCGCAAGGGUAGAAGGGGGCCUGACUGCGAGGCCUGCAAGCCGAGCAGGGGCGAAAGCCGGGCCUAG
        # UGAACCGGUGGUCCCGUGUGGAAGGGCCAUCGAUCAACGGAUAAAAGUUACCCCGGGGAUAACAGGCUGAUCUCCCCCGA
        # GCGUCCACAGCGGCGGGGAGGUUUGGCACCUCGAUGUCGGCUCGUCGCAUCCUGGGGCUGAAGAAGGUCCCAAGGGUUGG
        # GCUGUUCGCCCAUUAAAGCGGCACGCGAGCUGGGUUCAGAACGUCGUGAGACAGUUCGGUCUCUAUCCGCCACGGGCGCA
        # GGAGGCUUGAGGGGGGCUCUUCCUAGUACGAGAGGACCGGAAGGGACGCACCUCUGGUUUCCCAGCUGUCCCUCCAGGGG
        # CAUAAGCUGGGUAGCCAUGUGCGGAAGGGAUAACCGCUGAAAGCAUCUAAGCGGGAAGCCCGCCCCAAGAUGAGGCCUCC
        # CACGGCGUCAAGCCGGUAAGGACCCGGGAAGACCACCCGGUGGAUGGGCCGGGGGUGUAAGCGCCGCGAGGCGUUGAGCC
        # GACCGGUCCCAAUCGUCCGAGGUCUUGACCCCUCC"""
        # self.__seqB = "HHHHHHXXXGCGGGAGAACCCCGCCAAGGAACUCUGXXXXXXXXX"

        self.__seq_a = [
            ("MET", "1"),
            ("GLY", "2"),
            ("GLY", "3"),
            ("ILE", "4"),
            ("ARG", "5"),
            ("GLU", "6"),
            ("LYS", "7"),
            ("LYS", "8"),
            ("ALA", "9"),
            ("GLU", "10"),
            ("TYR", "11"),
            ("PHE", "12"),
            ("ALA", "13"),
            ("LYS", "14"),
            ("LEU", "15"),
            ("ARG", "16"),
            ("GLU", "17"),
            ("TYR", "18"),
            ("LEU", "19"),
            ("GLU", "20"),
            ("GLU", "21"),
            ("TYR", "22"),
            ("LYS", "23"),
            ("SER", "24"),
            ("LEU", "25"),
            ("PHE", "26"),
            ("VAL", "27"),
            ("VAL", "28"),
            ("GLY", "29"),
            ("VAL", "30"),
            ("ASP", "31"),
            ("ASN", "32"),
            ("VAL", "33"),
            ("SER", "34"),
            ("SER", "35"),
            ("GLN", "36"),
            ("GLN", "37"),
            ("MET", "38"),
            ("HIS", "39"),
            ("GLU", "40"),
            ("VAL", "41"),
            ("ARG", "42"),
            ("LYS", "43"),
            ("GLU", "44"),
            ("LEU", "45"),
            ("ARG", "46"),
            ("GLY", "47"),
            ("ARG", "48"),
            ("ALA", "49"),
            ("VAL", "50"),
            ("VAL", "51"),
            ("LEU", "52"),
            ("MET", "53"),
            ("GLY", "54"),
            ("LYS", "55"),
            ("ASN", "56"),
            ("THR", "57"),
            ("MET", "58"),
            ("VAL", "59"),
            ("ARG", "60"),
            ("ARG", "61"),
            ("ALA", "62"),
            ("ILE", "63"),
            ("ARG", "64"),
            ("GLY", "65"),
            ("PHE", "66"),
            ("LEU", "67"),
            ("SER", "68"),
            ("ASP", "69"),
            ("LEU", "70"),
            ("PRO", "71"),
            ("ASP", "72"),
            ("PHE", "73"),
            ("GLU", "74"),
            ("LYS", "75"),
            ("LEU", "76"),
            ("LEU", "77"),
            ("PRO", "78"),
            ("PHE", "79"),
            ("VAL", "80"),
            ("LYS", "81"),
            ("GLY", "82"),
            ("ASN", "83"),
            ("VAL", "84"),
            ("GLY", "85"),
            ("PHE", "86"),
            ("VAL", "87"),
            ("PHE", "88"),
            ("THR", "89"),
            ("ASN", "90"),
            ("GLU", "91"),
            ("PRO", "92"),
            ("LEU", "93"),
            ("THR", "94"),
            ("GLU", "95"),
            ("ILE", "96"),
            ("LYS", "97"),
            ("ASN", "98"),
            ("VAL", "99"),
            ("ILE", "100"),
            ("VAL", "101"),
            ("SER", "102"),
            ("ASN", "103"),
            ("ARG", "104"),
            ("VAL", "105"),
            ("ALA", "106"),
            ("ALA", "107"),
            ("PRO", "108"),
            ("ALA", "109"),
            ("ARG", "110"),
            ("ALA", "111"),
            ("GLY", "112"),
            ("ALA", "113"),
            ("VAL", "114"),
            ("ALA", "115"),
            ("PRO", "116"),
            ("GLU", "117"),
            ("ASP", "118"),
            ("ILE", "119"),
            ("TRP", "120"),
            ("VAL", "121"),
            ("ARG", "122"),
            ("ALA", "123"),
            ("VAL", "124"),
            ("ASN", "125"),
            ("THR", "126"),
            ("GLY", "127"),
            ("MET", "128"),
            ("GLU", "129"),
            ("PRO", "130"),
            ("GLY", "131"),
            ("LYS", "132"),
            ("THR", "133"),
            ("SER", "134"),
            ("PHE", "135"),
            ("PHE", "136"),
            ("GLN", "137"),
            ("ALA", "138"),
            ("LEU", "139"),
            ("GLY", "140"),
            ("VAL", "141"),
            ("PRO", "142"),
            ("THR", "143"),
            ("LYS", "144"),
            ("ILE", "145"),
            ("ALA", "146"),
            ("ARG", "147"),
            ("GLY", "148"),
            ("THR", "149"),
            ("ILE", "150"),
            ("GLU", "151"),
            ("ILE", "152"),
            ("VAL", "153"),
            ("SER", "154"),
            ("ASP", "155"),
            ("VAL", "156"),
            ("LYS", "157"),
            ("VAL", "158"),
            ("VAL", "159"),
            ("ASP", "160"),
            ("ALA", "161"),
            ("GLY", "162"),
            ("ASN", "163"),
            ("LYS", "164"),
            ("VAL", "165"),
            ("GLY", "166"),
            ("GLN", "167"),
            ("SER", "168"),
            ("GLU", "169"),
            ("ALA", "170"),
            ("SER", "171"),
            ("LEU", "172"),
            ("LEU", "173"),
            ("ASN", "174"),
            ("LEU", "175"),
            ("LEU", "176"),
            ("ASN", "177"),
            ("ILE", "178"),
            ("SER", "179"),
            ("PRO", "180"),
            ("PHE", "181"),
            ("THR", "182"),
            ("PHE", "183"),
            ("GLY", "184"),
            ("LEU", "185"),
            ("THR", "186"),
            ("VAL", "187"),
            ("VAL", "188"),
            ("GLN", "189"),
            ("VAL", "190"),
            ("TYR", "191"),
            ("ASP", "192"),
            ("ASN", "193"),
            ("GLY", "194"),
            ("GLN", "195"),
            ("VAL", "196"),
            ("PHE", "197"),
            ("PRO", "198"),
            ("SER", "199"),
            ("SER", "200"),
            ("ILE", "201"),
            ("LEU", "202"),
            ("ASP", "203"),
            ("ILE", "204"),
            ("THR", "205"),
            ("ASP", "206"),
            ("GLU", "207"),
            ("GLU", "208"),
            ("LEU", "209"),
            ("VAL", "210"),
            ("SER", "211"),
            ("HIS", "212"),
            ("PHE", "213"),
            ("VAL", "214"),
            ("SER", "215"),
            ("ALA", "216"),
            ("VAL", "217"),
            ("SER", "218"),
            ("THR", "219"),
            ("ILE", "220"),
            ("ALA", "221"),
            ("SER", "222"),
            ("ILE", "223"),
            ("SER", "224"),
            ("LEU", "225"),
            ("ALA", "226"),
            ("ILE", "227"),
            ("GLY", "228"),
            ("TYR", "229"),
            ("PRO", "230"),
            ("THR", "231"),
            ("LEU", "232"),
            ("PRO", "233"),
            ("SER", "234"),
            ("VAL", "235"),
            ("GLY", "236"),
            ("HIS", "237"),
            ("THR", "238"),
            ("LEU", "239"),
            ("ILE", "240"),
            ("ASN", "241"),
            ("ASN", "242"),
            ("TYR", "243"),
            ("LYS", "244"),
            ("ASP", "245"),
            ("LEU", "246"),
            ("LEU", "247"),
            ("ALA", "248"),
            ("VAL", "249"),
            ("ALA", "250"),
            ("ILE", "251"),
            ("ALA", "252"),
            ("ALA", "253"),
            ("SER", "254"),
            ("TYR", "255"),
            ("HIS", "256"),
            ("TYR", "257"),
            ("PRO", "258"),
            ("GLU", "259"),
            ("ILE", "260"),
            ("GLU", "261"),
            ("ASP", "262"),
            ("LEU", "263"),
            ("VAL", "264"),
            ("ASP", "265"),
            ("ARG", "266"),
            ("ILE", "267"),
            ("GLU", "268"),
            ("ASN", "269"),
            ("PRO", "270"),
            ("GLU", "271"),
            ("LYS", "272"),
            ("TYR", "273"),
            ("ALA", "274"),
            ("ALA", "275"),
            ("ALA", "276"),
            ("ALA", "277"),
            ("PRO", "278"),
            ("ALA", "279"),
            ("ALA", "280"),
            ("THR", "281"),
            ("SER", "282"),
            ("ALA", "283"),
            ("ALA", "284"),
            ("SER", "285"),
            ("GLY", "286"),
            ("ASP", "287"),
            ("ALA", "288"),
            ("ALA", "289"),
            ("PRO", "290"),
            ("ALA", "291"),
            ("GLU", "292"),
            ("GLU", "293"),
            ("ALA", "294"),
            ("ALA", "295"),
            ("ALA", "296"),
            ("GLU", "297"),
            ("GLU", "298"),
            ("GLU", "299"),
            ("GLU", "300"),
            ("GLU", "301"),
            ("SER", "302"),
            ("ASP", "303"),
            ("ASP", "304"),
            ("ASP", "305"),
            ("MET", "306"),
            ("GLY", "307"),
            ("PHE", "308"),
            ("GLY", "309"),
            ("LEU", "310"),
            ("PHE", "311"),
            ("ASP", "312"),
        ]

        self.__seq_b = [
            ("GLY", "3", "1"),
            ("ILE", "4", "1"),
            ("ARG", "5", "1"),
            ("GLU", "6", "1"),
            ("LYS", "7", "1"),
            ("LYS", "8", "1"),
            ("ALA", "9", "1"),
            ("GLU", "10", "1"),
            ("TYR", "11", "1"),
            ("PHE", "12", "1"),
            ("ALA", "13", "1"),
            ("LYS", "14", "1"),
            ("LEU", "15", "1"),
            ("ARG", "16", "1"),
            ("GLU", "17", "1"),
            ("TYR", "18", "1"),
            ("LEU", "19", "1"),
            ("GLU", "20", "1"),
            ("GLU", "21", "1"),
            ("TYR", "22", "1"),
            ("LYS", "23", "1"),
            ("SER", "24", "1"),
            ("LEU", "25", "1"),
            ("PHE", "26", "1"),
            ("VAL", "27", "1"),
            ("VAL", "28", "1"),
            ("GLY", "29", "1"),
            ("VAL", "30", "1"),
            ("ASP", "31", "1"),
            ("ASN", "32", "1"),
            ("VAL", "33", "1"),
            ("SER", "34", "1"),
            ("SER", "35", "1"),
            ("GLN", "36", "1"),
            ("GLN", "37", "1"),
            ("MET", "38", "1"),
            ("HIS", "39", "1"),
            ("GLU", "40", "1"),
            ("VAL", "41", "1"),
            ("ARG", "42", "1"),
            ("LYS", "43", "1"),
            ("GLU", "44", "1"),
            ("LEU", "45", "1"),
            ("ARG", "46", "1"),
            ("GLY", "47", "1"),
            ("ARG", "48", "1"),
            ("ALA", "49", "1"),
            ("VAL", "50", "1"),
            ("VAL", "51", "1"),
            ("LEU", "52", "1"),
            ("MET", "53", "1"),
            ("GLY", "54", "1"),
            ("LYS", "55", "1"),
            ("ASN", "56", "1"),
            ("THR", "57", "1"),
            ("MET", "58", "1"),
            ("VAL", "59", "1"),
            ("ARG", "60", "1"),
            ("ARG", "61", "1"),
            ("ALA", "62", "1"),
            ("ILE", "63", "1"),
            ("ARG", "64", "1"),
            ("GLY", "65", "1"),
            ("PHE", "66", "1"),
            ("LEU", "67", "1"),
            ("SER", "68", "1"),
            ("ASP", "69", "1"),
            ("LEU", "70", "1"),
            ("PRO", "71", "1"),
            ("ASP", "72", "1"),
            ("PHE", "73", "1"),
            ("GLU", "74", "1"),
            ("LYS", "75", "1"),
            ("LEU", "76", "1"),
            ("LEU", "77", "1"),
            ("PRO", "78", "1"),
            ("PHE", "79", "1"),
            ("VAL", "80", "1"),
            ("LYS", "81", "1"),
            ("GLY", "82", "1"),
            ("ASN", "83", "1"),
            ("VAL", "84", "1"),
            ("GLY", "85", "1"),
            ("PHE", "86", "1"),
            ("VAL", "87", "1"),
            ("PHE", "88", "1"),
            ("THR", "89", "1"),
            ("ASN", "90", "1"),
            ("GLU", "91", "1"),
            ("PRO", "92", "1"),
            ("LEU", "93", "1"),
            ("THR", "94", "1"),
            ("GLU", "95", "1"),
            ("ILE", "96", "1"),
            ("LYS", "97", "1"),
            ("ASN", "98", "1"),
            ("VAL", "99", "1"),
            ("ILE", "100", "1"),
            ("VAL", "101", "1"),
            ("SER", "102", "1"),
            ("ASN", "103", "1"),
            ("ARG", "104", "1"),
            ("VAL", "105", "1"),
            ("ALA", "106", "1"),
            ("ALA", "107", "1"),
            ("GLY", "184", "0"),
            ("LEU", "185", "1"),
            ("THR", "186", "1"),
            ("VAL", "187", "1"),
            ("VAL", "188", "1"),
            ("GLN", "189", "1"),
            ("VAL", "190", "1"),
            ("TYR", "191", "1"),
            ("ASP", "192", "1"),
            ("ASN", "193", "1"),
            ("GLY", "194", "1"),
            ("GLN", "195", "1"),
            ("VAL", "196", "1"),
            ("PHE", "197", "1"),
            ("PRO", "198", "1"),
            ("SER", "199", "1"),
        ]

    def runTest1(self):
        self.__authSeq = []
        for idx, tupL in enumerate(self.__seq_a):
            self.__authSeq.append((tupL[0], str(idx)))
        #
        self.__coorSeq = []
        for idx, tupL in enumerate(self.__seq_b):
            self.__coorSeq.append((tupL[0], str(idx), tupL[2], tupL[1]))
        #
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
        pA.setAuthSequence(self.__toList(self.__seqA))  # pylint: disable=no-member
        pA.addAlignSequence(self.__toList(self.__seqB))  # pylint: disable=no-member
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
        # pM.runTest2()


def __main():  # pragma: no cover
    pM = PseudoMultiAlignTest()
    pM.runTest1()
    # pM.runTest2()


if __name__ == "__main__":  # pragma: no cover
    __main()
