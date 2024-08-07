from DNAclasses2 import DNAStrc2
from Bio.Restriction import HindIII, EcoRI, PstI, EcoRV
import pytest


def test_str00():
    DNA = DNAStrc2("CAAGCTTCCC", 0, 0)
    assert str(DNA) == "CAAGCTTCCC\nGTTCGAAGGG"


def test_str02():
    DNA = DNAStrc2("CAAGCTTCCC", 0, 2)
    assert str(DNA) == "CAAGCTTCCC\nGTTCGAAG"


def test_str0n2():
    DNA = DNAStrc2("CAAGCTTCCC", 0, -2)
    assert str(DNA) == "CAAGCTTC\nGTTCGAAGGG"


def test_str20():
    DNA = DNAStrc2("CAAGCTTCCC", 2, 0)
    assert str(DNA) == "  AGCTTCCC\nGTTCGAAGGG"


def test_strn20():
    DNA = DNAStrc2("CAAGCTTCCC", -2, 0)
    assert str(DNA) == "CAAGCTTCCC\n  TCGAAGGG"


def test_str22():
    DNA = DNAStrc2("CAAGCTTCCC", 2, 2)
    assert str(DNA) == "  AGCTTCCC\nGTTCGAAG"


def test_str22n():
    DNA = DNAStrc2("CAAGCTTCCC", 2, -2)
    assert str(DNA) == "  AGCTTC\nGTTCGAAGGG"


def test_strn22():
    DNA = DNAStrc2("CAAGCTTCCC", -2, 2)
    assert str(DNA) == "CAAGCTTCCC\n  TCGAAG"


def test_str2n2n():
    DNA = DNAStrc2("CAAGCTTCCC", -2, -2)
    assert str(DNA) == "CAAGCTTC\n  TCGAAGGG"


def test_dscut():
    DNA = DNAStrc2("CAAGCTTCCC", 0, 0)
    assert DNA.dscut([HindIII]) == (DNAStrc2("CAAGCT"), DNAStrc2("AGCTTCCC"))


def test_dscut2n2n():
    DNA = DNAStrc2("CCCCCAAGCTTCCCAAGCTTCCCCCCC", -2, -2)
    assert DNA.dscut([HindIII]) == (DNAStrc2("CCCCCAAGCT"), DNAStrc2("AGCTTCCCAAGCT"), DNAStrc2("AGCTTCCCCCCC"))    # noqa


def test_dscut2enz():
    DNA = DNAStrc2("TTTTTAAGCTTTTTGAATTCTTTTTT", 0, 0)
    assert DNA.dscut([HindIII, EcoRI]) == (DNAStrc2("TTTTTAAGCT"), DNAStrc2("AGCTTTTTGAATT"), DNAStrc2("AATTCTTTTTT"))  # noqa


def test_dscut_circ_neg():
    DNA = DNAStrc2("TTTTTGAATTCTTTTT", 0, 0)
    assert DNA.dscut_circ([EcoRI]) == (DNAStrc2("AATTCTTTTTTTTTTGAATT"),)


def test_dscut_circ_pos():
    DNA = DNAStrc2("TTTTTCTGCAGTTTTT")
    assert DNA.dscut_circ([PstI]) == (DNAStrc2("TGCAGTTTTTTTTTTCTGCA"),)


def test_dscut_circ_zero():
    DNA = DNAStrc2("TTTTTGATATCTTTTT")
    assert DNA.dscut_circ([EcoRV]) == (DNAStrc2("ATCTTTTTTTTTTGAT"),)

def test_dscut_circ_neg_2enz():
    DNA = DNAStrc2("CCCCCAAGCTTCCCCCAAGCTTCCCCCC")
    assert DNA.dscut_circ([HindIII]) == (DNAStrc2("AGCTTCCCCCAAGCT"), DNAStrc2("AGCTTCCCCCCCCCCCAAGCT"))


def test_dscut_circ_neg_pos_2enz():
    DNA = DNAStrc2("CCCCCAAGCTTGGGGGCTGCAGCCCCCC")
    assert DNA.dscut_circ([PstI, HindIII]) == (DNAStrc2("AGCTTGGGGGCTGCA"), DNAStrc2("TGCAGCCCCCCCCCCCAAGCT"))

def test_dscut_circ_pos_neg_2enz():
    DNA = DNAStrc2("AAAAACTGCAGTTTTTAAGCTTCCCCC")
    assert DNA.dscut_circ([PstI, HindIII]) == (DNAStrc2("TGCAGTTTTTAAGCT"), DNAStrc2("AGCTTCCCCCAAAAACTGCA"))

def test_dscut_circ_pos_2enz():
    DNA = DNAStrc2("AAAAACTGCAGAAAAACTGCAGAAAAA")
    assert DNA.dscut_circ([PstI]) == (DNAStrc2("TGCAGAAAAACTGCA"), DNAStrc2("TGCAGAAAAAAAAAACTGCA"))