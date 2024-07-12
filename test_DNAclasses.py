from DNAclasses import DNAStrc
from Bio.Restriction import HindIII
import pytest

def test_assignment():
    DNA = DNAStrc("TTTTTTAAGCTTGGGGG")
    assert str(DNA) == "TTTTTTAAGCTTGGGGG"

def test_cut():
    DNA = DNAStrc("TTTTTTAAGCTTGGGGG")
    assert DNA.cut(HindIII) == (DNAStrc('TTTTTTA'), DNAStrc('AGCTTGGGGG'))

def test_cut2():
    DNA = DNAStrc("TTTTTTAAGCTTGGGGG")
    assert DNA.cut(HindIII) == DNA.cut2(HindIII)

def test_cut_circ():
    DNA = DNAStrc("CTTTTTTTTAAGCTTGGGGGAAGCTTCCCCCCCAAG")
    assert DNA.cut_circ(HindIII) == (DNAStrc('AGCTTTTTTTTA'), DNAStrc('AGCTTGGGGGA'), DNAStrc('AGCTTCCCCCCCA'))

def test_validate_types():
    DNA = DNAStrc("TTTTTTAAGCTTGGGGG")
    with pytest.raises(ValueError):
        DNA.cut("not_enzyme")

#write classes to test the circ function and the linear functions