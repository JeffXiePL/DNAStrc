from Bio.Seq import _SeqAbstractBaseClass, Seq, SequenceDataAbstractBaseClass
from Bio.Restriction import HindIII
from Bio.Restriction.Restriction import AbstractCut

class DNAStrc2(Seq):
    def __init__(self, seq1:Seq, seq2:Seq, shift:int = 0):
        #these two objects are not the samee, one is a 5' and the other 3'
        super().__init__(seq1)
        super().__init__(seq2)
        self.seq1 = seq1
        self.seq2 = seq2
        self.shift = shift

    def look(self):
        if self.shift >= 0:
            print(self.shift*Seq(" ") + self.seq1 + "\n" + self.seq2)
        else:
            print(self.seq1 + "\n" + -(self.shift)*Seq(" ") + self.seq2)

        #self.seq1 is different to seq1. 
        #self.seq1 is a property of the self, seq1 s(if used as a parameter) is a parameter that the user passes

    def cut(self, enzyme):
        if self.shift >= 0:
            seqs = [(self.shift*Seq(" ") + self.seq1), self.seq2]
        else:
            seqs = [self.seq1, (-(self.shift)*Seq(" ") + self.seq2)]
        if not isinstance(enzyme, AbstractCut):
            raise ValueError("Invalid restriction enzyme")
        else:
            locations1 = [1, *enzyme.search(seqs[0]), len(seqs[0])+1]
            #Add the shift values to enzyme.search output
            locations2 = [1, *enzyme.search(seqs[1]), len(seqs[1])+1]
            print(enzyme.search(seqs[1]))
            #Note: No spaces accounting for right overhangs. Issue?
            ends1 = list(zip(locations1, locations1[1:]))
            ends2 = list(zip(locations2, locations2[1:]))
            fragments1 = []
            fragments2 = []
            for i in range(len(ends1)):
                if ends1[i] == ends2[i]:
                    for start, end in ends1:
                        fragment1 = seqs[0][start-1:end-1]
                        fragments1.append(fragment1)
                    for start, end in ends2:
                        fragment2 = seqs[1][start-1:end-1]
                        fragments2.append(fragment2)
            cut_seqs = [fragments1, fragments2]
            return(cut_seqs)
        #change from 5'-3' to 3'-5' for seq[1]
        #understand @classmethod 
            
if __name__ == "__main__":
    my_seq = DNAStrc2(Seq("CCCAAGCTTCCC"), Seq("CCCAAGCTTCCC"), -3)
    print(my_seq.cut(HindIII))