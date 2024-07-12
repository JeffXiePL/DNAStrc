from Bio.Seq import Seq, complement
from Bio.Restriction import HindIII, PstI, EcoRV, EcoRI
from Bio.Restriction.Restriction import AbstractCut
from typing import SupportsIndex


class str_circ(str):
    def __new__(cls, data):
        return super().__new__(cls, data)

    def __getitem__(self, key: SupportsIndex | slice) -> str:
        if isinstance(key, int):
            return super().__getitem__(key)
        else:
            new_start = key.start if key.start else 0
            new_stop = key.stop - 1 if key.stop else 0
            if new_start > new_stop:
                self.data = self.data[key.start:] + self.data[ :key.stop]
            return super().__getitem__(slice(new_start, new_stop, key.step))


class DNAStrc2(Seq):
    def __init__(self, seq, shiftl: int = 0, shiftr: int = 0):
        super().__init__(seq)
        self.shiftl = shiftl
        self.shiftr = shiftr

    def __str__(self):
        seqstr = self._data.decode("ASCII")
        seqc = complement(seqstr)
        if self.shiftr >= 0:
            if self.shiftl < 0:
                return str(
                    seqstr 
                    + "\n"
                    + Seq(" ") * -self.shiftl + seqc[-self.shiftl:len(self) - self.shiftr]
                    )
            else:
                return str(
                    Seq(" ") * self.shiftl + seqstr[self.shiftl:] 
                    + "\n" 
                    + seqc[:len(self) - self.shiftr]
                    )
        else:
            if self.shiftl < 0:
                return str(
                    seqstr[:self.shiftr] 
                    + "\n" 
                    + Seq(" ") * -self.shiftl + seqc[-self.shiftl:]
                    )
            else:
                return str(
                    Seq(" ") * self.shiftl + seqstr[self.shiftl: self.shiftr] 
                    + "\n" 
                    + seqc
                    )

    def dscut(self, enzymes: list[AbstractCut]):
        seq = self._data.decode("ASCII")
        sites = []
        for enzyme in enzymes:
            sites += [(pos-1, enzyme.ovhg) for pos in enzyme.search(Seq(seq))]
        sites = sorted([(0, self.shiftl), *sites, (len(seq), self.shiftr)])
        site_pairs = [(sites[i][0], sites[i+1][0]) for i in range(len(sites)-1)]
        seqs = []
        ovhgs = []
        for start, end in site_pairs:
            for site in sites:
                if start == site[0]:
                    start_ovhg = site[1]
                elif end == site[0]:
                    end_ovhg = site[1]
            ovhgs.append((start_ovhg, end_ovhg))
            if end_ovhg < 0:
                seqs.append(seq[start:end - end_ovhg])
            elif start_ovhg > 0:
                seqs.append(seq[start - start_ovhg:end]) 
            else:
                seqs.append(seq[start:end])
        fragments = []
        for i in range(len(seqs)):
            fragments.append(DNAStrc2(seqs[i], ovhgs[i][0], ovhgs[i][1]))
        return tuple(fragments)

    def dscut_circ(self, enzymes: list[AbstractCut]):
        seq = str_circ(self._data.decode("ASCII"))
        sites = []
        for enzyme in enzymes:
            sites += [(pos-1, enzyme.ovhg) for pos in enzyme.search(Seq(seq))]
        sites = sorted([(0, self.shiftl), *sites, (len(seq), self.shiftr)])
        site_pairs = [(sites[i][0], sites[i+1][0]) for i in range(len(sites)-1)]
        seqs = []
        ovhgs = []
        for start, end in site_pairs:
            for site in sites:
                if start == site[0]:
                    start_ovhg = site[1]
                elif end == site[0]:
                    end_ovhg = site[1]
            ovhgs.append((start_ovhg, end_ovhg))
            if end_ovhg < 0:
                seqs.append(seq[start:end - end_ovhg])
            elif start_ovhg > 0:
                seqs.append(seq[start - start_ovhg:end]) 
            else:
                seqs.append(seq[start:end])
        fragments = []
        for i in range(len(seqs)):
            fragments.append(DNAStrc2(seqs[i], ovhgs[i][0], ovhgs[i][1]))
        return tuple(fragments)



if __name__ == "__main__":
    # my_seq = DNAStrc2("CCCCCAAGCTTCCCAAGCTTCCCCCCC", 2, 2)
    # new_seq = my_seq.dscut([HindIII])
    # my_seq = DNAStrc2("TTTTTCTGCAGTTTCTGCAGTTTTTT", 2, -2)
    # new_seq = my_seq.dscut(PstI)
    # my_seq = DNAStrc2("TTTTTGATATCTTTGATATCTTTTTT", -2, -2)
    # new_seq = my_seq.dscut(EcoRV)
    # my_seq = DNAStrc2("TTTTTAAGCTTTTTGAATTCTTTTTT", 2, 2)
    # new_seq = my_seq.dscut([EcoRI, HindIII])
    # my_seq = DNAStrc2("TTTTTAAGCTTTTTCTGCAGTTTTTT", 0, 0)
    # new_seq = my_seq.dscut([PstI, HindIII])
    my_seq = DNAStrc2("TTTTTGAATTCTTTTT", 0, 0)
    new_seq = my_seq.dscut_circ([EcoRI])
    print(new_seq)
    for i in new_seq:
        print(i)

