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
        elif key.start is not None and key.stop is not None:
            if key.start >= key.stop:
                return super().__getitem__(slice(key.start, len(self), key.step)) + super().__getitem__(slice(0, key.stop, key.step))
            else:
                return super().__getitem__(slice(key.start, key.stop, key.step))



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
        print(site_pairs)
        for start, end in site_pairs:
            for site in sites:
                if start == site[0]:
                    start_ovhg = site[1]
                elif end == site[0]:
                    end_ovhg = site[1]
            ovhgs.append((start_ovhg, end_ovhg))
            if end_ovhg < 0:
                if start_ovhg > 0:
                    seqs.append(seq[6:end - end_ovhg])
                else:
                    seqs.append(seq[start:end - end_ovhg])
            elif start_ovhg > 0:
                if end_ovhg < 0:
                    seqs.append(seq[start:end])
                else:
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
            if end == len(seq):
                if ovhgs[0][1] < 0:
                    seqs.append(seq[site_pairs[-1][0]:site_pairs[0][1]] +
                                seq[site_pairs[0][1]:site_pairs[0][1] - ovhgs[0][1]])
                    if ovhgs[-1][0] > 0:
                        seqs[-1] = seq[site_pairs[-1][0] - ovhgs[-1][0]:site_pairs[-1][0]] + seqs[-1]
                elif ovhgs[0][1] > 0:
                    if ovhgs[-1][0] < 0:
                        seqs.append(seq[site_pairs[-1][0]:site_pairs[0][1]])
                    else:
                        seqs.append(seq[site_pairs[-1][0] - ovhgs[-1][0]:site_pairs[-1][0]] +
                                    seq[site_pairs[-1][0]:site_pairs[0][1]])
                else:
                    seqs.append(seq[site_pairs[-1][0]:site_pairs[0][1]])
                del seqs[0]
                ovhgs[len(ovhgs)-1] = (ovhgs[-1][0], ovhgs[0][-1])
                del ovhgs[0]
            elif end != len(seq):
                if end_ovhg < 0:
                    if start_ovhg > 0:
                        seqs.append(seq[6:end - end_ovhg])
                    else:
                        seqs.append(seq[start:end - end_ovhg])
                elif start_ovhg > 0:
                    if end_ovhg < 0:
                        seqs.append(seq[start:end])
                    else:
                        seqs.append(seq[start - start_ovhg:end])
                else:
                    seqs.append(seq[start:end])
        fragments = []
        for i in range(len(seqs)):
            print(fragments)
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
    # my_seq = DNAStrc2("TTTTTCTGCAGTTTAAGCTTTTTTTT", 0, 0)
    # new_seq = my_seq.dscut([PstI, HindIII])
    my_seq = DNAStrc2("AAAAACTGCAGAAAAACTGCAGAAAAA", 0, 0)
    new_seq = my_seq.dscut_circ([HindIII, PstI])
    print(new_seq)
    for i in new_seq:
        print(i)
    # my_seq = str_circ("TTTTTGAATTCTTTTT")
    # print(my_seq[6:6])


