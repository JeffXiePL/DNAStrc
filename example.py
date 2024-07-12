
# Historical code
from Bio.Seq import Seq
from Bio.Restriction.Restriction import PstI

# def dscut(self, enzyme):
#     seq = self._data.decode("ASCII")
#     positions = [1, enzyme.search(self), len(seq)]
#     for i in range(len(positions)):
#         positions[i] -= 1
#     seq_start = seq[:positions[0] - enzyme.ovhg]
#     seq_end = seq[positions[-1]:]
#     seqs = []
#     for i in range(len(positions)-1):
#         seqs.append(seq[positions[i]:positions[i+1] - enzyme.ovhg])
#     fragment_start = DNAStrc2(seq_start, self.shiftl, enzyme.ovhg)
#     fragment_end = DNAStrc2(seq_end, enzyme.ovhg, self.shiftr)
#     fragments = []
#     for seq in seqs:
#         fragments.append(DNAStrc2(seq, enzyme.ovhg, enzyme.ovhg))
#     return (fragment_start, *fragments, fragment_end)

# my_seq = Seq("CAAGCTTCCCAAGCTTCCC")
# HindIII.search(my_seq)

# def dscut(self, enzyme):
#     seq = self._data.decode("ASCII")
#     positions = [1, *enzyme.search(self), len(seq)+1]
#     for i in range(len(positions)):
#         positions[i] -= 1
#     seqs = []
#     for start, end in list(zip(positions, positions[1:])):
#         if enzyme.ovhg < 0:
#             seqs.append(seq[start:end - enzyme.ovhg])
#         else:
#             if start == 0:
#                 seqs.append(seq[start:end])
#             else:
#                 seqs.append(seq[start - enzyme.ovhg:end])
#     ovhgs = [
#         self.shiftl,
#         *[
#             enzyme.ovhg
#             for _ in range(len(enzyme.search(self)))
#             ],
#         self.shiftr
#         ]
#     fragments = []
#     for i in range(len(positions)-1):
#         fragments.append(DNAStrc2(seqs[i], ovhgs[i], ovhgs[i+1]))
#     return tuple(fragments)

type(PstI.ovhg)



