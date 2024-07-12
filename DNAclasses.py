from Bio.Seq import Seq
from Bio.Restriction import HindIII
from Bio.Restriction.Restriction import AbstractCut

class DNAStrc(Seq):
    def __init__(self, seq:Seq, linear=True):
        super().__init__(seq)
        self.linear = linear
        #Due to property of Seq class, you cannot pass another compulsory object. Therefore linear=False (make optional)

    def cut(self, enzyme):
        if not isinstance(enzyme, AbstractCut):
            raise ValueError("Invalid restriction enzyme")
        return enzyme.catalyse(self)

    def cut_(self, enzyme):
        if not isinstance(enzyme, AbstractCut):
            raise ValueError("Invalid restriction enzyme")
        else:
            locations = [1] + enzyme.search(self) + [len(self)+1]
            ends = []
            for i in range(len(locations)-1):
                end_pair = [locations[i],locations[i+1]]
                ends.append(end_pair)
            fragments = []
            for start, end in ends:
                fragment = self[start-1:end-1]
                fragments.append(fragment)
            return tuple(fragments)
    
    def cut2(self, enzyme):
        if not isinstance(enzyme, AbstractCut):
            raise ValueError("Invalid restriction enzyme")
        else:
            locations = [1, *enzyme.search(self), len(self)+1]
            ends = zip(locations, locations[1:])
            fragments = []
            for start, end in ends:
                fragment = self[start-1:end-1]
                fragments.append(fragment)
            return tuple(fragments)
        
    def cut_circ(self, enzyme):
        if not isinstance(enzyme, AbstractCut):
            raise ValueError("Invalid restriction enzyme")
        else:
            positions = list(enzyme.search(self, linear = False))
            cut_in_ori = False
            if 1 in positions:
                positions.remove(1)
                cut_in_ori = True
            locations = [1, *positions, len(self)+1]
            print(locations)
            ends = zip(locations, locations[1:])
            fragments = []
            for start, end in ends:
                fragment = self[start-1:end-1]
                fragments.append(fragment)
            if cut_in_ori == True:
                return tuple(fragments)
            else:
                ori:DNAStrc = fragments[-1] + fragments[0]
                del fragments[0]
                del fragments[-1]
                fragments.insert(0, ori)
                return tuple(fragments)

    def cut3(self, enzyme):
        if self.linear:
            return self.cut2(enzyme)
        else:
            return self.cut_circ(enzyme)



if __name__ == "__main__":
    my_seq = DNAStrc("AGCTTTCCCTTCA")
    print(my_seq.cut_circ(HindIII))
