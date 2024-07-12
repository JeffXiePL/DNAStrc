from Bio.Seq import Seq

my_seq = Seq('ACGT')

print(my_seq._data)
print(my_seq._data.decode("ASCII"))