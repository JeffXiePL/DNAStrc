#%%
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import FeatureLocation

"""
for index, record in enumerate(SeqIO.parse("SPAPB1A10.09.gb", "genbank")):
    print(
        "index %i, ID = %s, length %i, with %i features"
        % (index, record.id, len(record.seq), len(record.features))
    )
record: SeqRecord
"""

#Defining a sequence
record1 = SeqRecord(
    seq = Seq("AAATTTTTAAAAGGGCGCGCAATCATAGGCATAA"),
    id="111.1",
    name="testseq",
    description="sequence for testing, DNA"
)

#Adding a translation
CDS = record1.seq[2:14]
CDS1 = record1.seq[18:27]
#Self note: How do I tell python that something is a method?
CDS_trans = CDS.translate()
CDS1_trans = CDS1.translate()
join = Seq("").join([CDS_trans,CDS1_trans])
record1.features.append(SeqFeature(
                         location = CompoundLocation([FeatureLocation(2,14),FeatureLocation(18,27)]),
                         strand=+1,
                         type = "CDS",
                         qualifiers={"translation":join}
                         ))
#Find ways to write two separate locations on a single seq feature edit



#Writing into .gb file
record1.annotations['molecule_type'] = 'DNA'
SeqIO.write(record1, "example2.gb","genbank")
#%%
