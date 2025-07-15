from Bio import SeqIO
from Bio.Seq import Seq

with open("../trnL/blast_seeds_output/trnL.fasta") as infile, open("../trnL/blast_seeds_output/trnL_with_Ns.fasta", "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        record.seq = Seq("N" * 10 + str(record.seq) + "N" * 10)
        SeqIO.write(record, outfile, "fasta")

