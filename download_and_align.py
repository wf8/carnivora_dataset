#! /usr/bin/python

import csv
from Bio import Entrez
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from time import sleep

# taxon,trnT-trnL_spacer,trnL_intron,trnL-trnF_spacer,psbA-trnH,rpl16

taxa = []
genes = []
accessions = []
records = []

print('Reading in csv file...')
with open('data/accessions.csv', 'rb') as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        if row[0] == 'taxon':
            for i, gene in enumerate(row):
                if i != 0:
                    genes.append(row[i])
                    accessions.append([])
                    records.append([])
        else:
            taxa.append(row[0])
            for i, accession in enumerate(row):
                if i != 0:
                    accessions[i - 1].append(row[i])

Entrez.email = 'freyman@berkeley.edu'

for i, gene in enumerate(genes):
    print('Downloading accessions for ' + gene + '...')
    for j, accession in enumerate(accessions[i]):
        if accession.strip() != '':
            handle = Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id=accession)
            record = SeqIO.read(handle, 'fasta')
            records[i].append(SeqRecord(Seq(str(record.seq), IUPAC.ambiguous_dna), id=taxa[j], description=""))
            handle.close()
            sleep(0.02)
    SeqIO.write(records[i], "sequences_unaligned/" + genes[i] + ".fasta", "fasta")

for i, gene in enumerate(genes):
    print("Aligning " + gene + " with MAFFT...")
    mafft_cline = MafftCommandline(input="sequences_unaligned/" + genes[i] + ".fasta")
    mafft_cline.set_parameter("--auto", True)
    mafft_cline.set_parameter("--adjustdirection", True)
    print(str(mafft_cline))
    stdout, stderr = mafft_cline()

    print("Writing " + gene + " alignment to FASTA file...")
    with open("sequences_aligned/" + genes[i] + ".fasta", "w") as handle:
        handle.write(stdout)

print("Done.")
