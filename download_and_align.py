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
            #records[i].append(SeqIO.read(handle, 'fasta'))
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

#print('Concatenating...')

#aligned_records = []
#for i, gene in enumerate(genes):
#    aligned_records.append( list(SeqIO.parse("sequences_aligned/" + genes[i] + ".fasta", 'fasta')) )

#def get_unknowns(length):
#    x = ''
#    for i in range(length):
#        x += '?'
#    return x

#final_sequences = []
#for i, taxon in enumerate(taxa):
#    sequence = ''
#    for j, accession in enumerate(accessions):
#        if accession[i].strip() == '':
#            sequence += get_unknowns(len(str(records[j][0].seq)))
#        else:
#            for record in records[j]:
#                if accession[i].strip() in record.description:
#                    sequence += str(record.seq)
#                    break
#    final_sequences.append(sequence)

#final_records = []
#for i in range(len(final_sequences)):
#    final_records.append(SeqRecord(Seq(final_sequences[i], IUPAC.ambiguous_dna), id=taxa[i]))

#print("Making final concatenated alignement file...")
#SeqIO.write(final_records, "carnivora_aligned.fasta", "fasta")

print("Done.")
