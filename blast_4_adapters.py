from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import re
import time

startTime = time.time()

# fasta = SeqIO.parse('blast_train.fasta', 'fasta')
fasta = SeqIO.parse('blast_train.fasta', 'fasta')

# print(record.) # id, name, description

# ungapped_alignment=True
#XML = open("bact_blast.xml", "w")
outfile = open('parsed_blast.csv', 'w')

for query in fasta:
    # what will happen with blank blast?
    # outfile.write(query.description)
    xml_filename = query.id + '.xml'
    blast = NCBIWWW.qblast('blastn', 'nt', query.seq, ungapped_alignment=True, expect=0.01)
    #blast = NCBIWWW.qblast("blastn", "nt", record.seq)
    xml = open(xml_filename, 'w')
    xml.write(blast.read())
    blast.close()
    xml.close()

    xml = open(xml_filename, 'r')
    records = NCBIXML.parse(xml)

    for record in records:
        species = set()
        number_of_hits = 0
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                # sbjct = consensus sequence
                # hsp.identities (?positives) == len of indent nucleotides
                # hsp.query_start == at which position of querry alignment starts (best == 1)
                if hsp.expect > 2e-26:
                    continue
                if not re.search('PREDICTED', alignment.title):
                    # alignment.title.split()[0] == gi|51854624|gb|AC148975.3| == hit-id
                    # may be defaultdict
                    species.add(str(alignment.title.split()[1]))
                    number_of_hits += 1
        species = list(species)
        species.sort()
        outfile.write(query.description + '\t' + str(number_of_hits) + '\t' + str(species) +'\n')
    xml.close()
outfile.close()

print("Finished. Elapsed time, minutes: " + str((time.time() - startTime) / 60))

#help(NCBIWWW)


# ftp://ftp.ncbi.nlm.nih.gov/pub/factsheets/HowTo_BLASTGuide.pdf
# http://www.ncbi.nlm.nih.gov/BLAST/Doc/urlapi.html

# maybe try 'refseq' instead of 'nt' (GenBank)
# expect == e-value == 0.01 (def=10) ?1.54773e-26
# entrez_query == 'eucaryotes (taxid:2759)'
# hitlist_size=50
# gapcosts=None
# nucl_penalty=-4 (def=-3), nucl_reward=1 (def=1) mismatch/match
# help(NCBIWWW)

