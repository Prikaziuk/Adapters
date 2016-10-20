#!/usr/bin/env python

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import re
import argparse
import sys
from collections import defaultdict
import logging
import os.path


parser = argparse.ArgumentParser(description="Script for adapters BLAST and parsing")
parser.add_argument('infile', help="Path to `*.fasta`", type=str)
parser.add_argument('-x', dest='xml_filename', help="Path to blast output `*.xml`", type=str, required=True)
parser.add_argument('-o', dest='outfile', help="Path to parsed output `*.txt`", type=str, required=True)


def check_file_path(path):
    """ Checks if there is such a file"""
    try:
        test_file = open(path, 'a+')
    except IOError:
        raise
    test_file.close()
    return path


def blast(fasta_file):
    """ Opens fasta and blast it with NCBIWWW.qblast()"""
    logging.info("Started to BLAST `{}`".format(os.path.basename(fasta_file)))
    with open(fasta_file, 'r') as fasta:
        fasta = fasta.read()
    # entrez_query = 'Eucarya (taxid:2759)'
    blast_results = NCBIWWW.qblast('blastn', 'nt', fasta, ungapped_alignment=True, expect=0.01)
    # to get more options of NCBIWWW.qblast() uncomment following line
    # help(NCBIWWW.qblast())
    logging.info("Finished to BLAST `{}`".format(os.path.basename(fasta_file)))
    return blast_results


def write_xml(xml_file, blast_results):
    """ Writes blast results to xml"""
    with open(xml_file, 'w') as xml:
        xml.write(blast_results.read())
    blast_results.close()
    logging.info("Blast result saved in `{}`".format(os.path.basename(xml_file)))
    return


def get_and_write_info(xml_file, out_file):
    """
    Parses xml file to BioBlast generator object
    writes : (str) : "query_name \t number_of_hits \t [species_that_hits]" : to outfile
    """
    with open(out_file, 'w') as outfile, open(xml_file, 'r') as xml:
        records = NCBIXML.parse(xml)
        for record in records:
            species = defaultdict(int)
            number_of_hits = 0
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.query_start != 1 or hsp.identities != record.query_length:
                        continue
                    # sbjct = consensus sequence
                    # hsp.identities (?positives) == len of indent nucleotides
                    # hsp.query_start == at which position of querry alignment starts (best == 1)
                    #if hsp.expect > 2e-26:
                        #continue
                    if not re.search('PREDICTED', alignment.title):
                        # alignment.title.split()[0] == gi|51854624|gb|AC148975.3| == hit-id
                        species[(str(alignment.title.split()[1]))] += 1
                        number_of_hits += 1
            species_count = ""
            for s in species:
                species_count += "{}:{}, ".format(s, str(species[s]))
            outfile.write(record.query + '\t' + str(number_of_hits) + '\t' + str(species_count[:-2]) + '\n')
    logging.info("Done!!! Output is here `{}`".format(os.path.abspath(out_file)))
    return


def main():
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', datefmt='%m/%d/%Y %H:%M:%S',
                        stream=sys.stdout, level=logging.DEBUG)
    fasta_file = check_file_path(args.infile)
    xml_file = check_file_path(args.xml_filename)
    out_file = check_file_path(args.outfile)
    blast_results = blast(fasta_file)
    write_xml(xml_file, blast_results)
    get_and_write_info(xml_file, out_file)


if __name__ == "__main__":
    main()
