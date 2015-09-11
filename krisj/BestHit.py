#!/usr/bin/env python

'''Class and functions for performing a Bidirectional/reciprocal best hit search between two
sets of protein sequences. The sets must be indexed by makeblastdb before the functions can run.'''

#
# Reciprocal Best Hit
# Kristian Jensen, Februar 2015
#

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from StringIO import StringIO
import re

def blast_query(query,db,cmd="blastp"):
    "Returns the STDOUT (XML format) of a blastp search of <query> against <db>."
    cline = NcbiblastpCommandline(db=db,cmd=cmd,outfmt=5,remote=False)
    # Argument: outfmt = 5 --- Output format is XML
    stdout, stderr = cline(stdin=query)
    return stdout

def parse_hit_description(description):
    "Extract the gene id from a sequence description"
    match = re.search(r"\[gene=(\w+)\]",description)
    if match is not None:
        gene_name = match.group(1)
    else:
        raise RuntimeError("Invalid sequence description format")
    return gene_name

class BestHit(object):
    '''Reciprocal best hit object. Initialize with the reference fasta file, the reference
    blast db and the query blast db (made from the query organism set of proteins).
    To run the RBH call the run method.'''
    def __init__(self,ref_fasta,ref_db,query_db):
        self.ref_fasta = SeqIO.parse(open(ref_fasta,"rU"),"fasta")
        self.ref_db = ref_db
        self.query_db = query_db
        
    def run(self,extractor_function=parse_hit_description):
        '''Generator that returns the reciprocal best hits. Pass a function as extractor_function that extracts
        the gene name from the sequence description. Default extracts [gene=<name>], as formatted by GenBank
        CDS downloads.'''
        for record in self.ref_fasta:
            ref_gene = extractor_function(record.description)
            blast_result = NCBIXML.read(StringIO(blast_query(record.format("fasta"),self.query_db)))
            try:
                best_hit = extractor_function(blast_result.alignments[0].hit_def)
                hit_seq = str(blast_result.alignments[0].hsps[0].sbjct)
                reciprocal_result = NCBIXML.read(StringIO(blast_query(hit_seq,self.ref_db)))
                reciprocal_best_hit = extractor_function(reciprocal_result.alignments[0].hit_def)
            except IndexError:
                best_hit = None
            
            if reciprocal_best_hit == ref_gene:
                yield (ref_gene, best_hit)
            else:
                yield (ref_gene, None)

if __name__ == "__main__":
    fasta_file = "Data/MG1655_CDS.fasta"
    ref_db = "Data/blast_db/MG1655"
    query_db = "Data/blast_db/DH1"
    
    bh = BestHit(fasta_file,ref_db,query_db)
    for ref, match in bh.run():
        print ref, match
        