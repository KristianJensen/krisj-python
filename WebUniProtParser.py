#!/usr/bin/env python

#
# Parses UniProt XML entries to a long-list for mapping variants and strains.
# Kristian Jensen, Januar 2015
#

import urllib, urllib2
import re
import sys
from Bio.SeqIO import UniprotIO
import Bio.SeqIO
import StringIO
from multiprocessing.pool import ThreadPool


def printTab(*args):
    print "\t".join((str(_) for _ in args))



def fetch_uniprot(acc):
    "Fetches a Uniprot record by accession number and returns a SeqRecord."
    base_url="http://www.uniprot.org/uniprot/"
    url=base_url+acc+".txt"
    req = urllib2.Request(url)
    response = urllib2.urlopen(req)
    swiss = response.read()
    record = Bio.SeqIO.read(StringIO.StringIO(swiss),"swiss")
    return record

def parse_natural_variants(record):
    for feature in record.features:
        if feature.type == 'VARIANT':
            name = record.name
            start_pos = feature.location.start + 1
            end_pos = feature.location.end
            desc = feature.qualifiers.get('description','')
            original = re.match(r"([ACGT-]+)\s" , desc).group(1)
            variation = re.search(r"->\s([ACGT-]+) ", desc).group(1)
            desc = re.search(r"\((.+)\)",desc).group(1) or desc
            printTab(name,start_pos,end_pos,desc,original,variation)

def parse_uniprot_from_web(acc,parse_function=parse_natural_variants):
    '''Given a UniProt accession number this function will fetch the record (as a SeqRecord)
    from the web and parse it with the function passed in arguments. Default function is parse_natural_variants.
    All errors will be excepted. Use with multiprocessing.pool.ThreadPool to speed up.'''
    try:
        record = fetch_uniprot(acc)
        parse_function(record)
        #sys.stderr.write(acc+" was successfully parsed\n")
    except:
        pass
        #sys.stderr.write(acc+" could NOT be parsed!\n")

if __name__ == "__main__":
    data_dir = "/Users/kristianjensen/Documents/DTU/Speciale/Data/"

    filename = "Ecoli_proteins_with_natural_variations.xml"
    
    p = ThreadPool(1000)
    accs = open("../Data/Ecoli_accessions.txt","r").read().splitlines()
    p.map(parse_uniprot_from_web,accs)
                