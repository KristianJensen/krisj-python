#!/usr/bin/env python

#
# Parses UniProt XML entries to a long-list for mapping variants and strains.
# Kristian Jensen, Januar 2015
#

from Bio.SeqIO import UniprotIO


def printTab(*args):
    print "\t".join((str(_) for _ in args))

data_dir = "/Users/kristianjensen/Documents/DTU/Speciale/Data/"

filename = "Ecoli_proteins_with_natural_variations.xml"

def parse_record(record):
    for feature in record.features:
        if feature.type == 'sequence variant':
            name = record.name
            start_pos = feature.location.start + 1
            end_pos = feature.location.end
            desc = feature.qualifiers.get('description','').split(";")[0]
            original = feature.qualifiers.get('original','')
            variation = feature.qualifiers.get('variation','')
            printTab(name,start_pos,end_pos,desc,original,variation)


if __name__ == "__main__":
    with open(data_dir+filename,"r") as xml_file:
        for record in UniprotIO.UniprotIterator(xml_file):
            parse_record(record)
                