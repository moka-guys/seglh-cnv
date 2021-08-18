#!/usr/bin/env python

__author__ = "david Brawand"
__copyright__ = ""
__credits__ = ["    "]
__license__ = "                 "
__version__ = "0.0.1"
__maintainer__ = ""
__email__ = ""
__status__ = "develeopment"
__doc__ = """
Summary line.
Extended description of function.

Parameters
----------
    arg1 : file
        Tab delimited segmental duplications (pseudogenes) from Blueprint
    arg2 : file
        RefGene BED12 file (RefSeq IDs)

Returns
-------
    int
        Description of return value
"""

import sys
import re


def getExons(s):
    exon_numbers = []
    for ex in map(lambda x: x.strip(), s.split(',')):
        exon_range = []
        if ex:
            m = re.match(r'(\d+)-(\d+)', ex)
            if m:
                exon_range = list(range(int(m.group(1)), int(m.group(2))+1))
            else:
                exon_range = [int(ex)]
        exon_numbers += exon_range
    return set(exon_numbers)


LEVELS = ["L", "H"]


def main(args):
    # read exons
    annotation = {}
    with open(args[0]) as fh:
        for line in fh:
            f = line.split()
            nm = re.match(r'(NM_\d+)\.\d+', f[2])
            exons = list(map(getExons, f[3:5]))
            for i in range(1, len(exons))[::-1]:
                exons[i-1] = exons[i-1] - exons[i]
            annotation[nm.group(1)] = exons
            # extract segmental duplication
            s = re.match(r'chr([^:]+):(\d+)-(\d+)', f[1])
            print('\t'.join(list(s.groups()) + ['SegDup']))

    # iterate over refGene BED12 and write exon bed
    with open(args[1]) as fh:
        for line in fh:
            f = line.split()
            m = re.match(r'([^_]+_\d+)', f[3])
            tx = m.group(1)
            # extract exon coordinates
            if tx in annotation.keys():
                order = -1 if f[5] == '-' else 1
                e_len = list(map(int, f[10].rstrip(',').split(',')))[::order]
                e_start = list(map(int, f[11].rstrip(',').split(',')))[::order]
                exons = []
                # build exons (UTR exons with size 0)
                for i, s in enumerate(e_start):
                    start = max(int(f[1])+s, int(f[6]))
                    end = min(int(f[1])+s+e_len[i], int(f[7]))
                    e = [f[0], start, end]
                    exons.append(e)
                # extract exons
                for i, l in enumerate(annotation[tx]):
                    name = LEVELS[i]
                    for e in sorted(list(l)):
                        exon_offset =e-1
                        bedline = exons[exon_offset] + [f'{e}{name}']
                        print('\t'.join(map(str, bedline)))

                # remove from annotaion if found
                del annotation[tx]

    print('Not found in BED file:', file=sys.stderr)
    print(annotation, file=sys.stderr)


if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)
