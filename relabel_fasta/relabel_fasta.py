#!/usr/bin/env python

DESCRIPTION = \
"""Replace FASTA labels with new labels <PREFIX>1, <PREFIX>2,
<PREFIX>3 ... (<PREFIX> is provided by the user via the command
line).

Can be used to label OTUs as OTU_1, OTU_2 etc.

This reimplements the functionality of the fasta_number.py utility
from https://drive5.com/python/fasta_number_py.html
"""

import argparse

def relabel_fasta(fp,prefix,include_size=False):
    """
    Relabel sequence records in a FASTA file

    Arguments:
      fp (File): file-like object opened for reading
        input FASTA data from
      prefix (str): prefix to use in new labels
      include_size (bool): if True then copy
        'size=...' records into new labels (default
        is not to copy the size)

    Yields: updated lines from the input FASTA.
    """
    # Iterate over lines in file
    nlabel = 0
    for line in fp:
        # Strip trailing newlines
        line = line.rstrip('\n')
        if not line:
            # Skip blank lines
            continue
        elif line.startswith('>'):
            # Deal with start of a sequence record
            nlabel += 1
            label = line[1:].strip()
            if include_size:
                # Extract size from the label
                try:
                    size = filter(
                        lambda x: x.startswith("size="),
                        label.split(';'))[0]
                except Exception as ex:
                    raise Exception("Couldn't locate 'size' in "
                                    "label: %s" % label)
                yield ">%s%d;%s" % (args.prefix,
                                    nlabel,
                                    size)
            else:
                yield ">%s%d" % (args.prefix,
                                 nlabel)
        else:
            # Echo the line to output
            yield line

if __name__ == "__main__":
    # Set up command line parser
    p = argparse.ArgumentParser(description=DESCRIPTION)
    p.add_argument("--needsize",
                   action="store_true",
                   help="include the size as part of the "
                   "output label ('size=...' must be present "
                   "in the input FASTA labels). Output labels "
                   "will be '<PREFIX><NUMBER>;size=<SIZE>'")
    p.add_argument("--nosize",
                   action="store_true",
                   help="don't include the size as part of "
                   "the output label (this is the default)")
    p.add_argument("fasta",
                   metavar="FASTA",
                   help="input FASTA file")
    p.add_argument("prefix",
                   metavar="PREFIX",
                   help="prefix to use for labels in output")
    # Process command line
    args = p.parse_args()
    # Relabel FASTA
    with open(args.fasta,'rU') as fasta:
        for line in relabel_fasta(fasta,
                                  args.prefix,
                                  include_size=args.needsize):
            print line

