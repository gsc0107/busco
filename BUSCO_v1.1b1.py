#!/bin/python

"""
BUSCO - Benchmarking sets of Universal Single-Copy Orthologs.

Copyright (C) 2015 E. Zdobnov lab: F. Simao Neto
<felipe.simao@unige.ch> based on code by R. Waterhouse.

BUSCO is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BUSCO is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Version 1.1c (may/15) - minor bug fixes

                       - lineages may be specified using the full path (e.g. -l /path/to/lineage)
                      - added threading support (signficant speed increase)
                      - added checks for the necessary programs before running BUSCO
"""

import sys
from busco_lib.busco_parser import busco_parser
from busco_lib.busco_blast import do_blast_step
from busco_lib.busco_augustus import do_augustus_step_3
from busco_lib.busco_hmmer import do_hmmer_step4
from busco_lib.busco_categorize import categorize
from busco_lib.busco_summary import summarise
from busco_lib.busco_retraining import retraining


def main():
    """
    Main calling function.
    :return:
    """

    args = busco_parser(sys.argv[1:])
    flank = 5000
    if args.mode in ('genome', 'blast'):  # scalled flanks
        genome_file = open(args.genome)
        size = 0
        for line in genome_file:
            if line.startswith('>'):
                pass
            else:
                size += len(line.strip())
        size /= 1000  # size in mb
        flank = int(size/50)  # proportional flank size
        if flank < 5000:
            flank = 5000
        elif flank > args.maxflank:
            flank = args.maxflank
    genome_dic, transdic, totalbuscos = do_blast_step(args, flank)
    do_augustus_step_3(args)
    transdic, score_dic = do_hmmer_step4(args, transdic)
    (score_dic, leng_dic, sd_dic,
     complete, frag, done, cc, fcc, mcc, unique) = categorize(args, score_dic)
    summarise(args, genome_dic, totalbuscos, mcc, cc, fcc, unique, score_dic, complete, frag)
    retraining(args, totalbuscos, genome_dic, score_dic)

if __name__ == "__main__":
    main()
