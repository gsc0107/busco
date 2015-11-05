#!/usr/bin/env python3

"""
This module contains the instructions needed to launch and analyse BLAST results.
"""

import os
import subprocess
from busco_lib.busco_utils import disentangle, gargantua
from collections import deque, namedtuple

__author__ = 'Luca Venturini'


class Tblastn(namedtuple("tblastn",
                         ["is_header", "name", "scaff",
                          "hitstart", "hitend", "posstart", "posend",
                          "e_val", "sizer"])):

    """Namedtuple derived directly from a TBlastN result line."""

    def __new__(cls, line):
        """
        Overload of the parent class. It will format the line and return a
         properly formatted named tuple object.
        :param cls:
        :param line: the line to turn into a named tuple.
        :type line: str

        :return:
        """

        attrs = dict.fromkeys(cls._fields, None)

        if line.startswith("#"):
            attrs["is_header"] = True
        else:
            attrs["is_header"] = False
            line = line.strip().split()
            attrs["name"], attrs["scaff"] = line[:2]
            attrs["sizer"] = int(line[3])
            attrs["hitstart"], attrs["hitend"] = int(line[6]), int(line[7])
            attrs["posstart"], attrs["posend"] = int(line[8]), int(line[9])
            # Reverse order if they are not correct
            attrs["posstart"], attrs["posend"] = sorted(attrs["posstart"], attrs["posend"])
            attrs["e_val"] = float(line[10])
        # Create the named tuple
        return super(Tblastn, cls).__new__(cls, **attrs)


def do_genome_blast(args, flank):

    # TODO: understand the flank parameter

    """
    This function calls and analyses the BLAST results according to the parameters provided by the
    command line options.
    :param args: the argparse Namespace.

    :param flank: a flanking value to add to the blast.
    :type flank: int
    :return:
    """

    print('*** Getting coordinates for candidate regions! ***')
    tblast_file = open('%s_tblastn' % args.abrev)  # open input file
    dic = {}
    coords = {}
    for tblast_line in tblast_file:
        tblast_line = Tblastn(tblast_line)
        if tblast_line.is_header is True:
            continue
        if tblast_line.name not in dic.keys():
            dic[tblast_line.name] = [tblast_line.scaff]
            coords[tblast_line.name] = {}
            coords[tblast_line.name][tblast_line.scaff] = [tblast_line.posstart,
                                                           tblast_line.posend,
                                                           deque([[tblast_line.hitstart,
                                                                   tblast_line.hitend]]),
                                                           tblast_line.sizer]
        # get just the top3 scoring regions
        elif tblast_line.scaff not in dic[tblast_line.name] and len(dic[tblast_line.name]) < 3:
            dic[tblast_line.name].append(tblast_line.scaff)
            coords[tblast_line.name][tblast_line.scaff] = [tblast_line.posstart,
                                                           tblast_line.posend,
                                                           deque([[tblast_line.hitstart,
                                                                   tblast_line.hitend]]),
                                                           tblast_line.sizer]
        # scaffold already checked, now update coordinates
        elif tblast_line.scaff in dic[tblast_line.name] and tblast_line.e_val < args.ev_cut:
            # starts before, and withing 50kb of current position
            if (tblast_line.posstart < coords[tblast_line.name][tblast_line.scaff][0] and
                    coords[tblast_line.name][tblast_line.scaff][0] - tblast_line.posstart <= 50000):

                coords[tblast_line.name][tblast_line.scaff][0] = tblast_line.posstart

                coords[tblast_line.name][tblast_line.scaff][2].append(
                    [tblast_line.hitstart, tblast_line.hitend])

            # ends after and within 50 kbs
            if (tblast_line.posend > coords[tblast_line.name][tblast_line.scaff][1] and
                    tblast_line.posend-coords[tblast_line.name][tblast_line.scaff][1] <= 50000):

                coords[tblast_line.name][tblast_line.scaff][1] = tblast_line.posend
                coords[tblast_line.name][tblast_line.scaff][3] = tblast_line.hitend
                coords[tblast_line.name][tblast_line.scaff][2].append(
                    [tblast_line.hitstart, tblast_line.hitend])

            # starts inside current coordinates
            elif (coords[tblast_line.name][tblast_line.scaff][0] < tblast_line.posstart
                    < coords[tblast_line.name][tblast_line.scaff][1]):
                # if ending inside, just add alignment positions to deque
                if tblast_line.posend < coords[tblast_line.name][tblast_line.scaff][1]:
                    coords[tblast_line.name][tblast_line.scaff][2].append([tblast_line.hitstart,
                                                                           tblast_line.hitend])
                # if ending after current coordinates, extend
                elif tblast_line.posend > coords[tblast_line.name][tblast_line.scaff][1]:
                    coords[tblast_line.name][tblast_line.scaff][2][1] = tblast_line.posend
                    coords[tblast_line.name][tblast_line.scaff][2].append([tblast_line.hitstart,
                                                                           tblast_line.hitend])
    with open('coordinates_{0}'.format(args.abrev), 'w') as out:
        for coord in coords:
            contest = {}
            maxi = 0
            for contig in coords[coord]:
                sizer = disentangle(coords[coord][contig][2])
                if contig not in contest:
                    contest[contig] = 0
                size = gargantua(sizer)
                contest[contig] = size
                if size > maxi:
                    maxi = size
            for contig in contest:
                print(coord, contig,
                      max(0, coords[coord][contig][0] - flank),
                      coords[coord][contig][1] + flank,
                      sep="\t", file=out)
    return dic


def do_transcriptome_blast(args):

    """
    This function calls BLAST and analyses the results when we are analysing a transcriptome.
    :param args: the argparse Namespace.
    :return:
    """

    print('*** Getting coordinates for candidate transcripts! ***')
    f = open('%s_tblastn' % args.abrev)  # open input file
    # transdic = None

    dic = {}
    transdic = {}
    maxi = None
    for i in f:  # get a dictionary of BUSCO matches vs candidate scaffolds
        if i.startswith('#'):
            pass
        else:
            line = i.strip().split()
            name = line[0]
            scaff = line[1]
            e_val = float(line[10])
            leng = int(line[3])
            if name not in dic.keys() and e_val <= args.ev:
                dic[name] = [scaff]
                maxi = leng
                transdic[scaff] = name
            elif (e_val <= args.ev and
                    scaff not in dic[name] and
                    len(dic[name]) < 3 and leng >= 0.7 * maxi):
                dic[name].append(scaff)
                transdic[scaff] = name

    scaff_list = []  # list of unique scaffolds with buscos matches
    for busco in dic:
        for scaff in dic[busco]:
            if scaff not in scaff_list:
                scaff_list.append(scaff)
    print('*** Extracting candidate transcripts! ***')
    f = open(args.genome)
    check = 0
    out = None
    for i in f:
        if i.startswith('>'):
            i = i.strip().split()
            i = i[0][1:]
            if i in scaff_list:
                out = open('%s%s_.temp' % (i, args.abrev), 'w')
                out.write('>%s\n' % i)
                check = 1
            else:
                check = 0
        elif check == 1:
            out.write(i)

    out.close()
    if os.path.exists('%stranslated_proteins' % args.mainout) is False:
        subprocess.call('mkdir %stranslated_proteins' % args.mainout,
                        shell=True)
    files = os.listdir('.')
    lista = []
    for entry in files:
        if entry.endswith(args.abrev+'_.temp'):
            lista.append(entry)

    print('Translating candidate transcripts !')
    for entry in lista:
        command = ["transeq", "-clean", "-frame", "6",
                   "-trim", "-sequence", entry,
                   "-outseq",
                   "{}.fas".format(args.mainout+'translated_proteins/'+entry.split(args.abrev)[0]+'_ts')]
        subprocess.call(command, shell=False)
    f2 = open('%s/scores_cutoff' % args.clade)  # open target scores file
    # Load dictionary of HMM expected scores and full list of groups
    score_dic = {}
    for i in f2:
        i = i.strip().split()
        try:
            score_dic[i[0]] = float(i[1])
        except (ValueError, IndexError, KeyError):
            pass
    totalbuscos = len(list(score_dic.keys()))
    return transdic, totalbuscos


def do_blast_step(args, flank):

    """
    Launcher of the BLAST analysis.
    :param args: The argparse Namespace
    :param flank: flank to be considered when analysing the genome.
    :return:
    """

    print('*** Running tBlastN ***')
    subprocess.call(["makeblastdb", "-in", args.genome,
                     "-dbtype", "nucl", "-out", args.abrev], shell=False)
    subprocess.call(["tblastn", "-num_threads", str(args.cpus),
                     "-query", "{0}/ancestral".format(args.clade),
                     "-db", args.abrev, "-out", "{0}_tblastn".format(args.abrev),
                     "-outfmt", "7"], shell=True)

    transdic = None
    totalbuscos = 0
    genome_dic = {}
    if args.mode in ("genome", "blast"):
        genome_dic = do_genome_blast(args, flank)
    elif args.mode in ("trans", "transcriptome"):
        transdic, totalbuscos = do_transcriptome_blast(args)

    return genome_dic, transdic, totalbuscos
