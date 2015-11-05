#!/usr/bin/env python3

"""
This module performs the categorization of the hits found in the rest of the analysis.
"""

import os
from busco_lib.busco_utils import measuring

__author__ = 'Luca Venturini'


def categorize(args, score_dic):

    # TODO: understand the parameters and the function of this method!
    """
    Main and sole function of the module.
    :param args: argparse Namespace
    :param score_dic:
    :type score_dic: dict
    :return:
    """

    # Categorizing genes found in Complete; multi-copy and partial hits
    leng_dic = {}
    sd_dic = {}
    complete = []
    frag = []
    done = []
    cc = []
    fcc = 0
    mcc = []
    unique = []
    if args.mode in ("genome", "report", "hmmer") or args.mode != "OGS":
        temp = os.listdir(os.path.join(args.mainout, 'hmmer_output'))
        files = []
        for i in temp:
            if i.endswith(('.out', '.1', '.2', '.3')):
                files.append(i)
        f = open(os.path.join(args.clade, '/lengths_cutoff'))
        for line in f:
            line = line.strip().split()
            leng_dic[line[0]] = float(line[3])
            sd_dic[line[0]] = float(line[2])
        for entry in files:
            f = open(os.path.join(args.mainout, 'hmmer_output', entry))
            hit_dic = {}
            for line in f:
                if line.startswith('# '):
                    pass
                else:
                    line = line.strip().split()
                    score = float(line[7])
                    group = line[3]
                    prot = line[0]
                    tlen = int(line[2])
                    qlen = int(line[5])
                    if tlen > 30*qlen:
                        pass
                    else:
                        if prot not in hit_dic.keys() and score >= score_dic[group]:
                            hit_dic[prot] = [[int(line[15]), int(line[16])]]
                        elif score >= score_dic[group]:
                            hit_dic[prot].append([int(line[15]), int(line[16])])
            length = measuring(hit_dic)
            try:		# get maximum length of the putative gene in question
                if len(length) == 1:
                    length = length[0]
                else:
                    length = max(length)+1
                sigma = abs(leng_dic[group] - length)/sd_dic[group]
                if sigma <= 2:
                    complete.append(entry)
                    cc.append(group)
                elif sigma > 2:
                    frag.append(entry)
            except (IndexError, TypeError):
                pass
        # check the multi hits
        for entry in complete:
            if entry.endswith('.out'):
                name = entry[:-4]
            else:
                name = entry[:-6]
            if name in done:
                if name not in mcc:
                    mcc.append(name)
                done.append(name)
        for i in cc:
            if i not in mcc:
                unique.append(i)
        for entry in frag:
            if entry.endswith('.out'):
                name = entry[:-4]
            else:
                name = entry[:-6]
            if name not in done and entry not in complete:
                done.append(name)
                fcc += 1

    if args.mode == 'OGS':
        complete = {}
        frag = {}
        done = []
        fcc = []
        temp = os.listdir(os.path.join(args.mainout, 'hmmer_output'))
        files = []
        for i in temp:
            if i.endswith(('.out', '.1', '.2', '.3')):
                files.append(i)
        f = open(os.path.join(args.clade, 'lengths_cutoff'))
        for line in f:
            line = line.strip().split()
            leng_dic[line[0]] = float(line[3])
            sd_dic[line[0]] = float(line[2])
        for entry in files:
            f = open(os.path.join(args.mainout, 'hmmer_output', entry))
            hit_dic = {}
            for line in f:
                if line.startswith('# '):
                    pass
                else:
                    line = line.strip().split()
                    score = float(line[7])
                    group = line[3]
                    prot = line[0]
                    tlen = int(line[2])
                    qlen = int(line[5])
                    if group not in complete:
                        complete[group] = []
                        frag[group] = []
                    if tlen > 30 * qlen:
                        pass
                    else:
                        if prot not in hit_dic.keys() and score >= score_dic[group]:
                            hit_dic[prot] = [[int(line[15]),
                                              int(line[16]),
                                              line[7]]]
                        elif score >= score_dic[group]:
                            hit_dic[prot].append([
                                int(line[15]),
                                int(line[16]),
                                line[7]])
            lengths = measuring(hit_dic)
            try:		# get maximum length of the putative gene in question
                if len(lengths) == 1:
                    length = lengths[0]
                    sigma = abs(leng_dic[group]-length)/sd_dic[group]
                    if sigma <= 2:
                        complete[group].append([
                            list(hit_dic.keys())[lengths.index(length)],
                            hit_dic[list(hit_dic.keys())[lengths.index(length)]][0][2],
                            length])
                    elif sigma > 2:
                        frag[group].append(list(hit_dic.keys())[lengths.index(length)])
                else:
                    for length in lengths:
                        # length=max(lengths)+1
                        sigma = abs(leng_dic[group]-length)/sd_dic[group]
                        if sigma <= 2:
                            complete[group].append([
                                list(hit_dic.keys())[lengths.index(length)],
                                hit_dic[list(hit_dic.keys())[lengths.index(length)]][0][2],
                                length])
                        elif sigma > 2:
                            frag[group].append(list(hit_dic.keys())[lengths.index(length)])
            except (IndexError, KeyError, ValueError, TypeError):
                pass
        # check the multi hits
        for entry in complete:
            if len(complete[entry]) == 0:
                pass
            elif len(complete[entry]) == 1:  # complete
                cc.append(entry)
            elif len(complete[entry]) > 1:
                mcc.append(entry)
        for entry in frag:
            if len(complete[entry]) != 0:
                pass
            elif frag[entry]:
                fcc.append(entry)

    # This is BAD .... how many variables?!?
    return score_dic, leng_dic, sd_dic, complete, frag, done, cc, fcc, mcc, unique
