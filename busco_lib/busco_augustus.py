#!/usr/bin/env python3

"""
This module contains the functions needed to perform Augustus.
"""

import time
import os
import subprocess
import multiprocessing

__author__ = 'Luca Venturini'


def process_command(command):
    """
    Simple wrapper to call a command with subprocess and record/print the time.
    :param command: the command to launch, as a string.
    :type command: str
    :return:
    """

    print(time.ctime(), "Started command: {0}".format(command))
    subprocess.call(command)
    print(time.ctime(), "Finished with: {0}".format(command))


def do_augustus_step_3(args):

    """
    Third step of the augustus modelling.

    :param args: the argparse Namespace.

    :return:
    """

    # Extract candidate contigs/scaffolds from genome assembly
    # (necessary because augustus doesn't handle multi-fasta files
    # when running on a specific target region)
    if args.mode in ("genome", "augustus"):
        # target_species=species_list[0]
        print('*** pre-Augustus scaffold extraction ***')
        coord = open('coordinates_%s' % args['abrev'])
        dic = {}
        scaff_list = []
        for i in coord:
            i = i.strip().split()
            if len(i) != 2:
                dic[i[0]] = [i[1], i[2], i[3]]
                if i[1] not in scaff_list:
                    scaff_list.append(i[1])
        f = open(args['genome'])
        check = 0
        out = None
        for i in f:
            if i.startswith('>'):
                i = i.split()
                i = i[0][1:]
                if i in scaff_list:
                    out = open('%s%s_.temp' % (i, args['abrev']), 'w')
                    out.write('>%s\n' % i)
                    check = 1
                else:
                    check = 0
            elif check == 1:
                assert out is not None
                out.write(i)
        out.close()
    # ################

    # ################

    # Step-4
    # Augustus search on candidate regions using the pre-built Block profiles (msa2prfl.pl)
        print('*** Running Augustus prediction ***')
        if os.path.exists('{}augustus'.format(args.mainout)) is False:
            os.makedirs("{}augustus".format(args.mainout))

        f = open('coordinates_%s' % args['abrev'])
        dic = {}
        for i in f:
            i = i.strip().split('\t')
            name = i[0]
            if name not in dic:
                dic[name] = [[i[1], i[2], i[3]]]  # scaffold,start and end
            elif name in dic:
                dic[name].append([i[1], i[2], i[3]])
        strings = []
        for i in dic:
            if len(dic[i]) > 1:
                for z in range(0, len(dic[i])):
                    command = ["augustus",
                               "--proteinprofile={clade}/{prot_profile}".format(
                                   clade=args.clade,
                                   prot_profile='prfl/'+i+'.prfl'),
                               "--predictionStart={0}".format(dic[i][z][1]),
                               "--predictionEnd={0}".format(dic[i][z][2]),
                               "--species={species}".format(species=args.target_species),
                               "\"{scaffold}\"".format(scaffold=dic[i][z][0] + args['abrev'] + '_.temp'),
                               ">",
                               "{output}".format(output=args.mainout+'augustus/'+i+'.out.'+str(z+1)),
                               "2>/dev/null"
                               ]
                    command = " ".join(command)
                    strings.append(command)
            else:
                command = [
                    "augustus",
                    "--proteinprofile={clade}/{prot_profile}".format(clade=args.clade,
                                                                     prot_profile='prfl/'+i+'.prfl'),
                    "--predictionStart={start_coord}".format(
                        start_coord=dic[i][0][1]),
                    "--predictionEnd={end_coord}".format(end_coord=dic[i][0][2]),
                    "--species={species}".format(species=args.target_species),
                    "\"{scaffold}\" > {output} 2>/dev/null".format(
                        scaffold=dic[i][0][0]+args['abrev']+'_.temp')
                ]
                command = " ".join(command)
                strings.append(command)

        pool = multiprocessing.Pool(processes=args.cpus)

        for command in strings:
            pool.apply_async(process_command, args=(command,))

        pool.close()
        pool.join()
