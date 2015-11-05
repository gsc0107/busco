#!/usr/bin/env python3

"""
This module contains the command line parser for BUSCO.
"""

import argparse
import subprocess
import os
import shutil

__author__ = 'Luca Venturini'


def to_cpu(string):

    """
    Quick function to transform a string into the number of CPU cores to be used.
    It maxes out at the maximum number of cores physically present on the machine.
    :param string:
    :return: counter
    :rtype: int
    """

    return min(os.cpu_count(), max(1, int(string)))


def to_evalue(string):

    """
    Simple function to convert a string into a positive float number (i.e. the evalue)
    :param string: the string to convert
    :return: the inferred number

    :rtype: float
    """

    evalue = float(string)
    if evalue < 0:
        raise TypeError("Evalue too low (less than 0)")
    return evalue


def _parser():

    """
    This function creates the argparse argument parser used by the program.
    :return: the parser.
    :rtype: argparse.ArgumentParser
    """

    usage = 'BUSCO_v1.0.py -in [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] [OTHER OPTIONS]'

    parser = argparse.ArgumentParser(description='''Welcome to the
    Benchmarking set of Universal Single Copy Orthologs (BUSCO).

    For further usage information, please check the README file provided with this distrubution.''',
                                     usage=usage)

    # genome assembly file
    parser.add_argument('-g', '--genome', '-in',
                        metavar='FASTA FILE',
                        type=str,
                        help='''Input file in fasta format.
                        Can be a genome, proteome or transcriptome.
                        Default analysis is run on the genome mode,
                        for other files please specify the mode with (-m [MODE])\n''')
    parser.add_argument('-c', '--cpu', metavar='N', type=to_cpu, default=1,
                        help='Number of threads/cores to use.')	 # Number of available threads
    # Four letter abbreviation for use with genome assembly
    parser.add_argument('-a', '--abrev',
                        '-o', metavar='output', required=True,
                        type=str, help='How to name output and temporary files.')
    parser.add_argument('--ev', '-e', '-ev', metavar='N',
                        type=to_evalue,
                        default=0.01,
                        help='E-value cutoff for BLAST searches. (Default: 0.01)')  # evalue option

    modes = ['all', 'blast',
             'hmmer', 'augustus',
             'parser', 'hmmer+',
             'OGS', 'transcriptome',
             'trans', 'ogs', 'genome']  # valid modes

    parser.add_argument('-m', '--mode', metavar='mode', type=str,
                        choices=modes,
                        default="all",
                        help='''which module to run the analysis to run.
                        Defaults to \'all\'''')
    parser.add_argument('-l', '--clade', '--lineage',
                        metavar='lineage', type=str,
                        help='Which BUSCO lineage to be used.')	 # lineage
    parser.add_argument('-f', action='store_true', default=False,
                        dest='force',
                        help='''Force rewrting of existing files.
                        Must be used when output files with the provided name already exist.''')
    valid_species = [
        'human', 'fly', 'generic',
        'arabidopsis', 'brugia', 'aedes',
        'tribolium', 'schistosoma', 'tetrahymena',
        'galdieria', 'maize', 'toxoplasma',
        'caenorhabditis', '(elegans)', 'aspergillus_fumigatus',
        'aspergillus_nidulans', '(anidulans)', 'aspergillus_oryzae',
        'aspergillus_terreus', 'botrytis_cinerea', 'candida_albicans',
        'candida_guilliermondii', 'candida_tropicalis', 'chaetomium_globosum',
        'coccidioides_immitis', 'coprinus', 'coprinus_cinereus',
        'coyote_tobacco', 'cryptococcus_neoformans_gattii',
        'cryptococcus_neoformans_neoformans_B',
        'cryptococcus_neoformans_neoformans_JEC21',
        '(cryptococcus)',
        'debaryomyces_hansenii',
        'encephalitozoon_cuniculi_GB',
        'eremothecium_gossypii',
        'fusarium_graminearum',
        '(fusarium)',
        'histoplasma_capsulatum',
        '(histoplasma)',
        'kluyveromyces_lactis',
        'laccaria_bicolor',
        'lamprey',
        'leishmania_tarentolae',
        'lodderomyces_elongisporus',
        'magnaporthe_grisea',
        'neurospora_crassa',
        '(neurospora)',
        'phanerochaete_chrysosporium',
        '(pchrysosporium)',
        'pichia_stipitis',
        'rhizopus_oryzae',
        'saccharomyces_cerevisiae_S288C',
        'saccharomyces_cerevisiae_rm11-1a_1',
        '(saccharomyces)',
        'schizosaccharomyces_pombe',
        'thermoanaerobacter_tengcongensis',
        'trichinella',
        'ustilago_maydis',
        '(ustilago)',
        'yarrowia_lipolytica',
        'nasonia',
        'tomato',
        'chlamydomonas',
        'amphimedon',
        'pneumocystis',
        'wheat',
        'chicken']

    parser.add_argument('-sp', '--species', default='generic',
                        metavar='species',
                        type=str,
                        choices=valid_species,
                        help='''Name of existing Augustus species gene finding metaparameters.
                        (Default: generic)''')
    parser.add_argument('-flank', '--flank', '-F',
                        metavar='flanks', type=int,
                        help='''Flanking sequence size for candidate regions.
                        If not provided, flank size is calculated based on genome size
                        with a range from 5 to 20 Kbp.''')
    parser.add_argument('-Z', '--size',
                        metavar='dbsize', type=int,
                        help='HMM library total size (Z). Important if using external datasets')
    parser.add_argument('--long',
                        action='store_true',
                        default=False,
                        dest='long',
                        help='''Optimization mode Augustus self-training (Default: Off).
                        It adds ~20h extra run time, but can improve results for
                        some non-model organisms''')

    return parser


def busco_parser(command_argv):

    """
    Wrapper around the parser. It calls the main function using the command_argv as sole argument.
    :param command_argv: the sys.argv list.
    :type command_argv: list
    :return: the Namespace from argparse.
    """

    parser = _parser()
    args = parser.parse_args(command_argv)  # parse the arguments
    args = check_args(args)
    return args


def check_args(args):

    """
    Function to check that the arguments from the command line pass sanity checks.
    :param args: the argparse Namespace.
    :return:
    """

    args.mainout = os.path.abspath(
        os.path.join(".", "run_{0}".format(args.abrev))
    )

    if os.path.exists(args.mainout) is False and args.abrev is not None:
        os.makedirs(args.mainout)
        # os.system('mkdir %s' % args.mainout)
    else:
        if args.force is False:
            print('''A run with that name already exists!
If are sure you wish to rewrite existing files please use the -f option''')
            raise SystemExit(0)
        else:
            os.removedirs(args.mainout)
            os.makedirs(args.mainout)

    valid_clade_info = {'arthropoda': 102785,
                        'metazoa': 91897,
                        'vertebrata': 143785,
                        'fungi': 174195,
                        'example': 102785,
                        'bacteria': 107114,
                        'eukaryota': 41317}
    args.maxflank = 20000
    # print(args.clade)
    if args.clade is not None:
        # clade = args.clade
        args.clade_name = args.clade.strip('/').split('/')[-1].lower()
        if args.clade_name in valid_clade_info:
            args.Z = valid_clade_info[args.clade_name]
        else:
            print('Using custom lineage data...')
            try:
                args.Z = args.dbsize
            except:
                error = 'Please indicate the size of the custom HMM database'
                raise SystemExit(error)
    else:
        err_msg = """Please indicate the full path to a BUSCO clade:
    Eukaryota,
    Metazoa,
    Arthropoda,
    Vertebrata,
    Fungi.
    Example: -l /path/to/Arthropoda"""
        raise SystemExit(err_msg)
    assert args.clade is not None

    mode_dict = {"all": "genome",
                 "genome": "genome",
                 "ogs": "OGS",
                 "OGS": "OGS",
                 "transcriptome": "trans",
                 "trans": "trans",
                 "hmmer": "hmmer",
                 "augustus": "augustus",
                 "hmmer+": "hmmer",
                 "parser": "parser"}

    args.mode = mode_dict[args.mode]
    #  -------------------------------- Check if necessary programs are acessible ---------------#

    if args.mode in ('genome', 'trans') and shutil.which('tblastn') is None:
        error = """Error: Blast is not accessible from the command-line.
        Please add it to the environment"""
        raise SystemExit(error)

    if args.mode in ('genome', 'trans', 'OGS') and shutil.which('hmmsearch') is None:
        error = """Error: HMMer is not accessible from the command-line.
        Please add it to the environment"""
        raise SystemExit(error)
    elif shutil.which('hmmsearch') is not None:
        hmmer_check = subprocess.check_output('hmmsearch -h', shell=True)
        hmmer_check = hmmer_check.decode('utf-8')
        hmmer_check = hmmer_check.split('\n')[1].split()[2]
        hmmer_check = float(hmmer_check[:3])
        if hmmer_check >= 3.1:
            pass
        else:
            print('Error: HMMer version detected is not unsupported, please use HMMer 3.1+')
            raise SystemExit

    if args.mode == 'genome' and shutil.which('augustus') is None:
        error = '''Error: Augustus is not accessible from the command-line,
        please add it to the environment'''
        raise SystemExit(error)

    if (args.mode == 'genome' and
            os.access(os.environ.get('AUGUSTUS_CONFIG_PATH'), os.W_OK) is False):
        error = """Error: Cannot write to Augustus directory,
        please make sure you have write permissions to {}
        """.format(os.environ.get('AUGUSTUS_CONFIG_PATH'))
        raise SystemError(error)

    if args.mode == 'trans' and shutil.which('transeq') is None:
        error = """Error: EMBOSS transeq is not accessible from the commandline,
                please add it to the environment"""
        raise SystemExit(error)

    return args
