#!/bin/python

# BUSCO - Benchmarking sets of Universal Single-Copy Orthologs.

# Copyright (C) 2015 E. Zdobnov lab: F. Simao Neto
# <felipe.simao@unige.ch> based on code by R. Waterhouse.

# BUSCO is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# BUSCO is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Version 1.1c (may/15) - minor bug fixes

#                        - lineages may be specified using the full path (e.g. -l /path/to/lineage)
#                       - added threading support (signficant speed increase)
#                       - added checks for the necessary programs before running BUSCO


# -------------------------------------------------------------------------------#

import os
import shutil
import argparse
from collections import deque
import time
import queue
import threading
import subprocess


start_time = time.time()

# ------------------------------ Argument parser START ----------------------------------------#

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
parser.add_argument('-c', '--cpu', metavar='N', type=str,
                    help='Number of threads/cores to use.')	 # Number of available threads
# Four letter abbreviation for use with genome assembly
parser.add_argument('-a', '--abrev',
                    '-o', metavar='output',
                    type=str, help='How to name output and temporary files.')
parser.add_argument('--ev', '-e', '-ev', metavar='N',
                    type=float,
                    help='E-value cutoff for BLAST searches. (Default: 0.01)')  # evalue option
parser.add_argument('-m', '--mode', metavar='mode', type=str,
                    help='''which module to run the analysis to run, valid modes are:
                    - 'all' (genome assembly),
                    - 'OGS' (gene set / proteome)
                    - 'Trans' (transcriptome).
                    Defaults to \'all\'''')
parser.add_argument('-l', '--clade', '--lineage',
                    metavar='lineage', type=str,
                    help='Which BUSCO lineage to be used.')	 # lineage
parser.add_argument('-f', action='store_true', default=False,
                    dest='force',
                    help='''Force rewrting of existing files.
                    Must be used when output files with the provided name already exist.''')
parser.add_argument('-sp', '--species', default='generic',
                    metavar='species',
                    type=str,
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


args = vars(parser.parse_args())  # parse the arguments
# print(args)
mainout = os.path.abspath(
    os.path.join(".", "run_{0}".format(args["abrev"]))
)

    # './run_%s/' % args['abrev']  # final output directory
if os.path.exists(mainout) is False and args['abrev'] is not None:
    os.makedirs(mainout)
    # os.system('mkdir %s' % mainout)
else:
    if args['force'] is False:
        print('''A run with that name already exists!
        If are sure you wish to rewrite existing files please use the -f option''')
        raise SystemExit
    else:
        os.removedirs(mainout)
        os.makedirs(mainout)

target_species = 'generic'
if args['species'] != 'generic':
    if args['species'] not in [
      'human',
      'fly',
      'generic',
      'arabidopsis',
      'brugia',
      'aedes',
      'tribolium',
      'schistosoma',
      'tetrahymena',
      'galdieria',
      'maize',
      'toxoplasma',
      'caenorhabditis',
      '(elegans)',
      'aspergillus_fumigatus',
      'aspergillus_nidulans',
      '(anidulans)',
      'aspergillus_oryzae',
      'aspergillus_terreus',
      'botrytis_cinerea',
      'candida_albicans',
      'candida_guilliermondii',
      'candida_tropicalis',
      'chaetomium_globosum',
      'coccidioides_immitis',
      'coprinus',
      'coprinus_cinereus',
      'coyote_tobacco',
      'cryptococcus_neoformans_gattii',
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
      'chicken']:
        print('Invalid gene predictor species parameters,'
              'please check the file \'Possible_species.txt\'')
        raise SystemExit
    else:
        target_species = args['species']
        print(target_species)

ev_cut = 0.01  # default e-value cuttof
try:
    if args['ev'] != ev_cut and args['ev'] is not None:
        print('WARNING: You are using a custom e-value cutoff')
        ev_cut = args['ev']
except:
    pass


valid_clade_info = {'arthropoda': 102785,
                    'metazoa': 91897,
                    'vertebrata': 143785,
                    'fungi': 174195,
                    'example': 102785,
                    'bacteria': 107114,
                    'eukaryota': 41317}
maxflank = 20000
# print(args['clade'])
clade = None
try:
    if args['clade'] is not None:
        clade = args['clade']
        clade_name = clade.strip('/').split('/')[-1].lower()
        print(clade_name)
        if clade_name in valid_clade_info:
            Z = valid_clade_info[clade_name]
        else:
            print('Using custom lineage data...')
            try:
                Z = args['dbsize']
            except:
                print('Please indicate the size of the custom HMM database using the (-Z integer)')
                raise SystemExit
except:
    err_msg = """Please indicate the full path to a BUSCO clade:
Eukaryota,
Metazoa,
Arthropoda,
Vertebrata,
Fungi.
Example: -l /path/to/Arthropoda"""
    print(err_msg)
    raise SystemExit
assert clade is not None

cpus = 1  # 1 core default
try:
    if args['cpu'] != cpus and args['cpu'] is not None:
        cpus = args['cpu']
except:
    pass

modes = ['all', 'blast',
         'hmmer', 'augustus',
         'parser', 'hmmer+',
         'OGS', 'transcriptome',
         'trans', 'ogs', 'genome']  # valid modes
mode = 'genome'  # unless otherwise specified, run on all (mode for genome assembly)
try:
    if args['mode'] is not None and args['mode'] in modes:
        mode = args['mode']
        if mode == 'ogs':
            mode = 'OGS'
        elif mode == 'all' or mode == 'genome':
            mode = 'genome'
        elif mode == 'transcriptome':
            mode = 'trans'
except:
    print('''Error: Unknown mode specified * %s *,
    please check the documentation for valid modes.''' % args['mode'])
    raise SystemExit

# -------------------------------- Check if necessary programs are acessible --------------------#


def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE) == 0

if mode in ('genome', 'trans') and cmd_exists('tblastn') is False:
    print('Error: Blast is not accessible from the command-line, please add it to the environment')
    raise SystemExit


if mode in ('genome', 'trans', 'OGS') and cmd_exists('hmmsearch') is False:
    print('Error: HMMer is not accessible from the command-line, please add it to the environment')
    raise SystemExit
elif cmd_exists('hmmsearch') is True:
    hmmer_check = subprocess.check_output('hmmsearch -h', shell=True)
    hmmer_check = hmmer_check.decode('utf-8')
    hmmer_check = hmmer_check.split('\n')[1].split()[2]
    hmmer_check = float(hmmer_check[:3])
    if hmmer_check >= 3.1:
        pass
    else:
        print('Error: HMMer version detected is not unsupported, please use HMMer 3.1+')
        raise SystemExit

if mode=='genome' and cmd_exists('augustus')==False:
    print('''Error: Augustus is not accessible from the command-line,
    please add it to the environment''')
    raise SystemExit

if mode == 'genome' and os.access(os.environ.get('AUGUSTUS_CONFIG_PATH'), os.W_OK) is False:
    print('Error: Cannot write to Augustus directory,',
          'please make sure you have write permissions to {}'.format(
        os.environ.get('AUGUSTUS_CONFIG_PATH')))

if mode == 'trans' and cmd_exists('transeq') is False:
    print('Error: EMBOSS transeq is not accessible from the commandline,'
          'please add it to the environment')
    raise SystemExit

# -----------------------------------------------checks over ------------------------------------#

# get clade specific parameters.
if mode in ('genome', 'blast'):  # scalled flanks
    f = open(args['genome'])
    size = 0
    for line in f:
        if line.startswith('>'):
            pass
        else:
            size += len(line.strip())
    size /= 1000  # size in mb
    flank = int(size/50)  # proportional flank size
    if flank < 5000:
        flank = 5000
    elif flank > maxflank:
        flank = maxflank

# ---------------------------------- Argument parser END ----------------------------------------#


def measuring(nested):
    if isinstance(nested, str):
        return '0'
    scaffolds = list(nested.keys())
    total_len = None
    if len(nested) == 1:
        total_len = [0]
        for hit in nested[scaffolds[0]]:
            total_len[0] += hit[1] - hit[0]
    elif len(nested) > 1:
        total_len = [0] * len(nested)
        for entry in range(0, len(scaffolds)):
            for hit in nested[scaffolds[entry]]:
                total_len[entry] += hit[1] - hit[0]
    try:
        return total_len
    except:
        pass


def extract(path, group):
    count = 0
    if group.endswith(('.1', '.2', '.3')):
        f = open('%saugustus/%s' % (path,group))
        out = open('%saugustus_proteins/%s.fas.%s' % (path,group[:-6],group[-1]), 'w')
    else:
        f = open('%saugustus/%s.out' % (path,group))
        out = open('%saugustus_proteins/%s.fas' % (path,group),'w')
    check = 0
    while True:
        line = f.readline()
        if not line:
            break
        if line.startswith('# start gene'):
            line = f.readline()
            line = line.split()
            places = [line[0], line[3], line[4]]
        elif line.startswith('# protein'):
            line = line.strip().split('[')
            count += 1
            out.write('>p%s[%s:%s-%s]\n' % (count,
                                            places[0],
                                            places[1],
                                            places[2]))
            if line[1][-1] == ']':
                line[1] = line[1][:-1]
            out.write(line[1])
            check = 1
        else:
            if line.startswith('# end'):
                check = 0
                out.write('\n')
            elif check == 1:
                line = line.split()[1]
                if line[-1] == ']':
                    line = line[:-1]
                out.write(line)
    out.close()


def disentangle(deck):
    structure = deque([deck.popleft()])
    try:
        while 1:
            temp = deck.popleft()
            start, end = temp[0], temp[1]
            for i in range(0, len(structure)):
                ds, de = structure[i][0], structure[i][1]
                if start < ds and end < ds:  # fully before
                    if i == 0:  # first entry, just appendleft
                        structure.appendleft(temp)
                        break
                    else:
                        new = structure[0:i];new.append(temp)
                        for z in range(i,len(structure)):
                            new.append(structure[z])
                        break
                # end overlaps inside, but the start is before
                elif start < ds < end < de:
                    structure[i][0]=start
                    break
                # start overlaps inside, but the end is after
                elif ds < start < de < end:
                    structure[i][1]=end
                    break
                elif start>de and end>de:  # fully after
                    # only if its the last entry can it be safely added to structure
                    if i == len(structure)-1:
                        structure.append(temp)
                # current structure is found fully inside the current entry
                elif start<ds and end>de:
                 structure[i]=temp
    except:
        return structure


def gargantua(deck):
    total = 0
    for entry in deck:
        total += entry[1] - entry[0]
    return total


def shrink(number):
    number *= 100
    if number >= 10:
        number = str(number)[:2]
    elif 0 < number < 10:
        number = str(number)[:3]
    return number

# ---------------------------BLAST steps START -------------------------------------------#

# Make a blast database and run tblastn
if mode in ('genome', 'blast', 'trans'):
    print('*** Running tBlastN ***')
    subprocess.call(["makeblastdb", "-in", args["genome"],
                     "-dbtype", "nucl", "-out", args["abrev"]], shell=False)
    subprocess.call(["tblastn", "-num_threads", str(cpus),
                     "-query", "{0}/ancestral".format(clade),
                     "-db", args["abrev"], "-out", "{0}_tblastn".format(args["abrev"]),
                     "-outfmt", "7"], shell=True)

# Get coordinates for a genome analysis
if mode in ("genome", "blast"):
    print('*** Getting coordinates for candidate regions! ***')
    f = open('%s_tblastn' % args['abrev'])  # open input file
    out = open('coordinates_%s' % args['abrev'], 'w')  # open Coordinates output file
    dic = {}
    coords = {}
    for i in f:
        if i.startswith('# '):
            pass
        else:
            line = i.strip().split()
            name = line[0]
            scaff = line[1]
            hitstart = int(line[6])
            hitend = int(line[7])
            postart = int(line[8])
            posend = int(line[9])
            e_val = float(line[10])
            sizer = int(line[3])
            if posend < postart:  # for minus-strand genes, invert coordinates for convenience
                temp = posend
                posend = postart
                postart = temp
            if name not in dic.keys():  # create new entry in dictionary for current BUSCO
                dic[name] = [scaff]
                coords[name] = {}
                coords[name][scaff] = [postart,
                                       posend,
                                       deque([[hitstart,hitend]]),
                                       sizer]
            # get just the top3 scoring regions
            elif scaff not in dic[name] and len(dic[name]) < 3:
                dic[name].append(scaff)
                coords[name][scaff] = [postart,
                                       posend,
                                       deque([[hitstart,hitend]]),
                                       sizer]
            # scaffold already checked, now update coordinates
            elif scaff in dic[name] and e_val < ev_cut:
                # starts before, and withing 50kb of current position
                if postart < coords[name][scaff][0] and coords[name][scaff][0] - postart <= 50000:
                        coords[name][scaff][0] = postart
                        coords[name][scaff][2].append([hitstart,
                                                       hitend])
                # ends after and within 50 kbs
                if posend > coords[name][scaff][1] and posend-coords[name][scaff][1] <= 50000:
                    coords[name][scaff][1] = posend
                    coords[name][scaff][3] = hitend
                    coords[name][scaff][2].append([hitstart, hitend])
                # starts inside current coordinates
                elif coords[name][scaff][0] < postart < coords[name][scaff][1]:
                    # if ending inside, just add alignment positions to deque
                    if posend < coords[name][scaff][1]:
                        coords[name][scaff][2].append([hitstart, hitend])
                    # if ending after current coordinates, extend
                    elif posend > coords[name][scaff][1]:
                        coords[name][scaff][2][1] = posend
                        coords[name][scaff][2].append([hitstart, hitend])

    for i in coords:
        contest = {}
        maxi = 0
        for contig in coords[i]:
            sizer = disentangle(coords[i][contig][2])
            if contig not in contest:
                contest[contig] = 0
            size=gargantua(sizer)
            contest[contig] = size
            if size > maxi:
                maxi = size
        for contig in contest:
            out.write('%s\t%s\t%s\t%s\n' % (i,
                                            contig,
                                            max(0, coords[i][contig][0] - flank),
                                            coords[i][contig][1] + flank))
    out.close()

# Get coordinates, candidate regions and translate sequences (transcriptome analysis)
transdic = None
totalbuscos = 0
if mode in ("transcriptome", "trans"):
    print('*** Getting coordinates for candidate transcripts! ***')
    f = open('%s_tblastn' % args['abrev'])  # open input file
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
            if name not in dic.keys() and e_val <= ev_cut:
                dic[name] = [scaff]
                maxi = leng
                transdic[scaff] = name
            elif (e_val <= ev_cut and
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
    f = open(args['genome'])
    check = 0
    out = None
    for i in f:
        if i.startswith('>'):
            i = i.strip().split()
            i=i[0][1:]
            if i in scaff_list:
                out=open('%s%s_.temp' % (i, args['abrev']), 'w')
                out.write('>%s\n' % i)
                check = 1
            else:
                check = 0
        elif check == 1:
            out.write(i)

    out.close()
    if os.path.exists('%stranslated_proteins' % mainout) is False:
        subprocess.call('mkdir %stranslated_proteins' % mainout,
                        shell=True)
    files=os.listdir('.')
    lista=[]
    for entry in files:
        if entry.endswith(args['abrev']+'_.temp'):
            lista.append(entry)

    print('Translating candidate transcripts !')
    for entry in lista:
        command = ["transeq", "-clean", "-frame", "6",
                   "-trim", "-sequence", entry,
                   "-outseq", "{}.fas".format(
                        mainout+'translated_proteins/'+entry.split(args['abrev'])[0]+'_ts'),
                   ]
        subprocess.call(command, shell=False)
    f2 = open('%s/scores_cutoff' % clade)  # open target scores file
    # Load dictionary of HMM expected scores and full list of groups
    score_dic = {}
    for i in f2:
        i = i.strip().split()
        try:
            score_dic[i[0]] = float(i[1])  # float [1] = mean value; [2] = minimum value
        except:
            pass
    totalbuscos = len(list(score_dic.keys()))

# ---------------------------AUGUSTUS steps START -------------------------------------------#

exitFlag = 0


class myThread(threading.Thread):

    def __init__(self, threadID, name, q):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q

    def run(self):
        print("Starting " + self.name)
        process_data(self.name, self.q)  # ,self.alpha)
        print("Exiting " + self.name)


def process_data(threadName, q):
    while not exitFlag:
        queueLock.acquire()
        if not workQueue.empty():
            data = q.get()
            queueLock.release()
            os.system('%s' % (data))
        else:
            queueLock.release()
        time.sleep(1)


# Step-3
# Extract candidate contigs/scaffolds from genome assembly
# (necessary because augustus doesn't handle multi-fasta files
# when running on a specific target region)
if mode in ("genome", "augustus"):
    # target_species=species_list[0]
    print('*** pre-Augustus scaffold extraction ***')
    coord = open('coordinates_%s' % args['abrev'])
    dic = {}
    scaff_list = []
    for i in coord:
        i = i.strip().split()
        if len(i)!=2:
            dic[i[0]] = [i[1],i[2],i[3]]
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
                out = open('%s%s_.temp' % (i,args['abrev']),'w')
                out.write('>%s\n' % (i))
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
    if os.path.exists('{}augustus'.format(mainout)) is False:
        os.makedirs("{}augustus".format(mainout))

    f = open('coordinates_%s' % args['abrev'])
    dic = {}
    for i in f:
        i = i.strip().split('\t')
        name = i[0]
        if name not in dic:
            dic[name] = [[i[1],i[2],i[3]]]  # scaffold,start and end
        elif name in dic:
            dic[name].append([i[1], i[2], i[3]])
    strings = []
    for i in dic:
        if len(dic[i])>1:
            for z in range(0,len(dic[i])):
                command = ["augustus",
                           "--proteinprofile={clade}/{prot_profile}".format(
                               clade=clade,
                               prot_profile='prfl/'+i+'.prfl'),
                           "--predictionStart={0}".format(dic[i][z][1]),
                           "--predictionEnd={0}".format(dic[i][z][2]),
                           "--species={species}".format(species=target_species),
                           "\"{scaffold}\"".format(scaffold=dic[i][z][0]+args['abrev']+'_.temp'),
                           ">",
                           "{output}".format(output=mainout+'augustus/'+i+'.out.'+str(z+1)),
                           "2>/dev/null"
                           ]
                command = " ".join(command)
                strings.append(command)
        else:
            # TODO: reprise from here
            command = [
                "augustus",
                "--proteinprofile={clade}/{prot_profile}".format(clade=clade,
                                                                 prot_profile='prfl/'+i+'.prfl'),
                "--predictionStart={start_coord}".format(
                    start_coord=dic[i][0][1]),
                "--predictionEnd={end_coord}".format(end_coord=dic[i][0][2]),
                "--species={species}".format(species=target_species),
                "\"{scaffold}\" > {output} 2>/dev/null".format(
                    scaffold=dic[i][0][0]+args['abrev']+'_.temp')
            ]
            command = " ".join(command)
            # command='augustus --proteinprofile=%(clade)s/%(prot_profile)s --predictionStart=%(start_coord)s --predictionEnd=%(end_coord)s --species=%(species)s \"%(scaffold)s\" > %(output)s 2>/dev/null' % \
            #         {'prot_profile':'prfl/'+i+'.prfl',
            #          'start_coord':dic[i][0][1],
            #          'end_coord':dic[i][0][2],
            #          'clade':clade,
            #          'species':target_species,
            #          'scaffold':dic[i][0][0]+args['abrev']+'_.temp',
            #          'output':mainout+'augustus/'+i+'.out'}
            strings.append(command)

    threadList = []
    for i in range(int(cpus)):
        threadList.append("Thread-%s" % str(i+1))

    nameList = list(dic.keys())
    queueLock = threading.Lock()
    workQueue = queue.Queue(len(strings))
    threads = []
    threadID = 1

    # Create new threads
    for tName in threadList:
            thread = myThread(threadID, tName, workQueue)
            thread.start()
            threads.append(thread)
            threadID += 1

        # Fill the queue
    queueLock.acquire()
    for word in strings:
            workQueue.put(word)
    queueLock.release()

    # Wait for queue to empty
    while not workQueue.empty():
            pass
    # Notify threads it's time to exit
    exitFlag = 1

# Wait for all threads to complete
    for t in threads:
            t.join()
    print("Exiting Main Thread")
    exitFlag = 0


# ---------------------------AUGUSTUS steps END -------------------------------------------#

# ---------------------------HMMER steps START -------------------------------------------#

if mode == 'genome' or mode == 'hmmer':  # should be augustus
    # STEP-1 EXTRACT AUGUSTUS PROTEINS
    print('*** Extracting predicted proteins ***')
    files=os.listdir(mainout+'augustus')
    count=0;check=0
    for i in files:
        os.system('sed -i \'1,3d\' %saugustus/%s' % (mainout,i))
    if os.path.exists(mainout+'augustus_proteins')==False:
        os.system('mkdir %saugustus_proteins' % mainout)

    for i in files:
        f=open(mainout+'augustus/'+i)
        if i.endswith('.out'):
            out=open('%saugustus_proteins/%s.fas' % (mainout,i[:-4]),'w')
        elif i.endswith(('.1','.2','.3')):
            out=open('%saugustus_proteins/%s.fas.%s' % (mainout,i[:-6],i[-1]),'w')
        count=0
        tr=0
        for line in f:
            if line.startswith('# start gene'):
                tr = 1
            elif tr == 1:
                line=line.split();places=[line[0],line[3],line[4]];tr=0
            elif line.startswith('# protein'):
                line=line.strip().split('[')
                count+=1
                out.write('>g%s[%s:%s-%s]\n' % (count,places[0],places[1],places[2]))
                if line[1][-1] == ']':
                    line[1] = line[1][:-1]
                out.write(line[1])
                check=1
            else:
                if line.startswith('# end'):
                    check = 0
                    out.write('\n')
                elif check == 1:
                    line = line.split()[1]
                if line[-1] == ']':
                    line = line[:-1]
                out.write(line)
        out.close()

# Run HMMer (genome mode)
if mode == 'genome' or mode == 'hmmer':
    print('*** Running HMMER to confirm orthology of predicted proteins ***')

    files = os.listdir(mainout+'augustus_proteins/')
    if os.path.exists(mainout+'hmmer_output') is False:
        os.makedirs(os.path.join(mainout,"hmmer_output"))
        #  os.system('mkdir %shmmer_output' % mainout)

    for i in files:
        if i.endswith('.fas'):
            f = open(mainout+'augustus_proteins/'+i)
            name = i[:-4]
            command = [
                "hmmsearch", "--domtblout",
                "{output_file}.out".format(output_file=os.path.join(mainout, 'hmmer_output', name)),
                "-Z",  "{db_size}".format(db_size=Z),
                "-o", "temp",
                "--cpu", str(cpus),
                "{group_file}.hmm".format(group_file=os.path.join(clade, 'hmms', name)),
                "{input_file}".format(input_file=os.path.join(
                    mainout, "augustus_proteins", i))
            ]
            command = " ".join(command)

        # os.system('hmmsearch --domtblout %(output_file)s.out -Z %(db_size)s -o temp --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' %
	     #  {'input_file':mainout+'/augustus_proteins/'+i,
        #    'db_size':Z,
        #    'cpu':cpus,
        #    'group_file':clade+'/hmms/'+name,
        #    'output_file':mainout+'hmmer_output/'+name})
        elif i.endswith(('.1','.2','.3')):
            f = open(mainout+'augustus_proteins/'+i)
            name = i[:-6]
            command = [
                "hmmsearch", "--domtblout",
                "{output_file}.out".format(
                    output_file=os.path.join(mainout, 'hmmer_output',
                                             "{0}.out.{1}".format(name, i[-1]))),
                "-Z",  "{db_size}".format(db_size=Z),
                "-o", "temp",
                "--cpu", str(cpus),
                "{group_file}.hmm".format(group_file=os.path.join(clade,'hmms', name)),
                "{input_file}".format(input_file=os.path.join(
                    mainout, "augustus_proteins", i
                ))
            ]
            command = " ".join(command)

            # os.system('hmmsearch --domtblout %(output_file)s -Z %(db_size)s -o temp --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' %
	         #    {'input_file':mainout+'/augustus_proteins/'+i,
            #      'db_size':Z,
            #      'cpu':cpus,
            #      'group_file':clade+'/hmms/'+name,
            #      'output_file':mainout+'hmmer_output/'+name+'.out.'+i[-1]})

# Run HMMer (transcriptome mode)
if mode == 'trans' or mode == 'transcriptome':
    print('*** Running HMMER to confirm transcript orthology ***')
    files = os.listdir(os.path.join(mainout, 'translated_proteins'))
    if os.path.exists(os.path.join(mainout, 'hmmer_output')) is False:
        os.makedirs(os.path.join(mainout, 'hmmer_output'))
    group = ''
    grouplist = []
    for i in files:
        if i.endswith('.fas'):
            f = open('%stranslated_proteins/%s' % (mainout,i))
            name = i[:-7]
            group = transdic[name]
            if group not in grouplist:
                grouplist.append(group)
                command = ["hmmsearch",
                           "--domtblout", "{0}".format(
                                    os.path.join(mainout, 'hmmer_output', group)),
                           "-Z", str(Z),
                           "--cpu", str(cpus),
                           "{group_file}.hmm".format(
                                group_file=os.path.join(
                                    clade, 'hmms', group)),
                           "{input_file}".format(
                               input_file=os.path.join(mainout, "translated_proteins", i))
                           ]
                command = " ".join(command)
                subprocess.call(command, shell=True)
            else:
                grouplist.append(group)
                command = ["hmmsearch", "--domtblout",
                           "{0}.out.{1}".format(
                               os.path.join(mainout, "hmmer_output", group),
                               grouplist.count(group)),
                           "-Z", str(Z),
                           "-o", "temp", "--cpu", str(cpus),
                           "{group_file}".format(
                               group_file=os.path.join(clade, 'hmms', group)),
                           "{input_file}".format(
                               input_file=os.path.join(mainout, "translated_proteins", i))
                           ]
                command = " ".join(str(_) for _ in command)
                subprocess.call(command, shell=True)

                # os.system('hmmsearch --domtblout %(output_file)s.out.%(count)s -Z %(db_size)s -o temp --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' %
	            # {'input_file': mainout+'/translated_proteins/'+i,
                 # 'db_size':Z,
                 # 'cpu':cpus,
                 # 'group_file':clade+'/hmms/'+group,
                 # 'output_file':mainout+'hmmer_output/'+group,
                 # 'count':str(grouplist.count(group))})

# OGS/Proteome module
if mode == 'OGS':
    if os.path.exists(os.path.join(mainout, 'hmmer_output')) is False:
        os.makedirs(os.path.join(mainout, 'hmmer_output'))
        # os.system('mkdir %shmmer_output' % mainout)
    files = os.listdir(os.path.join(clade, '/hmms'))
    f2 = open(os.path.join(clade, 'scores_cutoff'))  # open target scores file
    # Load dictionary of HMM expected scores and full list of groups
    score_dic = {}
    for i in f2:
        i = i.strip().split()
        try:
            score_dic[i[0]] = float(i[1]) 	# [1] = mean value; [2] = minimum value
        except:
            pass
    totalbuscos = len(list(score_dic.keys()))
    for i in files:
        name = i[:-4]
        if name in score_dic:
            command = ["hmmsearch", "--domtblout",
                       "{0}.out".format(
                           os.path.join(mainout, "hmmer_output", name)),
                       "-o", "temp", "-Z", str(Z),
                       "--cpu", str(cpus),
                       "{group_file}.hmm".format(
                           group_file=os.path.join(clade, "hmms", name)),
                       "{input_file}".format(
                           input_file=args['genome'])
                       ]
            command = " ".join(command)
            subprocess.call(command, shell=True)

            # os.system('hmmsearch --domtblout %(output_file)s.out -o temp  -Z %(db_size)s --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' %
	        # {'input_file':args['genome'],
             # 'db_size':Z,
             # 'cpu':cpus,
             # 'group_file':clade+'/hmms/'+name,
             # 'output_file':mainout+'hmmer_output/'+name})


# ##*******get list to be re-run
if mode == 'genome' or mode == 'hmmer':
    print('*** Parsing HMMER results ***')
    # Open the output file; if no name was specified the default name will be used
    f2 = open(os.path.join(clade, 'scores_cutoff'))	# open target scores file
    # Load dictionary of HMM expected scores and full list of groups
    score_dic = {}
    for i in f2:
        i = i.strip().split()
        try:
            score_dic[i[0]] = float(i[1])
        except:
            pass
    totalbuscos = len(list(score_dic.keys()))
    f = open('coordinates_%s' % args['abrev'])
    dic = {}
    for i in f:
        i = i.strip().split('\t')
        name = i[0]
        if name not in dic:
            dic[name] = [[i[1], i[2], i[3]]]  # scaffold,start and end
        elif name in dic:
            dic[name].append([i[1], i[2], i[3]])
# ##*********


# ###Make summary

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
if mode == 'genome' or mode != 'OGS' or mode == 'report' or mode == 'hmmer':
    temp=os.listdir('%s/hmmer_output' % mainout)
    files=[]
    for i in temp:
        if i.endswith(('.out','.1','.2','.3')):
            files.append(i)
    f=open('%s/lengths_cutoff' % clade)
    for line in f:
        line = line.strip().split()
        leng_dic[line[0]] = float(line[3])
        sd_dic[line[0]] = float(line[2])
    for entry in files:
        f = open('%s/hmmer_output/%s' % (mainout,entry))
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
                if tlen>30*qlen:
                    pass
                else:
                    if prot not in hit_dic.keys() and score >= score_dic[group]:
                        hit_dic[prot] = [[int(line[15]), int(line[16])]]
                    elif score >= score_dic[group]:
                        hit_dic[prot].append([int(line[15]),int(line[16])])
        length = measuring(hit_dic)
        try:		# get maximum length of the putative gene in question
            if len(length) == 1:
                length = length[0]
            else:
                length = max(length)+1
            sigma=abs(leng_dic[group]-length)/sd_dic[group]
            if sigma <= 2:
                complete.append(entry)
                cc.append(group)
            elif sigma > 2:
                frag.append(entry)
        except:
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

if mode == 'OGS':
    complete = {}
    frag = {}
    done = []
    fcc = []
    temp = os.listdir(os.path.join(mainout, 'hmmer_output'))
    files = []
    for i in temp:
        if i.endswith(('.out', '.1', '.2', '.3')):
            files.append(i)
    f = open(os.path.join(clade, 'lengths_cutoff'))
    for line in f:
        line = line.strip().split()
        leng_dic[line[0]] = float(line[3])
        sd_dic[line[0]] = float(line[2])
    for entry in files:
        f = open('%s/hmmer_output/%s' % (mainout,entry))
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
                prediction = line[0]
                if group not in complete:
                    complete[group]=[];frag[group]=[]
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
        except:
            pass
  # check the multi hits
    for entry in complete:
        if len(complete[entry]) == 0:
            pass
        elif len(complete[entry]) == 1: # complete
            cc.append(entry)
        elif len(complete[entry]) > 1:
            mcc.append(entry)
    for entry in frag:
        if len(complete[entry]) != 0:
            pass
        elif frag[entry] != []:
            fcc.append(entry)

# summarize results, print and write to output files
summary = open('short_summary_'+args['abrev'], 'w')
if mode == 'OGS':
    print('Total complete BUSCOs found in assembly (<2 sigma) :  {0}\t({1} duplicated).'.format(
          len(set(cc))+len(set(mcc)),
          len(mcc))
         )
    print('Total BUSCOs partially recovered (>2 sigma) :  {0}'.format((len(fcc))))
else:
    print('Total complete BUSCOs found in assembly (<2 sigma) :  {0}\t({1} duplicated).' .format(
          len(set(unique)),
          len(mcc))
    )
    print('Total BUSCOs partially recovered (>2 sigma) :  {0}'.format(fcc))
print('Total groups searched: {0}'.format(totalbuscos))
try:
    if mode != 'OGS':
        print('Total BUSCOs not found:  {0}'.format(
            totalbuscos-(len(set(cc))+fcc)))
    else:
        print('Total BUSCOs not found: {0}'.format(
            totalbuscos-(len(set(cc))+len(set(mcc))+len(fcc))))
except:
    print('Total BUSCOs not found:  {0}'.format(totalbuscos-(len(set(cc)) + len(fcc))))


summary.write('# Summarized BUSCO benchmarking for file: {0}\n'.format(args["genome"]))
summary.write('#BUSCO was run in mode: {0}\n\n'.format(mode))
if mode != 'OGS' and mode != 'trans':
    summary.write('Summarized benchmarks in BUSCO notation:\n')
    summary.write('\tC:{0}%[D:{1}%],F:{2}%,M:{3}%,n:{4}\n\n'.format(
        shrink((len(set(cc))+len(set(mcc)))/totalbuscos),
        shrink(len(set(mcc))/totalbuscos),
        shrink(fcc/totalbuscos),
        shrink((totalbuscos-(len(set(cc))+fcc))/totalbuscos),
        totalbuscos))
elif mode == 'OGS':
    summary.write('Summarized benchmarks in BUSCO notation:\n')
    summary.write('\tC:{0}%[D:{1}%],F:{2}%,M:{3}%,n:{4}\n\n'.format(
        shrink((len(set(cc))+len(set(mcc)))/totalbuscos),
        shrink(len(set(mcc))/totalbuscos),
        shrink(len(fcc)/totalbuscos),
        shrink((totalbuscos-(len(set(cc))+len(set(mcc))+len(fcc)))/totalbuscos),
        totalbuscos))
elif mode == 'trans':
    summary.write('Summarized benchmarks in BUSCO notation:\n')
    summary.write('\tC:{0}%[D:{1}%],F:{2}%,M:{3}%,n:{4}\n\n'.format(
        shrink(len(set(cc))/totalbuscos),
        shrink(len(set(mcc))/totalbuscos),
        shrink(fcc/totalbuscos),
        shrink((totalbuscos-(len(set(cc))+fcc))/totalbuscos),
        totalbuscos))

summary.write('Representing:\n')
if mode != 'trans' and mode != 'OGS':
    summary.write('\t{0}\tComplete Single-copy BUSCOs\n'.format(len(set(cc))))
    summary.write('\t{0}\tComplete Duplicated BUSCOs\n'.format(len(set(mcc))))
elif mode == 'OGS':
    summary.write('\t{0}\tComplete Single-copy BUSCOs\n'.format(
        len(set(cc)) + len(set(mcc))
    ))
    summary.write('\t{0}\tComplete Duplicated BUSCOs\n'.format(len(set(mcc))))
elif mode == 'trans':
    summary.write('\t{0}\tComplete Single-copy BUSCOs\n'.format(
        len(set(cc)) - len(set(mcc))))
    summary.write('\t{0}\tComplete Duplicated BUSCOs\n'.format(len(set(mcc))))
if mode != 'OGS':
    summary.write('\t{0}\tFragmented BUSCOs\n'.format(fcc))
    summary.write('\t{0}\tMissing BUSCOs\n'.format(
        totalbuscos - (len(set(cc)) + fcc)
    ))
elif mode == 'OGS':
    summary.write('\t{0}\tFragmented BUSCOs\n'.format(len(fcc)))
    summary.write('\t{0}\tMissing BUSCOs\n'.format(
        totalbuscos - (len(set(cc)) + len(set(mcc)) + len(fcc))
    ))

summary.write('\t{0}\tTotal BUSCO groups searched\n'.format(totalbuscos))
summary.close()
summary = open('full_table_%s' % args['abrev'],'w')
# write correct header
if mode == 'genome' or mode == 'report':
    summary.write('# BUSCO_group\tStatus\tScaffold\tStart\tEnd\tBitscore\tLength\n')
elif mode == 'OGS':
    summary.write('# BUSCO_group\tStatus\tGene\tBitscore\tLength\n')
elif mode == 'trans' or mode == 'transcriptome':
    summary.write('# BUSCO_group\tStatus\tTranscript\tBitscore\tLength\n')

temp = os.listdir(os.path.join(mainout, 'hmmer_output' ))
done = []
files = []
for i in temp:
    if i.endswith(('.out','.1','.2','.3')):
        files.append(i)

for i in files:
    if i.endswith('.out'):
        name = i[:-4]
        marker = 0
    elif i.endswith(('.1','.2','.3')):
        name = i[:-6]
        marker = int(i[-1])-1
    f = open(os.path.join(mainout, 'hmmer_output/{0}'.format(i)))
    score = []
    hit_dic = {}
    for line in f:
        if line.startswith('# '):
            pass
        else:
            line = line.strip().split()
            score.append(float(line[7]))
            group = line[3]
            prot = line[0]
            tlen = int(line[2])
            qlen = int(line[5])
            prediction = line[0]
            if prot not in hit_dic.keys() and float(line[7]) >= score_dic[group]:
                hit_dic[prot] = [
                    [int(line[15]), int(line[16]), line[7]]
                ]
            elif float(line[7])>=score_dic[group]:
                hit_dic[prot].append(
                    [int(line[15]), int(line[16]), line[7]]
                )
    length = measuring(hit_dic)
    if mode == 'genome' or mode == 'report' or mode == 'hmmer':
        if hit_dic == {}:
            pass
        elif i in complete and name not in mcc:
            summary.write('{0}\tComplete\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                name,
                dic[group][marker][0],
                dic[group][marker][1],
                dic[group][marker][2],
                max(score),
                max(length)+1))
        elif i in complete and name in mcc:
            summary.write('{0}\tDuplicated\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                name,
                dic[group][marker][0],
                dic[group][marker][1],
                dic[group][marker][2],
                max(score),
                max(length) + 1))
        elif i in frag and name not in cc and name not in done:
            summary.write('{0}\tFragmented\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                name,
                dic[group][marker][0],
                dic[group][marker][1],
                dic[group][marker][2],
                max(score),
                max(length)+1))
    elif mode == 'OGS':
        if hit_dic == {}:
            pass
        elif name in complete and name not in mcc and name in cc:
            summary.write('{0}\tComplete\t{1}\t{2}\t{3}\n'.format(
                name,
                complete[name][0][0],
                max(score),
                max(length)+1))
        elif name in mcc:
            for entry in complete[name]:
                summary.write('{0}\tDuplicated\t{1}\t{2}\t{3}\n'.format(
                    name,
                    entry[0],
                    entry[1],
                    entry[2]+1))
        elif name in fcc and name not in cc:
            summary.write('{0}\tFragmented\t{1}\t{2}\t{3}\n'.format(
                name,
                frag[name][0],
                max(score),
                max(length)+1))
    elif mode == 'trans' or mode == 'Transcriptome':
        if hit_dic == {}:
            pass
        elif i in complete and name not in mcc:
            summary.write('{0}\tComplete\t{1}\t{2}\t{3}\n'.format(
                name,
                dic[group][marker],
                max(score),
                max(length)+1))
        elif i in complete and name in mcc:
            summary.write('{0}\tDuplicated\t{1}\t{2}\t{3}\n'.format(
                name,
                dic[group][marker],
                max(score),
                max(length)+1))
        elif i in frag and name not in cc and name not in done:
            summary.write('{0}\tFragmented\t{1}\t{2}\t{3}\n'.format(
                name,
                dic[group][marker],
                max(score),
                max(length)+1))
    f.close()
summary.close()


f = open('full_table_%s' % args['abrev'],'r')
lista = []
for i in f:
    i = i.strip().split()
    if i[0] not in lista:
        lista.append(i[0])
f.close()

# get final list of missing buscos
out = open('missing_buscos_list_{0}'.format(args['abrev']),'w')
f = open('full_table_{0}'.format(args['abrev']),'a')
for i in score_dic.keys():
    if i in lista:
        pass
    else:
        out.write(i+'\n')
        f.write('%s\tMissing\n' % (i))
out.close()
f.close()

# ######retraining
exitFlag = 0

if mode == 'genome' or mode == 'genome' or mode == 'hmmer':
    if os.path.exists(os.path.join(mainout, 'selected')) is False:
        os.makedirs(os.path.join(mainout, 'selected'))
    if os.path.exists(os.path.join(mainout, 'gffs')) is False:
        os.makedirs(os.path.join(mainout, 'gffs'))
    if os.path.exists(os.path.join(mainout, 'gb')) is False:
        os.makedirs(os.path.join(mainout, 'gb'))

    f = open('full_table_%s'.format(args['abrev']))
    lista = []
    re_run = []
    for line in f:
        status = line.split()[1]
        if status == 'Complete':
            lista.append(line.split()[0])
        elif status == 'Missing' or status == 'Fragmented':
            re_run.append(line.split()[0])

        files = os.listdir('%s/hmmer_output' % mainout)
        chosen = []
        for i in files:
            if i.endswith('.out'):
                name=i[:-4]
                if name in lista:
                    os.system('cp %s/hmmer_output/%s %s/selected/' % (mainout,i,mainout))
                    chosen.append(i)

    for entry in chosen:
        f=open('%s/selected/%s' % (mainout,entry))
        out=open('%s/gffs/%s' % (mainout,entry),'w');choicy=''
        for line in f:
            if line.startswith('# '):
                pass
            elif choicy=='':
                choicy=line.split()[0].split('[')[0]
        f.close()
        f=open('%s/augustus/%s' % (mainout,entry));check=0
        for line in f:
            if line.startswith('# start gene'):
                name = line.strip().split()[-1]
                if name == choicy:
                    check = 1
                elif line.startswith('# '):
                    check = 0
                elif check == 1:
                    out.write(line)
    f.close()
    for entry in chosen:
        f=open(os.path.join(mainout, 'gffs', entry))
        subprocess.call(
            [os.path.join(os.environ["AUGUSTUS_CONFIG_PATH"], "..", "scripts", "gff2gbSmallDNA.pl"),
             f.name, args["genome"], "1000",
             os.path.join(mainout, "gb", "{0}.raw.gb".format(entry[:-4]))
             ], shell=False)
        # os.system('$AUGUSTUS_CONFIG_PATH/../scripts/gff2gbSmallDNA.pl %s/gffs/%s %s 1000 %s/gb/%s.raw.gb' %
        #           (mainout,entry,args['genome'],mainout,entry[:-4]))

        print('Training augustus gene predictor')
        # create new species config file from template
        subprocess.call(
            [os.path.join(os.environ["AUGUSTUS_CONFIG_PATH"],
                          "..", "scripts", "new_species.pl"),
             "--species={0}".format(args["abrev"])],
            shell=False)
        # Create training set by catting files
        with open("training_set_{0}".format(args['abrev']), 'w') as gb_out:
            for gb_name in iter(x for x in os.listdir(os.path.join(mainout, "gb")) if
                          x.endswith(".gb")):
                with open(gb_name, "rb") as gb_file:
                    shutil.copyfileobj(gb_file, gb_out)

        # Do training .. without optimisation
        subprocess.call(["etraining", "--species={0}".format(args["abrev"]),
                         "training_set_{0}".format(args['abrev'])], shell=False)

        # os.system('$AUGUSTUS_CONFIG_PATH/../scripts/new_species.pl --species=%s' % (args['abrev']))
        # os.system('cat %sgb/*.gb > training_set_%s' % (mainout,args['abrev']))
        # os.system('etraining --species=%s training_set_%s' % (args['abrev'],args['abrev'])) # train on new training set (complete single copy buscos)

    # train on new training set (complete single copy buscos)
    if args['long']:
        print('Optimizing augustus metaparameters, this may take around 20 hours')
        subprocess.call([
            os.path.join(os.environ["AUGUSTUS_CONFIG_PATH"], "..",
                         "scripts", "optimize_augustus.pl"),
            "--species={0}".format(args["abrev"]),
            "training_set_{0}".format(args['abrev'])], shell=False)
        # os.system('$AUGUSTUS_CONFIG_PATH/../scripts/optimize_augustus.pl --species=%s training_set_%s' % (args['abrev'],args['abrev']))
        # train on new training set (complete single copy buscos)
        subprocess.call([
            "etraining",
            "--species={0}".format(args["abrev"]),
            "training_set_{0}".format(args['abrev'])], shell=False)

        # os.system('etraining --species=%s training_set_%s' % (args['abrev'],args['abrev']))

    print('*** Re-running failed predictions with different constraints, total number {0} ***'.format(
          len(re_run)))
    done = []
    target_species = args['abrev']
    strings=[]
    hammers=[]
    seds = []
    ripped = []
    for item in re_run:
        if item not in dic: # no coordinates found
            pass
        elif len(dic[item])>1: # more than one target coordinate
            count=0
            for entry in dic[item]:
                count+=1;
                command='augustus --proteinprofile=%(clade)s/%(prot_profile)s.prfl --predictionStart=%(start_coord)s --predictionEnd=%(end_coord)s --species=%(species)s \"%(scaffold)s\" > %(output)s 2>/dev/null' % {'prot_profile':'prfl/'+item,'start_coord':entry[1],'end_coord':entry[2],'clade':clade,'species':target_species,'scaffold':entry[0]+args['abrev']+'_.temp','output':mainout+'augustus/'+item+'.out.'+str(count)}

                strings.append(command)

                command='sed -i \'1,3d\' %(group_name)s' % {'group_name':mainout+'augustus/'+item+'.out.'+str(count)};

                seds.append(command)

                ripped.append(item+'.out.'+str(count));

                command='hmmsearch --domtblout %(output_file)s -Z %(db_size)s -o temp --cpu %(cpu)s %(group_file)s.hmm %(input_file)s ' %   {'input_file':mainout+'/augustus_proteins/'+item+'.fas.'+str(count),'db_size':Z,'cpu':1,'group_file':clade+'/hmms/'+item,'output_file':mainout+'hmmer_output/'+item+'.out.'+str(count)}
                hammers.append(command)
        elif len(dic[item])==1:
            entry=dic[item][0]
            try:
                command='augustus --proteinprofile=%(clade)s/%(prot_profile)s.prfl --predictionStart=%(start_coord)s --predictionEnd=%(end_coord)s --species=%(species)s \"%(scaffold)s\" > %(output)s 2>/dev/null' %    {'prot_profile':'prfl/'+item,'start_coord':entry[1],'end_coord':entry[2],'clade':clade,'species':target_species,'scaffold':entry[0]+args['abrev']+'_.temp','output':mainout+'augustus/'+item+'.out'}
                strings.append(command)

                command='sed -i \'1,3d\' %(group_name)s' % {'group_name':mainout+'augustus/'+item+'.out'}
                seds.append(command)

                ripped.append(item);name=item

                command='hmmsearch --domtblout %(output_file)s.out -Z %(db_size)s -o temp --cpu %(cpu)s %(group_file)s.hmm %(input_file)s.fas' % {'input_file':mainout+'/augustus_proteins/'+name,'db_size':Z,'cpu':1,'group_file':clade+'/hmms/'+name,'output_file':mainout+'hmmer_output/'+name}
                hammers.append(command)
            except:
                pass
  # missing(mainout,args['abrev'],'missing_buscos_list_')

    print('Starting to run Augustus again....')
    queueLock=threading.Lock()
    workQueue=queue.Queue(len(strings))
    threads=[]
    threadID=1

    needed=len(strings)
    mark=0

  # Create new threads
    for tName in threadList:
        mark+=1
        thread = myThread(threadID, tName, workQueue)
        thread.start()
        threads.append(thread)
        threadID+=1
        if mark>=needed:
            break

  # Fill the queue
    queueLock.acquire()
    for word in strings:
        workQueue.put(word)
    queueLock.release()

  # Wait for queue to empty
    while not workQueue.empty():
      # print(workQueue2.get())
        pass
  # Notify threads it's time to exit
    exitFlag=1

# Wait for all threads to complete
    for t in threads:
        t.join()
    print("Exiting Main Thread")


    print('Starting to run SED....')
    exitFlag=0

    queueLock=threading.Lock()
    workQueue=queue.Queue(len(seds))
    threads=[]
    threadID=1

    needed=len(seds)
    mark=0

    # Create new threads
    for tName in threadList:
        mark+=1
        thread = myThread(threadID, tName, workQueue)
        thread.start()
        threads.append(thread)
        threadID+=1
        if mark>=needed:
            break

    queueLock.acquire()
    for word in seds:
        workQueue.put(word)
    queueLock.release()

    # Wait for queue to empty
    while not workQueue.empty():
        pass
    # Notify threads it's time to exit
    exitFlag = 1

    # Wait for all threads to complete
    for t in threads:
        t.join()
    print("Exiting Main Thread")

    print('Starting to run EXTRACT....')
    for entry in ripped:
        extract(mainout,entry)

    print('Starting to run HMMER....')
    exitFlag=0

    queueLock=threading.Lock()
    workQueue=queue.Queue(len(seds))
    threads=[]
    threadID=1

    needed=len(seds)
    mark=0

    # Create new threads
    for tName in threadList:
        mark+=1
        thread = myThread(threadID, tName, workQueue)
        thread.start()
        threads.append(thread)
        threadID+=1
        if mark>=needed:
            break

    queueLock.acquire()
    for word in hammers:
        workQueue.put(word)
    queueLock.release()

  # Wait for queue to empty
    while not workQueue.empty():
        pass
  # Notify threads it's time to exit
    exitFlag = 1

# Wait for all threads to complete
    for t in threads:
        t.join()
    print("Exiting Main Thread")


# ##retraining and running over
# clean up temporary files
if mode != 'OGS':
    for filename in iter(fileno for fileno in os.listdir(".") if
                         (fileno.endswith("{0}_.temp".format(args["abrev"])) or
                          fileno.endswith("{0}.nsq".format(args["abrev"])) or
                          fileno.endswith("{0}.nin".format(args["abrev"])) or
                          fileno.endswith("{0}.nhr".format(args["abrev"]))
                          )):
        os.remove(filename)

    shutil.move("{0}_tblastn".format(args["abrev"]),
                "run_{0}".format(args["abrev"]))
    shutil.move("short_summary_{0}".format(args["abrev"]),
                "run_{0}".format(args["abrev"]))

    # os.system('rm *%s_.temp' % args['abrev'])
    # os.system('rm %s.nsq %s.nin %s.nhr'  % (args['abrev'],args['abrev'],args['abrev']))
    # os.system('mv %s_tblastn run_%s' % (args['abrev'],args['abrev']))
    # os.system('mv short_summary_%s run_%s' % (args['abrev'],args['abrev']))
    if mode != 'trans':
        shutil.move("coordinates_{0}".format(args["abrev"]),
                    "run_{0}".format(args["abrev"]))
        # os.system('mv coordinates_%s run_%s' % (args['abrev'],args['abrev']))
    shutil.move("missing_buscos_list_{0}".format(args["abrev"]),
                "run_{0}".format(args["abrev"]))
    shutil.move("full_table_{0}".format(args["abrev"]),
                "run_{0}".format(args["abrev"]))

    # os.system('mv full_table_%s run_%s' % (args['abrev'],args['abrev']))
else:
    shutil.move("missing_buscos_list_{0}".format(args["abrev"]),
                "run_{0}".format(args["abrev"]))
    shutil.move("full_table_{0}".format(args["abrev"]),
                "run_{0}".format(args["abrev"]))
    shutil.move("short_summary_%s".format(args["abrev"]),
                "run_{0}".format(args["abrev"]))

    # os.system('mv missing_buscos_list_%s run_%s' % (args['abrev'],args['abrev']))
    # os.system('mv full_table_%s run_%s' % (args['abrev'],args['abrev']))
    # os.system('mv short_summary_%s run_%s' % (args['abrev'],args['abrev']))
# Report run time per step
print('Total running time:  ', time.time() - start_time, "seconds")

# parse results and write final summary
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
if mode in ("genome", "hmmer"):
# if mode=='genome' or mode=='genome' or mode=='hmmer':
    temp = os.listdir(os.path.join(mainout, 'hmmer_output'))
    files = []
    for i in temp:
        if i.endswith(('.out','.1','.2','.3')):
            files.append(i)
    f = open(os.path.join(clade, 'lengths_cutoff'))
    for line in f:
        line = line.strip().split()
        leng_dic[line[0]] = float(line[3])
        sd_dic[line[0]] = float(line[2])
    for entry in files:
        f = open('%s/hmmer_output/%s' % (mainout,entry))
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
                if tlen>30*qlen:
                    pass
                else:
                    if prot not in hit_dic.keys() and score>=score_dic[group]:
                        hit_dic[prot]=[[int(line[15]),int(line[16])]]
                    elif score>=score_dic[group]:
                        hit_dic[prot].append([int(line[15]),int(line[16])])
        length = measuring(hit_dic)
        try:		# get maximum length of the putative gene in question
            if len(length) == 1:
                length = length[0]
            else:
                length = max(length) + 1
            sigma = abs(leng_dic[group] - length) / sd_dic[group]
            if sigma <= 2:
                complete.append(entry)
                cc.append(group)
            elif sigma > 2:
                frag.append(entry)
        except:
            pass
  # check the multi hits
    for entry in complete:
        if entry.endswith('.out'):
            name=entry[:-4]
        else:
            name=entry[:-6]
        if name in done:
            if name not in mcc:
                mcc.append(name)
        done.append(name)
    for i in cc:
        if i not in mcc:
            unique.append(i)
    for entry in frag:
        if entry.endswith('.out'):
            name=entry[:-4]
        else:
            name=entry[:-6]
        if name not in done and entry not in complete:
            done.append(name)
            fcc+=1

    # summarize results, print and write to output files
    summary=open('short_summary_'+args['abrev'], 'w')
    print('Total complete BUSCOs found in assembly (<2 sigma) :  {0}\t({1} duplicated).'.format(
          len(set(unique)), len(mcc)))
    print('Total BUSCOs partially recovered (>2 sigma) :  {0}'.format(fcc))
    print('Total groups searched: {0}'.format(totalbuscos))
    try:
        print('Total BUSCOs not found:  {0}'.format(
            totalbuscos - (len(set(cc)) + fcc)))
    except:
        print('Total BUSCOs not found:  %s'.format(
            totalbuscos - (len(set(cc)) + len(fcc))))

    summary.write('# Summarized BUSCO benchmarking for file: %s\n#BUSCO was run in mode: %s\n\n' % (args['genome'],mode))
    summary.write('Summarized benchmarks in BUSCO notation:\n\tC:%s%%[D:%s%%],F:%s%%,M:%s%%,n:%s\n\n' % (shrink(len(set(cc))/totalbuscos),shrink(len(set(mcc))/totalbuscos),shrink(fcc/totalbuscos),shrink((totalbuscos-(len(set(cc))+fcc))/totalbuscos),totalbuscos))

    summary.write('Representing:\n')
    summary.write('\t%s\tComplete Single-Copy BUSCOs\n' % (len(set(cc))))
    summary.write('\t%s\tComplete Duplicated BUSCOs\n' % (len(set(mcc))))
    summary.write('\t%s\tFragmented BUSCOs\n' % (fcc))
    summary.write('\t%s\tMissing BUSCOs\n' % (totalbuscos-(len(set(cc))+fcc)))

    summary.write('\t%s\tTotal BUSCO groups searched\n' % (totalbuscos))
    summary.close()
    summary=open('full_table_'+args['abrev'],'w')
    summary.write('# BUSCO_group\tStatus\tScaffold\tStart\tEnd\tBitscore\tLength\n')

    temp=os.listdir('%shmmer_output' % mainout)
    done=[]
    files=[]
    for i in temp:
        if i.endswith(('.out','.1','.2','.3')):
            files.append(i)

    for i in files:
        if i.endswith('.out'):
            name=i[:-4];marker=0
        elif i.endswith(('.1','.2','.3')):
            name=i[:-6];marker=int(i[-1])-1
        f=open('%shmmer_output/%s'% (mainout,i));score=[]
        hit_dic={}
        for line in f:
            if line.startswith('# '):
                pass
            else:
                line=line.strip().split()
                score.append(float(line[7]))
                group = line[3]
                prot = line[0]
                tlen = int(line[2])
                qlen = int(line[5])
                prediction = line[0]
                if prot not in hit_dic.keys() and float(line[7])>=score_dic[group]:
                    hit_dic[prot]=[
                        [int(line[15]), int(line[16]), line[7]]]
                elif float(line[7])>=score_dic[group]:
                    hit_dic[prot].append(
                        [int(line[15]),int(line[16]),line[7]]
                    )
        length=measuring(hit_dic)
        if hit_dic=={}:
            pass
        elif i in complete and name not in mcc:
            summary.write('{0}\tComplete\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                name,
                dic[group][marker][0],
                dic[group][marker][1],
                dic[group][marker][2],
                max(score),
                max(length)+1))
        elif i in complete and name in mcc:
            summary.write('{0}\tDuplicated\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                name,
                dic[group][marker][0],
                dic[group][marker][1],
                dic[group][marker][2],
                max(score),
                max(length)+1))
        elif i in frag and name not in cc and name not in done:
            summary.write('{0}\tFragmented\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                name,
                dic[group][marker][0],
                dic[group][marker][1],
                dic[group][marker][2],
                max(score),
                max(length)+1))
    summary.close()


    f=open('full_table_%s' % args['abrev'],'r');lista=[]
    for i in f:
        i=i.strip().split()
        if i[0] not in lista:
            lista.append(i[0])
    out = open('missing_buscos_list_%s' % args['abrev'],'w')	# get final list of missing buscos
    f = open('full_table_%s' % args['abrev'],'a')
    for i in score_dic.keys():
        if i in lista:
            pass
        else:
            out.write('{0}\n'.format(i))
            f.write('{0}\tMissing\n'.format(i))
    out.close()
    f.close()

    shutil.move('short_summary_{0}'.format(args['abrev']),
                'run_{0}'.format(args['abrev']))
    shutil.move("missing_buscos_list_{0}".format(args['abrev']),
                'run_{0}'.format(args['abrev']))
    shutil.move("full_table_{0}".format(args['abrev']),
                'run_{0}'.format(args['abrev']))
    shutil.move("training_set_{0}".format(args['abrev']),
                'run_{0}'.format(args['abrev']))


    # os.system('mv short_summary_%s run_%s' % (args['abrev'],args['abrev']))
    # os.system('mv missing_buscos_list_%s run_%s' % (args['abrev'],args['abrev']))
    # os.system('mv full_table_%s run_%s' % (args['abrev'],args['abrev']))
    # os.system('mv training_set_%s run_%s' % (args['abrev'],args['abrev']))
