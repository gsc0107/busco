#!/usr/bin/env python3

"""
This module is focussed on retraining Augustus for our species of interest.
"""

from busco_lib.busco_utils import extract, measuring, shrink
import shutil
import os
import subprocess
import time
import multiprocessing

__author__ = 'Luca Venturini'


def launch_command(command):
    """
    Function to launch a command using subprocess.

    :param command: the command to launch
    :type command: str
    :return:
    """

    if not isinstance(command, str):
        if isinstance(command, list):
            command = " ".join([str(_) for _ in command])
        else:
            raise TypeError("Invalid command: {0}".format(command))

    launcher = subprocess.call(command, shell=True)
    if launcher > 0:
        raise subprocess.CalledProcessError(launcher, command)
    return launcher


def retraining(args, totalbuscos, genome_dic, score_dic):

    # TODO: improve this abysmal stub.
    """
    This method performs the real retraining.
    :param args:
    :param totalbuscos:
    :param genome_dic:
    :param score_dic:
    :return:
    """

    # ######retraining

    start_time = time.ctime()
    if args.mode == 'genome' or args.mode == 'genome' or args.mode == 'hmmer':
        if os.path.exists(os.path.join(args.mainout, 'selected')) is False:
            os.makedirs(os.path.join(args.mainout, 'selected'))
        if os.path.exists(os.path.join(args.mainout, 'gffs')) is False:
            os.makedirs(os.path.join(args.mainout, 'gffs'))
        if os.path.exists(os.path.join(args.mainout, 'gb')) is False:
            os.makedirs(os.path.join(args.mainout, 'gb'))
    
        f = open('full_table_%s'.format(args['abrev']))
        lista = []
        re_run = []
        chosen = []
        for line in f:
            status = line.split()[1]
            if status == 'Complete':
                lista.append(line.split()[0])
            elif status == 'Missing' or status == 'Fragmented':
                re_run.append(line.split()[0])
    
            files = os.listdir('%s/hmmer_output' % args.mainout)
            # chosen = []  # It does not make any sense to delete it at each line!
            for i in files:
                if i.endswith('.out'):
                    name = i[:-4]
                    if name in lista:
                        shutil.copy(os.path.join(args.mainout, "hmmer_output", i),
                                    os.path.join(args.mainout, "selected", i))
                        chosen.append(i)
    
        for entry in chosen:
            f = open('%s/selected/%s' % (args.mainout, entry))
            out = open('%s/gffs/%s' % (args.mainout, entry), 'w')
            choicy = ''
            for line in f:
                if line.startswith('# '):
                    pass
                elif choicy == '':
                    choicy = line.split()[0].split('[')[0]
            f.close()
            f = open('%s/augustus/%s' % (args.mainout, entry))
            check = 0
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
            f = open(os.path.join(args.mainout, 'gffs', entry))
            subprocess.call(
                [os.path.join(os.environ["AUGUSTUS_CONFIG_PATH"], "..", "scripts", "gff2gbSmallDNA.pl"),
                 f.name, args["genome"], "1000",
                 os.path.join(args.mainout, "gb", "{0}.raw.gb".format(entry[:-4]))
                 ], shell=False)
            # os.system('$AUGUSTUS_CONFIG_PATH/../scripts/gff2gbSmallDNA.pl %s/gffs/%s %s 1000 %s/gb/%s.raw.gb' %
            #           (args.mainout,entry,args['genome'],args.mainout,entry[:-4]))
    
            print('Training augustus gene predictor')
            # create new species config file from template
            subprocess.call(
                [os.path.join(os.environ["AUGUSTUS_CONFIG_PATH"],
                              "..", "scripts", "new_species.pl"),
                 "--species={0}".format(args["abrev"])],
                shell=False)
            # Create training set by catting files
            with open("training_set_{0}".format(args['abrev']), 'w') as gb_out:
                for gb_name in iter(_ for _ in os.listdir(os.path.join(args.mainout, "gb")) if _.endswith(".gb")):
                    with open(gb_name, "rb") as gb_file:
                        shutil.copyfileobj(gb_file, gb_out)
    
            # Do training .. without optimisation
            subprocess.call(["etraining", "--species={0}".format(args["abrev"]),
                             "training_set_{0}".format(args['abrev'])], shell=False)
    
        # train on new training set (complete single copy buscos)
        if args['long']:
            print('Optimizing augustus metaparameters, this may take around 20 hours')
            subprocess.call([
                os.path.join(os.environ["AUGUSTUS_CONFIG_PATH"], "..",
                             "scripts", "optimize_augustus.pl"),
                "--species={0}".format(args["abrev"]),
                "training_set_{0}".format(args['abrev'])], shell=False)
            subprocess.call([
                "etraining",
                "--species={0}".format(args["abrev"]),
                "training_set_{0}".format(args['abrev'])], shell=False)
    
        print('*** Re-running failed predictions with different constraints, total number {0} ***'.format(
              len(re_run)))
        target_species = args['abrev']
        strings = []
        hammers = []
        seds = []
        ripped = []
        for item in re_run:
            if item not in genome_dic:  # no coordinates found
                pass
            elif len(genome_dic[item]) > 1:  # more than one target coordinate
                count = 0
                for entry in genome_dic[item]:
                    count += 1
                    command = ['augustus',
                               '--proteinprofile={clade}/%(prot_profile)s.prfl'.format(
                                   os.path.join(args.clade,
                                                "{0}.prfl".format(os.path.join("prfl", item)))
                               ),
                               '--predictionStart={0}'.format(entry[1]),
                               '--predictionEnd={1}'.format(entry[2]),
                               '--species={0}'.format(target_species),
                               '\"{0}\"'.format(entry[0]+args.abrev+"_.temp"),
                               '>', os.path.join(args.mainout, "augustus",
                                                 "{0}.out.{1}".format(item, count)),
                               '2>/dev/null']
                    command = " ".join(str(_) for _ in command)
                    strings.append(command)

                    command = ["sed", "-i", "'1,3d'",
                               os.path.join(args.mainout, "augustus",
                                            "{0}.out.{1}".format(item, count))]

                    seds.append(" ".join(str(_) for _ in command))
    
                    ripped.append(item+'.out.'+str(count))

                    command = ["hmmsearch", "--domtblout",
                               os.path.join(args.mainout,
                                            'hmmer_output',
                                            "{0}.out.{1}".format(item, count)),
                               "-Z", args.Z,
                               "-o", "temp", "--cpu", 1,
                               "{hmm}".format(
                                   hmm=os.path.join(args.clade, "hmms", "{0}.hmm".format(item))),
                               "{input}".format(
                                   input=os.path.join(args.mainout, "augustus_proteins",
                                                      "{0}.fas.{1}".format(item, count)))
                               ]
                    command = " ".join(str(_) for _ in command)
                    hammers.append(command)
            elif len(genome_dic[item]) == 1:
                entry = genome_dic[item][0]
                try:
                    command = ["augustus",
                               "--proteinprofile={0}".format(
                                   os.path.join(args.clade,
                                                os.path.join("prlf", "{0}.prfl".format(item)))),
                               "--predictionStart={0}".format(entry[1]),
                               "--predictionEnd={0}", format(entry[2]),
                               "--species={0}".format(target_species),
                               "\"{0}\"".format(entry[0]+args['abrev']+'_.temp'),
                               ">",
                               "{output}".format(
                                   output=os.path.join(args.mainout, "augustus", "{0}.out".format(item))),
                               "2>/dev/null"]
                    command = " ".join(str(_) for _ in command)

                    strings.append(command)

                    command = ["sed", "-i", "'1,3d'",
                               os.path.join(args.mainout, "augustus",
                                            "{0}.out".format(item))]

                    seds.append(" ".join(str(_) for _ in command))
    
                    ripped.append(item)
                    name = item

                    command = ["hmmsearch", "--domtblout",
                               "{output}".format(
                                   output=os.path.join(
                                       args.mainout, "hmmer_output", "{0}.out".format(name))),
                               "-Z", args.Z,
                               "-o", "temp", "--cpu", 1,
                               os.path.join(args.clade, "hmms", "{0}.hmm".format(name)),
                               "{input_file}".format(
                                   input_file=os.path.join(args.mainout, "augustus_proteins",
                                                           "{0}.fas".format(name)))
                               ]
                    command = " ".join(str(_) for _ in command)
                    hammers.append(command)
                # TODO: Verify that these exceptions cover all necessary cases.
                except (OSError, IndexError):
                    pass
    
        print('Starting to run Augustus again....')

        pool = multiprocessing.Pool(processes=args.cpu)

        # Run augustus

        augustus_runner = pool.map_async(launch_command,
                                         strings)
        # Get all results before proceeding
        augustus_runner.get()

        sedder = pool.map_async(launch_command,
                                seds)

        # Get all SED results before proceeding
        sedder.get()
    
        print('Starting to run EXTRACT....')
        for entry in ripped:
            extract(args.mainout, entry)
    
        print('Starting to run HMMER....')
    
        # Get all HMMER results
        hammer_user = pool.map_async(launch_command, hammers)
        hammer_user.get()

        print("Exiting Main Thread")
    
    # ##retraining and running over
    # clean up temporary files
    if args.mode != 'OGS':
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
        if args.mode != 'trans':
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
    if args.mode in ("genome", "hmmer"):
        # if args.mode=='genome' or args.mode=='genome' or args.mode=='hmmer':
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
            f = open('%s/hmmer_output/%s' % (args.mainout, entry))
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
                    if tlen > 30 * qlen:
                        pass
                    else:
                        if prot not in hit_dic.keys() and score >= score_dic[group]:
                            hit_dic[prot] = [[int(line[15]), int(line[16])]]
                        elif score >= score_dic[group]:
                            hit_dic[prot].append([int(line[15]), int(line[16])])
            length = measuring(hit_dic)
            # get maximum length of the putative gene in question
            try:
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
            # TODO: verify that these exceptions cover the obvious cases
            except (TypeError, IndexError):
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
    
        # summarize results, print and write to output files
        summary = open('short_summary_'+args['abrev'], 'w')
        print('Total complete BUSCOs found in assembly (<2 sigma) :  {0}\t({1} duplicated).'.format(
              len(set(unique)), len(mcc)))
        print('Total BUSCOs partially recovered (>2 sigma) :  {0}'.format(fcc))
        print('Total groups searched: {0}'.format(totalbuscos))
        try:
            print('Total BUSCOs not found:  {0}'.format(
                totalbuscos - (len(set(cc)) + fcc)))
        # TODO: Verify that these exceptions cover all possible cases
        except (TypeError, IndexError):
            print('Total BUSCOs not found:  %s'.format(
                totalbuscos - (len(set(cc)) + len(fcc))))
    
        summary.write('# Summarized BUSCO benchmarking for file: %s\n' % args['genome'])
        summary.write('#BUSCO was run in args.mode: %s\n\n' % args.mode)
        summary.write('Summarized benchmarks in BUSCO notation:\n')
        summary.write('\tC:%s%%[D:%s%%],F:%s%%,M:%s%%,n:%s\n\n' % (shrink(len(set(cc))/totalbuscos),
                                                                   shrink(len(set(mcc))/totalbuscos),
                                                                   shrink(fcc/totalbuscos),
                                                                   shrink((totalbuscos-(len(set(cc))+fcc))/totalbuscos),
                                                                   totalbuscos))
    
        summary.write('Representing:\n')
        summary.write('\t%s\tComplete Single-Copy BUSCOs\n' % len(set(cc)))
        summary.write('\t%s\tComplete Duplicated BUSCOs\n' % len(set(mcc)))
        summary.write('\t%s\tFragmented BUSCOs\n' % fcc)
        summary.write('\t%s\tMissing BUSCOs\n' % (totalbuscos - (len(set(cc))+fcc)))
    
        summary.write('\t%s\tTotal BUSCO groups searched\n' % totalbuscos)
        summary.close()
        summary = open('full_table_' + args['abrev'], 'w')
        summary.write('# BUSCO_group\tStatus\tScaffold\tStart\tEnd\tBitscore\tLength\n')
    
        temp = os.listdir('%shmmer_output' % args.mainout)
        done = []
        files = []
        for i in temp:
            if i.endswith(('.out', '.1', '.2', '.3')):
                files.append(i)
    
        for i in files:
            if i.endswith('.out'):
                name = i[:-4]
                marker = 0
            elif i.endswith(('.1', '.2', '.3')):
                name = i[:-6]
                marker = int(i[-1]) - 1
            f = open('%shmmer_output/%s' % (args.mainout, i))
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
                    # tlen = int(line[2])
                    # qlen = int(line[5])
                    # prediction = line[0]
                    if prot not in hit_dic.keys() and float(line[7]) >= score_dic[group]:
                        hit_dic[prot] = [
                            [int(line[15]), int(line[16]), line[7]]]
                    elif float(line[7]) >= score_dic[group]:
                        hit_dic[prot].append(
                            [int(line[15]), int(line[16]), line[7]])
            length = measuring(hit_dic)
            if hit_dic == {}:
                pass
            elif i in complete and name not in mcc:
                summary.write('{0}\tComplete\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                    name,
                    genome_dic[group][marker][0],
                    genome_dic[group][marker][1],
                    genome_dic[group][marker][2],
                    max(score),
                    max(length)+1))
            elif i in complete and name in mcc:
                summary.write('{0}\tDuplicated\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                    name,
                    genome_dic[group][marker][0],
                    genome_dic[group][marker][1],
                    genome_dic[group][marker][2],
                    max(score),
                    max(length)+1))
            elif i in frag and name not in cc and name not in done:
                summary.write('{0}\tFragmented\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                    name,
                    genome_dic[group][marker][0],
                    genome_dic[group][marker][1],
                    genome_dic[group][marker][2],
                    max(score),
                    max(length)+1))
        summary.close()
    
        f = open('full_table_{0}'.format(args.abrev), 'r')
        lista = []
        for i in f:
            i = i.strip().split()
            if i[0] not in lista:
                lista.append(i[0])
        out = open('missing_buscos_list_%s' % args['abrev'], 'w')  # get final list of missing buscos
        f = open('full_table_%s' % args['abrev'], 'a')
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
