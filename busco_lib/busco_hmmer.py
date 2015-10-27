import os
import subprocess

__author__ = 'Luca Venturini'


def extract_proteins(args):

    print('*** Extracting predicted proteins ***')
    files = os.listdir(os.path.join(args.mainout, 'augustus'))
    check = 0
    for i in files:
        # WTH is this?
        subprocess.call('sed -i \'1,3d\' %saugustus/%s' % (args.mainout, i))
    if os.path.exists(os.path.join(args.mainout, 'augustus_proteins')) is False:
        os.makedirs(os.path.join(args.mainout, 'augustus_proteins'))
        # os.system('mkdir %saugustus_proteins' % mainout)

    for i in files:
        augustus_file = open(os.path.join(args.mainout, 'augustus', i))
        if i.endswith('.out'):
            out = open(os.path.join(args.mainout,
                                    "augustus",
                                    "{0}.fas".format(i[:-4])),
                       "w")
        elif i.endswith(('.1', '.2', '.3')):
            out = open(os.path.join(args.mainout,
                                    "augustus_proteins",
                                    "{0}.fas.{1}".format(i[:-6], i[-1])),
                       "w")
        count = 0
        tr = 0
        for line in augustus_file:
            if line.startswith('# start gene'):
                tr = 1
            elif tr == 1:
                line = line.split()
                places = [line[0], line[3], line[4]]
                tr = 0
            elif line.startswith('# protein'):
                line = line.strip().split('[')
                count += 1
                print(">g{count}[{chrom}:{start}-{end}]".format(count=count,
                                                                chrom=places[0],
                                                                start=places[1],
                                                                end=places[2]),
                      file=out)
                if line[1][-1] == ']':
                    line[1] = line[1][:-1]
                print(line[1].rstrip(), file=out)
                check = 1
            else:
                if line.startswith('# end'):
                    check = 0
                    print(file=out)
                elif check == 1:
                    line = line.split()[1]
                if line[-1] == ']':
                    line = line[:-1]
                print(line.rstrip(), file=out)
        out.close()


def do_hmmer_step4(args, transdic):

    # Run HMMer (genome mode)
    score_dic = {}
    if args.mode in ('genome', 'hmmer'):
        print('*** Running HMMER to confirm orthology of predicted proteins ***')
        extract_proteins(args)

        files = os.listdir(os.path.join(args.mainout, 'augustus_proteins/'))
        if os.path.exists(os.path.join(args.mainout, 'hmmer_output')) is False:
            os.makedirs(os.path.join(args.mainout, "hmmer_output"))

        for filename in files:
            if filename.endswith('.fas'):
                # f = open(os.path.join(args.mainout, 'augustus_proteins', filename))
                name = filename[:-4]
                command = [
                    "hmmsearch", "--domtblout",
                    "{output_file}.out".format(output_file=os.path.join(args.mainout,
                                                                        'hmmer_output',
                                                                        name)),
                    "-Z",  "{db_size}".format(db_size=args.Z),
                    "-o", "temp",
                    "--cpu", args.cpus,
                    "{group_file}.hmm".format(group_file=os.path.join(args.clade, 'hmms', name)),
                    "{input_file}".format(input_file=os.path.join(
                        args.mainout, "augustus_proteins", filename))
                ]
                command = " ".join(str(_) for _ in command)
                subprocess.call(command)

            elif filename.endswith(('.1', '.2', '.3')):
                # f = open(os.path.join(args.mainout,
                #                       'augustus_proteins',
                #                       i))
                name = filename[:-6]
                command = [
                    "hmmsearch", "--domtblout",
                    "{output_file}.out".format(
                        output_file=os.path.join(args.mainout, 'hmmer_output',
                                                 "{0}.out.{1}".format(name, filename[-1]))),
                    "-Z",  "{db_size}".format(db_size=args.Z),
                    "-o", "temp",
                    "--cpu", args.cpus,
                    "{group_file}.hmm".format(group_file=os.path.join(args.clade,
                                                                      'hmms', name)),
                    "{input_file}".format(input_file=os.path.join(
                        args.mainout, "augustus_proteins", filename
                    ))
                ]
                command = " ".join(str(_) for _ in command)
                subprocess.call(command)

    # Run HMMer (transcriptome mode)
    if args.mode in ('trans', 'transcriptome'):
        print('*** Running HMMER to confirm transcript orthology ***')
        files = os.listdir(os.path.join(args.mainout, 'translated_proteins'))
        if os.path.exists(os.path.join(args.mainout, 'hmmer_output')) is False:
            os.makedirs(os.path.join(args.mainout, 'hmmer_output'))
        grouplist = []
        for filename in files:
            if filename.endswith('.fas'):
                # f = open(os.path.join(args.mainout,
                #                       'translated_proteins',
                #                       filename))
                name = filename[:-7]
                group = transdic[name]
                if group not in grouplist:
                    grouplist.append(group)
                    command = ["hmmsearch",
                               "--domtblout", "{0}".format(
                                        os.path.join(args.mainout,
                                                     'hmmer_output',
                                                     group)),
                               "-Z", args.Z,
                               "--cpu", args.cpus,
                               "{group_file}.hmm".format(
                                    group_file=os.path.join(
                                        args.clade, 'hmms', group)),
                               "{input_file}".format(
                                   input_file=os.path.join(args.mainout,
                                                           "translated_proteins",
                                                           filename))
                               ]
                    command = " ".join(str(_) for _ in command)
                    subprocess.call(command, shell=True)
                else:
                    grouplist.append(group)
                    command = ["hmmsearch", "--domtblout",
                               "{0}.out.{1}".format(
                                   os.path.join(args.mainout,
                                                "hmmer_output",
                                                group),
                                   grouplist.count(group)),
                               "-Z", args.Z,
                               "-o", "temp", "--cpu", args.cpus,
                               "{group_file}".format(
                                   group_file=os.path.join(args.clade, 'hmms', group)),
                               "{input_file}".format(
                                   input_file=os.path.join(args.mainout,
                                                           "translated_proteins", filename))
                               ]
                    command = " ".join(str(_) for _ in command)
                    subprocess.call(command, shell=True)

    if args.mode == 'OGS':
        if os.path.exists(os.path.join(args.mainout, 'hmmer_output')) is False:
            os.makedirs(os.path.join(args.mainout, 'hmmer_output'))
            # os.system('mkdir %shmmer_output' % mainout)
        files = os.listdir(os.path.join(args.clade, '/hmms'))
        f2 = open(os.path.join(args.clade, 'scores_cutoff'))  # open target scores file
        # Load dictionary of HMM expected scores and full list of groups
        score_dic = {}
        for i in f2:
            i = i.strip().split()
            try:
                score_dic[i[0]] = float(i[1]) 	# [1] = mean value; [2] = minimum value
            except:
                pass
        totalbuscos = len(list(score_dic.keys()))
        for filename in files:
            name = filename[:-4]
            if name in score_dic:
                command = ["hmmsearch", "--domtblout",
                           "{0}.out".format(
                               os.path.join(args.mainout, "hmmer_output", name)),
                           "-o", "temp", "-Z", args.Z,
                           "--cpu", args.cpus,
                           "{group_file}.hmm".format(
                               group_file=os.path.join(args.clade, "hmms", name)),
                           "{input_file}".format(
                               input_file=args['genome'])
                           ]
                command = " ".join(str(_) for _ in command)
                subprocess.call(command, shell=True)


    # ##*******get list to be re-run
    if args.mode in ('genome', 'hmmer'):
        print('*** Parsing HMMER results ***')
        # Open the output file; if no name was specified the default name will be used
        f2 = open(os.path.join(args.clade, 'scores_cutoff'))	# open target scores file
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
    return transdic, score_dic
    # ##*********
