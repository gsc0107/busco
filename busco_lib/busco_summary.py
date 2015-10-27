import os
from busco_lib.busco_utils import shrink, measuring

__author__ = 'Luca Venturini'


def summarise(args, genome_dic, totalbuscos, mcc, cc, fcc, unique, score_dic, complete, frag):

    # summarize results, print and write to output files
    summary = open('short_summary_'+args['abrev'], 'w')
    if args.mode == 'OGS':
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
        if args.mode != 'OGS':
            print('Total BUSCOs not found:  {0}'.format(
                totalbuscos-(len(set(cc))+fcc)))
        else:
            print('Total BUSCOs not found: {0}'.format(
                totalbuscos-(len(set(cc))+len(set(mcc))+len(fcc))))
    except:
        print('Total BUSCOs not found:  {0}'.format(totalbuscos-(len(set(cc)) + len(fcc))))
    
    
    summary.write('# Summarized BUSCO benchmarking for file: {0}\n'.format(args["genome"]))
    summary.write('#BUSCO was run in mode: {0}\n\n'.format(args.mode))
    if args.mode != 'OGS' and args.mode != 'trans':
        summary.write('Summarized benchmarks in BUSCO notation:\n')
        summary.write('\tC:{0}%[D:{1}%],F:{2}%,M:{3}%,n:{4}\n\n'.format(
            shrink((len(set(cc))+len(set(mcc)))/totalbuscos),
            shrink(len(set(mcc))/totalbuscos),
            shrink(fcc/totalbuscos),
            shrink((totalbuscos-(len(set(cc))+fcc))/totalbuscos),
            totalbuscos))
    elif args.mode == 'OGS':
        summary.write('Summarized benchmarks in BUSCO notation:\n')
        summary.write('\tC:{0}%[D:{1}%],F:{2}%,M:{3}%,n:{4}\n\n'.format(
            shrink((len(set(cc))+len(set(mcc)))/totalbuscos),
            shrink(len(set(mcc))/totalbuscos),
            shrink(len(fcc)/totalbuscos),
            shrink((totalbuscos-(len(set(cc))+len(set(mcc))+len(fcc)))/totalbuscos),
            totalbuscos))
    elif args.mode == 'trans':
        summary.write('Summarized benchmarks in BUSCO notation:\n')
        summary.write('\tC:{0}%[D:{1}%],F:{2}%,M:{3}%,n:{4}\n\n'.format(
            shrink(len(set(cc))/totalbuscos),
            shrink(len(set(mcc))/totalbuscos),
            shrink(fcc/totalbuscos),
            shrink((totalbuscos-(len(set(cc))+fcc))/totalbuscos),
            totalbuscos))
    
    summary.write('Representing:\n')
    if args.mode != 'trans' and args.mode != 'OGS':
        summary.write('\t{0}\tComplete Single-copy BUSCOs\n'.format(len(set(cc))))
        summary.write('\t{0}\tComplete Duplicated BUSCOs\n'.format(len(set(mcc))))
    elif args.mode == 'OGS':
        summary.write('\t{0}\tComplete Single-copy BUSCOs\n'.format(
            len(set(cc)) + len(set(mcc))
        ))
        summary.write('\t{0}\tComplete Duplicated BUSCOs\n'.format(len(set(mcc))))
    elif args.mode == 'trans':
        summary.write('\t{0}\tComplete Single-copy BUSCOs\n'.format(
            len(set(cc)) - len(set(mcc))))
        summary.write('\t{0}\tComplete Duplicated BUSCOs\n'.format(len(set(mcc))))
    if args.mode != 'OGS':
        summary.write('\t{0}\tFragmented BUSCOs\n'.format(fcc))
        summary.write('\t{0}\tMissing BUSCOs\n'.format(
            totalbuscos - (len(set(cc)) + fcc)
        ))
    elif args.mode == 'OGS':
        summary.write('\t{0}\tFragmented BUSCOs\n'.format(len(fcc)))
        summary.write('\t{0}\tMissing BUSCOs\n'.format(
            totalbuscos - (len(set(cc)) + len(set(mcc)) + len(fcc))
        ))
    
    summary.write('\t{0}\tTotal BUSCO groups searched\n'.format(totalbuscos))
    summary.close()
    summary = open('full_table_%s' % args['abrev'],'w')
    # write correct header
    if args.mode == 'genome' or args.mode == 'report':
        summary.write('# BUSCO_group\tStatus\tScaffold\tStart\tEnd\tBitscore\tLength\n')
    elif args.mode == 'OGS':
        summary.write('# BUSCO_group\tStatus\tGene\tBitscore\tLength\n')
    elif args.mode == 'trans' or args.mode == 'transcriptome':
        summary.write('# BUSCO_group\tStatus\tTranscript\tBitscore\tLength\n')
    
    temp = os.listdir(os.path.join(args.mainout, 'hmmer_output' ))
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
        f = open(os.path.join(args.mainout, 'hmmer_output/{0}'.format(i)))
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
        if args.mode == 'genome' or args.mode == 'report' or args.mode == 'hmmer':
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
                    max(length) + 1))
            elif i in frag and name not in cc and name not in done:
                summary.write('{0}\tFragmented\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                    name,
                    genome_dic[group][marker][0],
                    genome_dic[group][marker][1],
                    genome_dic[group][marker][2],
                    max(score),
                    max(length)+1))
        elif args.mode == 'OGS':
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
        elif args.mode == 'trans' or args.mode == 'Transcriptome':
            if hit_dic == {}:
                pass
            elif i in complete and name not in mcc:
                summary.write('{0}\tComplete\t{1}\t{2}\t{3}\n'.format(
                    name,
                    genome_dic[group][marker],
                    max(score),
                    max(length)+1))
            elif i in complete and name in mcc:
                summary.write('{0}\tDuplicated\t{1}\t{2}\t{3}\n'.format(
                    name,
                    genome_dic[group][marker],
                    max(score),
                    max(length)+1))
            elif i in frag and name not in cc and name not in done:
                summary.write('{0}\tFragmented\t{1}\t{2}\t{3}\n'.format(
                    name,
                    genome_dic[group][marker],
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