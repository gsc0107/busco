from collections import deque

__author__ = 'Luca Venturini'


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
