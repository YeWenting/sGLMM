__author__ = 'Haohan Wang'

import numpy as np

def generatePositions_mice():
    snpIDs = [line.strip() for line in open('../data/real/mice/snpID.txt')]
    text = [line.strip() for line in open('../data/real/mice/mapfile.txt')][1:]
    pos = [[] for i in range(20)]

    snpD = {}

    m = []

    for line in text:
        items = line.split()
        sid = items[0]
        c = items[1]
        p = items[2]

        if c != 'X':
            snpD[sid] = (int(c), int(p))
        else:
            snpD[sid] = (20, int(p))

    for sid in snpIDs:
        c, p = snpD[sid]
        pos[c-1].append(p)

    np.save('../data/real/mice/position', pos)

def generateCausals_mice():
    causal = [[] for i in range(20)]

    text = [line.strip() for line in open('../data/real/mice/immunity.csv')]
    genes = []
    for line in text[1:]:
        items = line.split(',')
        genes.append(items[0])

    text = [line.strip() for line in open('../data/real/mice/MRK_List1.rpt')]
    geneInfo = {}
    for line in text[1:]:
        items = line.split('\t')
        gid = items[0]
        c = items[1]
        s = items[3]
        e = items[4]
        if gid!='NULL':
            if c!='UN':
                try:
                    if c!='X':
                        geneInfo[gid] = (int(c), int(s), int(e))
                    else:
                        geneInfo[gid] = (20, int(s), int(e))
                except:
                    pass


    for gid in genes:
        if gid in geneInfo:
            c, s, e = geneInfo[gid]
            causal[c-1].append((s, e, gid))
        else:
            print gid

    # print causal

    np.save('../data/real/mice/causal_immunity', causal)

def getPositions_mice():
    return np.load('../data/real/mice/position.npy')

def getCausals_mice(i):
    if 49<=i<=57:
        return np.load('../data/real/mice/causal_glucose.npy')
    if  78<=i<=85:
        return np.load('../data/real/mice/causal_immunity.npy')
    if i>=86:
        return np.load('../data/real/mice/causal_insulin.npy')


if __name__ == '__main__':
    generatePositions_mice()
    generateCausals_mice()