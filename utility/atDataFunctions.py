__author__ = 'Haohan Wang'

import numpy as np

def generatePositions_at():
    text = [line.strip() for line in open('../data/real/at/genomeInformation.txt')]
    pos = [[] for i in range(5)]
    for line in text[1:]:
        items = line.split()
        pos[int(items[0])-1].append(int(items[1]))

    np.save('../data/real/at/position', pos)

def generateCausals_at():
    causal = []
    for i in range(5):causal.append([])
    text = [line.strip() for line in open('../data/real/at/ft.txt')]
    for line in text:
        items = line.split()
        causal[int(items[2])-1].append((int(items[0]), int(items[1])))

    np.save('../data/real/at/causal', causal)

def getCausals_Dict():
    text = [line.strip() for line in open('../data/real/at/phenotypeIndex.txt')]
    pheno2file = {}
    for line in text:
        items = line.split()
        fNAME = int(items[0][:-1])
        for item in items[1:]:
            pheno2file[int(item)] = fNAME
    return pheno2file

def getCausals_real(ind, pheno2file, phenoID):
    causal = []
    for i in range(5):causal.append([])
    fileName = pheno2file[phenoID]
    text = [line.strip() for line in open('../data/real/at/'+str(fileName)+'.txt')]
    for line in text:
        items = line.split()
        causal[int(items[2])-1].append((int(items[0]), int(items[1]), items[4]))
    np.save('../data/real/at/causal_'+str(ind), causal)

def getPositions_at():
    return np.load('../data/real/at/position.npy')

def getCausals_at(i):
    return np.load('../data/real/at/causal_'+str(i)+'.npy')

def matchPhenoPositions():
    text = [line.strip() for line in open('../data/real/at/phenotypeNames.txt')]
    names = [line.strip() for line in open('../data/real/at/phenotype_published_raw_organized.tsv')][0].split('\t')[2:]
    ind = []
    for n in names:
        n = n.strip()
        ind.append(text.index(n))
    print ind

def getI2P():
    text = [line.strip() for line in open('../data/real/at/phenotype_published_raw_organized.tsv')][0].split('\t')[2:]
    pind = []
    for name in text:
        e = name.find('_')
        p = int(name[:e])
        pind.append(p)
    return pind


if __name__ == '__main__':
    # generatePositions_at()
    # generateCausals_at()
    # matchPhenoPositions()
    p2f = getCausals_Dict()
    ind = [82, 83, 84, 65, 66, 67, 79, 80, 81, 88, 89, 90, 7, 100, 104, 105, 106, 91, 0, 1, 2, 3, 4, 5, 6, 39, 40,
           32, 33, 34, 73, 74, 75, 101, 102, 98, 99, 48, 49, 50, 51, 52, 71, 72]
    c = -1
    pind = getI2P()
    for i in ind:
        c += 1
        getCausals_real(i, p2f, pind[c])