__author__ = 'Haohan Wang'

import numpy as np

def generateLabels_alz():
    threshold = 1
    snps = []
    text = [line.strip() for line in open('../data/real/alz/marker.txt')]
    for line in text[1:]:
        items = line.split()
        snps.append(items[1])

    text = [line.strip() for line in open('../data/real/alz/gwas_alzheimer.tsv')]
    snpD = {}
    for line in text[1:]:
        items = line.split('\t')
        try:
            pvalue = float(items[27])
        except:
            pvalue = 1
        snp = items[21]
        if pvalue <= threshold:
            snpD[snp] = 0

    label = np.zeros(len(snps))
    for s in snpD:
        if s in snps:
            label[snps.index(s)] = 1

    print sum(label)

    np.save('../data/real/alz/label', label)


def getLabels_alz():
    return np.load('../data/real/alz/label.npy')


if __name__ == '__main__':
    generateLabels_alz()