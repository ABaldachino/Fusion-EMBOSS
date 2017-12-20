#!/usr/bin/python3
# coding: utf-8

nb_SQ = taille_SQ = GC = GC1 = GC2 = GC3 = pGC = pGC1 = pGC2 = pGC3 = 0

code = {'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
'TGC': 'C', 'TGT': 'C',
'GAC': 'D', 'GAT': 'D',
'GAA': 'E', 'GAG': 'E',
'TTC': 'F', 'TTT': 'F',
'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
'CAC': 'H', 'CAT': 'H',
'ATA': 'I', 'ATC': 'I', 'ATT': 'I',
'AAA': 'K', 'AAG': 'K',
'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'TTA': 'L', 'TTG': 'L',
'ATG': 'M',
'AAC': 'N', 'AAT': 'N',
'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
'CAA': 'Q', 'CAG': 'Q',
'AGA': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
'AGC': 'S', 'AGT': 'S', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
'TGG': 'W',
'TAC': 'Y', 'TAT': 'Y',
'TAA': '*', 'TAG': '*', 'TGA': '*'}

codon_nbr = {'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0,
'TGC': 0, 'TGT': 0,
'GAC': 0, 'GAT': 0,
'GAA': 0, 'GAG': 0,
'TTC': 0, 'TTT': 0,
'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0,
'CAC': 0, 'CAT': 0,
'ATA': 0, 'ATC': 0, 'ATT': 0,
'AAA': 0, 'AAG': 0,
'CTA': 0, 'CTC': 0, 'CTG': 0, 'CTT': 0, 'TTA': 0, 'TTG': 0,
'ATG': 0,
'AAC': 0, 'AAT': 0,
'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0,
'CAA': 0, 'CAG': 0,
'AGA': 0, 'AGG': 0, 'CGA': 0, 'CGC': 0, 'CGG': 0, 'CGT': 0,
'AGC': 0, 'AGT': 0, 'TCA': 0, 'TCC': 0, 'TCG': 0, 'TCT': 0,
'ACA': 0, 'ACC': 0, 'ACG': 0, 'ACT': 0,
'GTA': 0, 'GTC': 0, 'GTG': 0, 'GTT': 0,
'TGG': 0,
'TAC': 0, 'TAT': 0,
'TAA': 0, 'TAG': 0, 'TGA': 0}

codon_fraction = codon_frequence = dictembl = {}

#########################################################################################
#PARSERS
#########################################################################################

def parseEmboss(emboss_file):
    ''' Prend un fichier emboss et retourne un dictionnaire avec pour keys les noms des enzyme de restriction et pour
    valeurs une liste contenant le site de restriction, et les positions de coupe '''
    fichier = open(emboss_file, 'r')
    dictenz = {}

    ## Place dans un dictionnaire chaque enzyme (key) avec le site, la first cut 5', 3',
    # et la seconde cut 5', et 3'
    for line in fichier:
        if line[0] != "#":  # Ne prend pas en compte les commentaires
            coupe = line.split("\t")    # Decoupe le str en une liste
            coupe [8] = coupe[8].replace("\n", "") # retire le \n present
            dictenz[coupe[0]] = [coupe[1], coupe[5], coupe[6], coupe[7], coupe[8]]
    return dictenz

def parseEMBL(emblfile):
    ''' Lis un fichier embl et place le nom du ficher, comme clef,
     et la séquence nuléotidique, comme item, dans un dictionnaire global '''

    fichier = open(emblfile, "r")
    isseq = False
    seq = ''

    for line in fichier:
        ## Isseq == True, marque le debut de la sequence nucleotidique
        if re.search("SQ *Sequence", line) != None: 
            isseq = True

        ## Stock toute les lignes de nt dans un dict (key = SQ)
        elif isseq == True:
            seq += (re.sub(' *[0-9]*\n$', '', line)).replace(" ", "")
            # seq = seq.replace(" ", "")

    dictembl[emblfile] = seq


