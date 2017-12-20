#!/usr/bin/python3
# coding: utf-8

import re 
import argparse

parser = argparse.ArgumentParser(description="Restrict report position of restriction site on a DNA sequence, by printing it by default, and sorted by hit postion")
parser.add_argument("sequence", nargs="+", help="Nucleotide sequence in a embl format")
parser.add_argument("-e", "--enzyme", help="Use Restrict program, need the restriction enzymes in a file in the 37e emboss format")
parser.add_argument("-o", "--outputfile", help="The outputfile to use other than default (*_output.cusp/restrict)")
parser.add_argument("-a", "--alphabetic_sort", action="store_true", help="Sort the result by alphabetique order of the enzymes name")
parser.add_argument("-r", "--reverse_sort", action="store_true", help="Reverse the order of result")
parser.add_argument("-s", "--size", help="Minimum recognition site length")
args = parser.parse_args()


nb_SQ, taille_SQ, GC, GC1, GC2, GC3, pGC, pGC1, pGC2, pGC3 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
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

AA_nbr = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 
'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 
'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0, 
'*': 0}

codon_fraction, codon_frequence, dict_seq, dict_enz = {}, {}, {}, {}

#########################################################################################
#PARSERS
#########################################################################################

def parse_Emboss(emboss_file):
    ''' Exrait les sites de restrictions et les stocke dans un dictionnaire global dict_enz'''

    fichier = open(emboss_file, 'r')
    dict_enz = {}

    ## Place dans un dictionnaire chaque enzyme (key) avec le site, la first cut 5', 3',
    # et la seconde cut 5', et 3'
    for line in fichier:
        if line[0] != "#":  # Ne prend pas en compte les commentaires
            coupe = line.split("\t")    # Decoupe le str en une liste
            coupe [8] = coupe[8].replace("\n", "") # retire le \n present
            dict_enz[coupe[0]] = [coupe[1], coupe[5], coupe[6], coupe[7], coupe[8]]
    return dict_enz

def parse_EMBL(embl_file):
    ''' Exrait la séquence nucléotidique d'un fichier EMBL et la stocke dans un dictionnaire global dict_seq'''

    fichier = open(embl_file, "r")
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

    dict_seq[embl_file] = seq.upper()

def parse_FASTA(fasta_file):
    ''' Exrait les séquences nucléotidiques d'un fichier (multi)FASTA et les stocke dans un dictionnaire global dict_seq'''

    f = open(fasta_file,"r")
    isseq = False
    for line in f:
    #Pour chaque ligne chercher si elle commence par un ">"
        if line[0] == ">" and isseq == False:
            isseq = True
            id = re.search("^>[^|]|([^|])", line)
            name = id.group(1)
            seq = ''
            #Et prendre toutes les lignes suivantes commençant par une base ATCG.
        elif line[0] in 'ATCG' and isseq == True:
            seq = seq+line
            seq = seq.strip('\n')
            # Si on enchaîne sur une nouvelle séquence, lancer la recherche de séquence codante et le comptage puis réinitialiser SQ
        elif line[0] == ">" and isseq == True:
            dict_seq[name] = seq.upper()
            id = re.search("^>[^|]|([^|])", line)
            name = id.group(1)
            seq = ''
            # Si on tombe sur une ligne vide alors traiter la séquence et réinitialiser isseq
        elif line[0] == "" and isseq == True:
            dict_seq[name] = seq.upper()
            isseq = False
    #Si on arrive à la fin du fichier et que l'on a une séquence en attente, traiter la séquence.
    if isseq == True:
        isseq = False
    f.close()

##############################################################################################
#Traitement
##############################################################################################

def codante(seq):
    ''' Retourne une séquence codante à partir d'une séquence donnée '''

    #Recherche d'un ATG
    posATG = seq.find('ATG')
    seq = seq[posATG:]
    #Recherche d'un stop sur le même cadre de lecture
    stops = ['TAA','TAG','TGA']
    for ntp in range(0, len(SQ), 3):
        codon = seq[ntp: ntp + 3]
        if codon in stops:
            seq = seq[: ntp + 3]
            break
    return seq

def comptage(seq):
    ''' Compte différents paramètres dans la séquence donnée'''

    global nb_SQ, taille_SQ, GC, GC1, GC2, GC3
    nb_SQ = nb_SQ + 1
    taille_SQ = taille_SQ + len(SQ)
    for ntp in SQ:
        if ntp == 'G' or ntp == 'C':
            GC += 1
    for ntp in range(0, len(seq), 3):
        codon = seq[ntp: ntp + 3]
        if codon[0] == 'G' or codon[0] == 'C':
            GC1 += 1
        if codon[1] == 'G' or codon[1] == 'C':
            GC2 += 1
        if codon[2] == 'G' or codon[2] == 'C':
            GC3 += 1
        codon_nbr[codon] = codon_nbr[codon]+1
        AA_nbr[code[codon]] = AA_nbr[code[codon]]+1

def usage():
    ''' Calcule l'usage du code et les pourcentages de GC '''

    # Calcul des pourcentages en GC.
    global nb_SQ, taille_SQ, GC, GC1, GC2, GC3, pGC, pGC1, pGC2, pGC3
    pGC = (GC * 100) / taille_SQ
    pGC1 = (GC1 * 100) / (taille_SQ / 3)
    pGC2 = (GC2 * 100) / (taille_SQ / 3)
    pGC3 = (GC3 * 100) / (taille_SQ / 3)
    # Calcul de l'usage du code
    for codon in codon_nbr:
        if AA_nbr[code[codon]] != 0:
            codon_fraction[codon] = codon_nbr[codon] / AA_nbr[code[codon]]
        else:
            codon_fraction[codon] = 0
    rapport = 1000 / (taille_SQ / 3)
    for codon in codon_nbr:
        codon_frequence[codon] = codon_nbr[codon] * rapport
    for codon in codon_nbr:
        codon_nbr[codon] = int(codon_nbr[codon])

def impression_cusp():
    fo=open("output.cusp",'w')
    fo.write('#CdsCount:%i' %nb_SQ + '\n')
    fo.write('\n')
    fo.write("#Coding GC %.2f" % pGC+"%\n")
    fo.write("#1st letter GC %.2f" % pGC1+"%\n")
    fo.write("#2nd letter GC %.2f" % pGC2+"%\n")
    fo.write("#3rd letter GC %.2f" % pGC3+"%\n")
    fo.write('\n')
    fo.write('#Codon AA Fraction Frequency Number\n')
    for codon in code:
        fo.write("%s" %(codon)+"    ")
        fo.write("%s  " %(code[codon]))
        fo.write( "% 8.3f" %(codon_fraction[codon])+" ")
        fo.write( "% 9.3f" %(codon_frequence[codon])+" ")
        fo.write( "% 6i" %(codon_nbr[codon]))
        fo.write("\n")
    fo.close()

if __name__ == '__main__':

    if args.enzyme != None:  
        pass # faire restrict
    elif args.enzyme == None and args.alphabetic_sort != False or args.reverse_sort != False or args.size != None:
        print("Error, those options (alphabetic_sort,reverse_sort and size) are for restrict not cusp")
    else:
        for seq in dict_seq:
            seq_codante = codante(dict_seq[seq])
            comptage(seq_codante)
        usage()
        impression_cusp()


