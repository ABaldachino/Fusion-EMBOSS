#!/usr/bin/python3
# coding: utf-8

import re
import argparse
import operator


#########################################################################################
## PARSERS
#########################################################################################

def parse_Emboss(emboss_file):
    ''' 
    Read a emboss file and stock the information about restriction enzyme in a 
    global dictionary with name of enzyme as key and a list of characteristics 
    as values.

    :parm emboss_file: name of the file to read
    :type emboss_file: string
    :return: None
    '''

    fichier = open(emboss_file, 'r')
    ## Place dans un dictionnaire chaque enzyme (key) avec le site, la first cut 5', 3',
    # et la seconde cut 5', et 3'
    for line in fichier:
        if line[0] != "#":  # Ne prend pas en compte les commentaires
            coupe = line.split("\t")    # Decoupe le str en une liste
            coupe[8] = coupe[8].replace("\n", "") # retire le \n present
            dict_enz[coupe[0]] = [coupe[1], coupe[5], coupe[6], coupe[7], coupe[8]]

def parse_EMBL(embl_file):
    ''' 
    Read a embl file, extract the DNA sequence and stock it in a global dictionary 
    with it's id as key and the sequence as value.

    :parm embl_file: name of the file to read
    :type embl_file: string
    :return: None
    '''

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

def parse_FASTA(fasta):
    ''' 
    Read a fasta file, extract the DNA sequence and stock it in a global dictionary 
    with it's id as key and the sequence as value.

    :parm fasta: name of the file to read
    :type fasta: string
    :return: None
    '''

    file = open(fasta, "r")

    for line in file:
        if line[0] == ">":
            dict_seq[line[1:-1]] = ""
            seq_actuel = line[1:-1]
        elif line == "\n":
            pass
        else:
            dict_seq[seq_actuel] += line[:-1]

    file.close()

##############################################################################################
## RESTRICT
##############################################################################################

def reverse_dna(sequence):
    ''' 
    Return the cDNA of the given sequence

    :parm sequence: Sequence nucléotidique
    :type sequence: string
    :return: cDNA
    :rtype: string 
    '''
    
    sequence_reverse = ""
    for nt in sequence.upper():
        if nt == "A":
            sequence_reverse += "T"
        elif nt == "T":
            sequence_reverse += "A"
        elif nt == "G":
            sequence_reverse += "C"
        else:
            sequence_reverse += "G"
    return sequence_reverse


def restrict_search(dict_enz, sequence, size_restriction_site=4):
    ''' 
    Find the restriction site in the sequence

    :parm sequence: DNA sequence
    :type sequence: string
    :parm dict_enz: Dictionary of restriction enzyme with their names as keys and their characteristics as value in a list
    :type dict_enz: dict
    :parm size_restriction_site: minimum recognition size to be take in account
    :type size_restriction_site: int
    :return: list of hit
    :rtype: list
    '''

    # Initiation valeurs
    lresult =[]
    strand = "+"
    lseq = []

    # Lève un attribute error si size n'est pas >= à 1 
    if size_restriction_site < 1:
        print("size_restriction_site should be >= 1")
        raise AttributeError
    
    # Prend sequence et en fait une liste contenant le 5'3' et 3'5'
    lseq.append(sequence.upper()) # Attribue la sequence dans une liste
    lseq.append(reverse_dna(lseq[0])) # Ajoute le brin complementaire à la liste

    for seq in lseq:
        for key in dict_enz: # Prend les enzymes de restriction une à une
            enz = dict_enz[key]
            if len(enz[0]) >= int(size_restriction_site): # Verifie que la taille du site de restriction est suffisante
                exp_re = ""

                # Transformation du site de restriction à une expression regulière
                for ntp in enz[0].upper():
                    exp_re += IUPAC_dna[ntp]

                # Recherche
                all_match = re.finditer(exp_re, seq)

                # Formatage des resultats
                for match in all_match:
                    lresult.append([match.start(), match.end(), key, enz[0], strand])
        strand = "-"

    return lresult

def sortrestrict(lresult, reverse, tri_alphabetique):
    ''' 
    Sort result by position of hit or by alphabetical order

    :parm lresult: List of result of search restrict
    :type lresult: list
    :parm reverse: Sort the hit in reverse order 
    :type reverse: Booleen
    :parm tri_alphabetique: sort the result b alphabetical order
    :type tri_alphabetique: Booleen
    :return: list of hit
    :treturn: list
    '''

    if tri_alphabetique == True:
        tri_alphabetique = 2
    else:
        tri_alphabetique = 1

    result_sorted = sorted(lresult, key = operator.itemgetter(tri_alphabetique), reverse = reverse)
    return result_sorted

def count_hit(lresult):
    ''' 
    Count the number of hit by sequence

    :parm lresult: list of hit found by restrictsearch
    :type lresult: list
    :return: number of hit
    :rtype: int
    '''

    compteur = 0
    for hit in lresult:
        compteur += 1
    return compteur

def save_restrict(dict_result, outputfile):
    '''
    Save the result of restrict search in a file, can save result of 
    different sequence in one file

    :parm dict_result: dictionary of result found by restrictsearch in multiple (or not) sequence
    :type dic_result: dict
    :parm outputfile: Where save the result
    :type outputfile: string
    :return: None
    '''

    fichier = open(outputfile, 'w')
    for key in dict_result:
        lresult = dict_result[key]
        
        nb_hit = count_hit(lresult)
        fichier.write("\n" + "#"*60 + "\n")
        fichier.write("# Name sequence : " + str(key) +"\n")
        fichier.write("# nb_hit : " + str(nb_hit) + "\n")
        fichier.write("#"*60 + "\n" + "Start  End  Enzyme  Restriction site  Strand" + "\n")
        for hit in result_sorted:
            fichier.write(str(hit[0]) + ",     " + str(hit[1]) + ",  " + str(hit[2]) + ",        " + str(hit[3].upper()) + ",    " + hit[4] + "\n")


##############################################################################################
## CUSP
##############################################################################################

def codante(seq):
    ''' 
    Return an ORF from the DNA sequence given == 
    Retourne une séquence codante à partir d'une séquence donnée ?

    :parm seq: A DNA sequence
    :type seq: string
    '''

    #Recherche d'un ATG
    posATG = seq.find('ATG')
    seq = seq[posATG:]
    #Recherche d'un stop sur le même cadre de lecture
    stops = ['TAA','TAG','TGA']
    for ntp in range(0, len(seq), 3):
        codon = seq[ntp: ntp + 3]
        if codon in stops:
            seq = seq[: ntp + 3]
            break
    return seq

def comptage(seq):
    ''' Compte différents paramètres dans la séquence donnée'''

    global nb_SQ, taille_SQ, GC, GC1, GC2, GC3
    nb_SQ = nb_SQ + 1
    taille_SQ = taille_SQ + len(seq)
    for ntp in seq:
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

def impression_cusp(outputfile):
    ''' 
    Save the result of cusp in a table in the outpufile

    :parm outputfile: where save the file
    :type outputfile: string
    '''
    fo=open(outputfile,'w')
    fo.write('#CdsCount:%i' %nb_SQ + '\n')
    fo.write('\n')
    fo.write("#Coding GC %.2f" % pGC+"%\n")
    fo.write("#1st letter GC %.2f" % pGC1+"%\n")
    fo.write("#2nd letter GC %.2f" % pGC2+"%\n")
    fo.write("#3rd letter GC %.2f" % pGC3+"%\n")
    fo.write('\n')
    fo.write('#Codon AA Fraction Frequency Number\n')
    for codon in ordre_code:
        fo.write("%s" %(codon)+"    ")
        fo.write("%s  " %(code[codon]))
        fo.write( "% 8.3f" %(codon_fraction[codon])+" ")
        fo.write( "% 9.3f" %(codon_frequence[codon])+" ")
        fo.write( "% 6i" %(codon_nbr[codon]))
        fo.write("\n")
    fo.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Restrict report position of restriction site on a DNA sequence, by printing it by default, and sorted by hit postion")
    parser.add_argument("sequence", nargs="+", type=str, help="Nucleotide sequence in a embl format")
    parser.add_argument("-e", "--enzyme", type=str, help="Use Restrict program, need the restriction enzymes in a file in the 37e emboss format")
    parser.add_argument("-o", "--outputfile", help="The outputfile to use other than default (*_output.cusp/restrict)")
    parser.add_argument("-a", "--alphabetic_sort", action="store_true", help="Sort the result by alphabetique order of the enzymes name")
    parser.add_argument("-r", "--reverse_sort", action="store_true", help="Reverse the order of result")
    parser.add_argument("-s", "--size", default=4, type=int, help="Minimum recognition site length (default=4)")
    args = parser.parse_args()
    

    IUPAC_dna = {'A':'A', 'T':'T', 'G':'G', 'C':'C', 'R':'AG', 'M':'AC', 
                'W':'AT', 'Y':'CT', 'S':'CG', 'K':'GT', 'V':'ACG', 
                'H':'ACT', 'D':'AGT', 'B':'CGT', 'N':'AGCT'}

    ordre_code=['GCA', 'GCC', 'GCG', 'GCT', 'TGC', 'TGT', 'GAC', 'GAT', \
    'GAA', 'GAG', 'TTC', 'TTT', 'GGA', 'GGC', 'GGG', 'GGT', 'CAC', 'CAT', \
    'ATA', 'ATC', 'ATT', 'AAA', 'AAG', 'CTA', 'CTC', 'CTG', 'CTT', 'TTA', \
    'TTG', 'ATG', 'AAC', 'AAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CAA', 'CAG', \
    'AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT', 'AGC', 'AGT', 'TCA', 'TCC', \
    'TCG', 'TCT', 'ACA', 'ACC', 'ACG', 'ACT', 'GTA', 'GTC', 'GTG', 'GTT', \
    'TGG', 'TAC', 'TAT', 'TAA', 'TAG', 'TGA']

    code = {'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', \
    'TGC': 'C', 'TGT': 'C', \
    'GAC': 'D', 'GAT': 'D', \
    'GAA': 'E', 'GAG': 'E', \
    'TTC': 'F', 'TTT': 'F', \
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', \
    'CAC': 'H', 'CAT': 'H', \
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', \
    'AAA': 'K', 'AAG': 'K', \
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'TTA': 'L', 'TTG': 'L', \
    'ATG': 'M', \
    'AAC': 'N', 'AAT': 'N', \
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', \
    'CAA': 'Q', 'CAG': 'Q', \
    'AGA': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', \
    'AGC': 'S', 'AGT': 'S', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', \
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', \
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', \
    'TGG': 'W', \
    'TAC': 'Y', 'TAT': 'Y', \
    'TAA': '*', 'TAG': '*', 'TGA': '*'}

    codon_nbr = {'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0, \
    'TGC': 0, 'TGT': 0, \
    'GAC': 0, 'GAT': 0, \
    'GAA': 0, 'GAG': 0, \
    'TTC': 0, 'TTT': 0, \
    'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0, \
    'CAC': 0, 'CAT': 0, \
    'ATA': 0, 'ATC': 0, 'ATT': 0, \
    'AAA': 0, 'AAG': 0, \
    'CTA': 0, 'CTC': 0, 'CTG': 0, 'CTT': 0, 'TTA': 0, 'TTG': 0, \
    'ATG': 0, \
    'AAC': 0, 'AAT': 0, \
    'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0, \
    'CAA': 0, 'CAG': 0, \
    'AGA': 0, 'AGG': 0, 'CGA': 0, 'CGC': 0, 'CGG': 0, 'CGT': 0, \
    'AGC': 0, 'AGT': 0, 'TCA': 0, 'TCC': 0, 'TCG': 0, 'TCT': 0, \
    'ACA': 0, 'ACC': 0, 'ACG': 0, 'ACT': 0, \
    'GTA': 0, 'GTC': 0, 'GTG': 0, 'GTT': 0, \
    'TGG': 0, \
    'TAC': 0, 'TAT': 0, \
    'TAA': 0, 'TAG': 0, 'TGA': 0}


    AA_nbr = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, \
    'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, \
    'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, \
    'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0, \
    '*': 0}
    
    nb_SQ, taille_SQ, GC, GC1, GC2, GC3, pGC, pGC1, pGC2, pGC3 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    codon_fraction, codon_frequence, dict_seq, dict_enz, dict_result = {}, {}, {}, {}, {}



    ##Ouverture des séquences et stockage dans un dictionnaire global, determine si il s'agit d'un fichier
    for seq in args.sequence:
        if re.search("[fa|fasta]$", seq) != None:
            parse_FASTA(seq)
        elif re.search("embl$", seq) != None:
            parse_EMBL(seq)

    if args.enzyme != None:

        parse_Emboss(args.enzyme) # Creation dict_enz
        for key in dict_seq: # Prend les clés qui sont inconnu une à une
            print(key)
            seq = dict_seq[key] # S'en sert pour obtenir la sequence
            result = restrict_search(dict_enz, seq, args.size) # Cherche les sites de restriction
            result_sorted = sortrestrict(result, args.reverse_sort, args.alphabetic_sort) # Les tris
            dict_result[key] = result_sorted # Place le résultat dans un dictionnaire global

        ## Enregistrement des resultats
        if args.outputfile != None: # Là où a defini l'utilisateur
            save_restrict(dict_result, args.outputfile)
        else: # Par défaut dans *.restrict
            out = re.sub("[^.]*$", "", args.sequence[0]) + "restrict"
            save_restrict(dict_result, out)

    elif args.enzyme == None and args.alphabetic_sort != False or args.reverse_sort != False or args.size != 4:
        print("Error, those options (alphabetic_sort,reverse_sort and size) are for restrict not cusp")
    
    else:
        for seq in dict_seq:
            seq_codante = codante(dict_seq[seq])
            comptage(seq_codante)
        usage()

        ## Enregistrement des resultats
        if args.outputfile != None: # Là où a defini l'utilisateur
            impression_cusp(args.outputfile)
        else: # Par défaut dans *.restrict
            out = re.sub("[^.]*$", "", args.sequence[0]) + "cusp"
            impression_cusp(out)
