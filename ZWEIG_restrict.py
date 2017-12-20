#############################
# Program : Restrict.py
# Author : Nathanael ZWEIG
#############################

import re
import sys
import operator
import argparse

parser = argparse.ArgumentParser(description="Restrict report position of restriction site on a DNA sequence, by printing it by default, and sorted by hit postion")
parser.add_argument("sequence", help="Nucleotide sequence in a embl format")
parser.add_argument("-e", "-enzyme", help="Use Restrict program, need the restriction enzymes in a file in the 37e emboss format")
parser.add_argument("-a", "-alphabetic_sort", help="Sort the result by alphabetique order of the enzymes name")
parser.add_argument("-r", "-reverse_sort", action="store_true", help="Reverse the order of result")
parser.add_argument("-s", help="Minimum recognition site length")
args = parser.parse_args()

restrict_help = "\
Restrict report position of restriction site on a DNA sequence, by printing it by default, and sorted by hit postion  \n\
    Standard qualifers:\n\
    [-enzymes]         The restriction enzymes in a file in the 37e emboss format\n\
    [-sequence]        Nucleotide sequence in a embl format\n\
    \n\
    Optionnal qualifers:\n\
    -a                  Sort the result by alphabetique order of the enzymes name \n\
    -r                  Reverse the order of result \n\
    -o                  Save the report in a file, will ask you for the destination \n\
                        (*.restrict.txt by default) \n\
    -t                  Test the result and print the hit with the restriction site for manual verification if necessary \n\
    -s                  Minimum recognition site length"


equ = {'A':'A', 'T':'T', 'G':'G', 'C':'C', 'R':'AG', 'M':'AC', 'W':'AT', 'Y':'CT', 'S':'CG', 'K':'GT', 'V':'ACG', \
       'H':'ACT', 'D':'AGT', 'B':'CGT', 'N':'AGCT'}

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
    ''' Lis un fichier embl et place l'identifiant (key = ID), l'organism (key = organism),
    la sequence nucleotidique (key = SQ) et la traduction (key = translation) dans un dictionnaire qu'elle retourne '''

    fichier = open(emblfile, "r")
    dictembl = {"SQ":""}
    stilltranslation = True
    isseq = False

    for line in fichier:

        ## Stock l'identifiant embl dans un dict (key = ID)
        if re.search("FT *source *", line) != None:
            line = re.sub("FT *source *", "", line)
            line = re.sub("\n$", "", line)
            dictembl["ID"] = line

        ## Stock l'organisme dans un dict (key = organism)
        elif re.search("/organism", line) != None:
            line = re.sub('FT */organism="', "", line)
            line = re.sub('"\n$', "", line)
            dictembl["organism"] = line

        ## Demarre le stockage de la traduction dans le dict (key = translation) 
        elif re.search("/translation=", line) != None:
            line = re.sub('FT */translation="', '', line)
            line = re.sub('\n$', '', line)
            dictembl["translation"] = line

        ## Continue le stockage de la traduction et y met fin
        elif "translation" in dictembl.keys()\
        and re.search("FT *", line) != None\
        and stilltranslation == True:
            line = re.sub('FT *', '', line)
            line = re.sub('\n$', '', line)
            if re.search('"$', line) != None:
                line = re.sub('"$', '', line)
                stilltranslation = False
            dictembl["translation"] += line

        ## Isseq == True, marque le debut de la sequence nucleotidique
        elif re.search("SQ *Sequence", line) != None: 
            isseq = True

        ## Stock toute les lignes de nt dans un dict (key = SQ)
        elif isseq == True:
            line = re.sub(' *[0-9]*\n$', '', line)
            line = line.replace(" ", "")
            if re.search("//$", line) != None:
                break
            dictembl["SQ"] += line

    return dictembl

def reverse_dna(sequence):
    ''' Retourne le brin complémentaire du brin d'adn donne '''
    
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


def restrictsearch(dictenz, dictembl, size_restriction_site):
    '''recherche les site de restriction dans la sequence donnee'''

    ## Initiation valeurs
    compteur_seq = 0
    compteur_nt = 0
    lrestult =[]
    hit = True
    strand = "+"
    lsequence = []
    
    lsequence.append(dictembl['SQ'].upper()) # Attribue la sequence dans une liste
    lsequence.append(reverse_dna(lsequence[0])) # Ajoute le brin complementaire à la liste

    for sequence in lsequence:
        for key in dictenz: # Prend les enzymes de restriction une à une
            enz = dictenz[key]
            
            if len(enz[0]) >= int(size_restriction_site): # Ne prend pas les enzymes de restriction avec un site de restriction inferieur ou egal a size_restriction
                while compteur_seq < len(sequence) - len(enz[0]): # Boucle pour scanner la sequence

                    for nt in sequence[compteur_seq:compteur_seq+len(enz[0])]: # Prend les nt un à un sur la longueur
                        # de la target de l'enzyme de restriction
                        if nt not in equ[enz[0][compteur_nt].upper()]: # Si le nt correspondant n'est pas dans la table des
                            hit = False # equivalance Hit devient faux

                        compteur_nt += 1 # Incremente le compteur de position de nt dans la portion etudie

                    ## Stockage du resultat si la seq correspond
                    if hit == True:
                        lrestult.append([compteur_seq, compteur_seq+len(enz[0]), key, enz[0], strand])

                    ## Reinitialisation valeurs et incrementation
                    else:
                        hit = True
                    compteur_nt = 0
                    compteur_seq += 1
                compteur_seq = 0
        strand = "-"

    return lrestult

def sortrestrict(lrestult, reverse, tri_alphabetique):
    ''' Tri les enzymes par position de hit croissante par défaut, peut-être triee par ordre alphabetique et en ordre decroissant '''

    result_sorted = sorted(lrestult, key = operator.itemgetter(tri_alphabetique), reverse = reverse)
    return result_sorted

def count_hit(lrestult):
    ''' Verifie le nombre de hit trouve par restrict.py '''
    
    compteur = 0
    for hit in lrestult:
        compteur += 1
    return compteur

def affichage_result(result):
    ''' Affiche 'proprement' le résultat '''
    
    print("Start ", "End ", "Enzyme", "Restriction", "Strand")
    for hit in result:
        print(hit[0], hit[1], hit[2], hit[3].upper(), hit[4])

def enregistrement_result(result_sorted, outputfile):
    ''' Enregistre le resultat dans un fichier '''

    fichier = open(outputfile, 'w')
    nb_hit = count_hit(result_sorted)

    fichier.write("#"*40 + "\n" + "# Commandline :" + "\n")
    for arg in sys.argv[1:]:
        fichier.write("# " + "      " + str(arg) + "\n")
    fichier.write("# nb_hit : " + str(nb_hit) + "\n")
    fichier.write("#"*40 + "\n" + "Start  End  Enzyme  Restriction site  Strand" + "\n")
    for hit in result_sorted:
        fichier.write(str(hit[0]) + ",  " + str(hit[1]) + ",    " + str(hit[2]) + ",    " + str(hit[3].upper()) + ",    " + hit[4] + "\n")


def main():

    ## Lecture des arguments donnes
    if "--help" in sys.argv or "-h" in sys.argv: # Ouverture de la documentation
        print(restrict_help)
        sys.exit(1)

    reverse = False
    tri_alphabetique = 1
    embossfile = sys.argv[1]
    emblfile = sys.argv[2]

    if "-r" in sys.argv: # Active la fonction reverse de sorted
        reverse = True

    if "-a" in sys.argv: # Active le tri des hit par ordre alphabetique
        tri_alphabetique = 2

    if "-o" in sys.argv and "-t" in sys.argv: # Retourne une erreurs si on utilise -t et -o en meme temps
        print("-t argument is not compatible with -o")
        sys.exit(1)

    if "-o" in sys.argv: # Obtention de l'outputfile
        outputfile = input("Outputfile (default name = *.restrict.txt) : ")
        if outputfile == '':
            outputfile = str(emblfile) + ".restrict.txt"
            print("Outputfile =", outputfile)

    if "-s" in sys.argv: # Obtention de la taille minimum d'un site de restriction pour restrictsearch
        size_restriction = input("Minimum recognition site length (default = 4) : ")
        if size_restriction == '':
            size_restriction = 4
    else:
        size_restriction = 4


    ## Execution des fonctions necessaires au programme
    emboss = parseEmboss(embossfile)
    dictembl = parseEMBL(emblfile)
    result = restrictsearch(emboss, dictembl, size_restriction)
    result_sorted = sortrestrict(result, reverse, tri_alphabetique)



    ## Fait appel au fonction d'affichage, d'enregistrement si l'argument '-o' est present
    if "-o" not in sys.argv and "-t" not in sys.argv:
        affichage_result(result_sorted)
    elif "-t" in sys.argv:
        test_restrict(result_sorted, dictembl)
    else:
        enregistrement_result(result_sorted, outputfile)

def test_restrict(result_sorted, dictembl):
    '''Test restrict /!\ la methode est très proche de celle dans restrictsearch(), il peut-être interessant de verifier certaine\
    valeurs a la main'''

    confirm = True
    compteur_nt = 0
    sequence_reverse = reverse_dna(dictembl["SQ"])

    for hit in result_sorted:

        if hit[4] == "+":
            for nt in dictembl["SQ"][hit[0]:hit[1]].upper():
                if nt not in equ[hit[3][compteur_nt].upper()].upper():
                    confirm = False
                compteur_nt += 1
            print(str(hit[3]), ", seq :", str(dictembl["SQ"][hit[0]:hit[1]].upper()), "confirm :", confirm)
            confirm = True
            compteur_nt = 0
        else:
            for nt in sequence_reverse[hit[0]:hit[1]].upper():
                if nt not in equ[hit[3][compteur_nt].upper()].upper():
                    confirm = False
                compteur_nt += 1
            print(str(hit[3]), ", seq :", str(sequence_reverse[hit[0]:hit[1]].upper()), "confirm :", confirm)
            confirm = True
            compteur_nt = 0

if __name__ == '__main__':

    print(args.enzyme)
    if args.enzyme != None:
        print("faire restrict")
    else:
        print("faire cusp")

    reverse = False
    tri_alphabetique = 1
    embossfile = args.enzymes
    emblfile = args.sequence