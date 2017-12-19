import sys

#!ATTENTION!: Ce programme présuppose que les fichiers entrés sont en format txt embl ou fasta
#et que les séquences dans ces fichiers sont monocistroniques ou ne présentent pas d'introns.

#Définition des variables et des dictionnaires globaux:
#Nombre de séquences
nb_SQ = 0
#somme de la taille des séquences traduites
taille_SQ=0
#compteurs de GC dans la séquence traduite et à chaque position des codons
GC=0
GC1=0
GC2=0
GC3=0
#code génétique classique
code={}
code['GCA']='A'
code['GCC']='A'
code['GCG']='A'
code['GCT']='A'
code['TGC']='C'
code['TGT']='C'
code['GAC']='D'
code['GAT']='D'
code['GAA']='E'
code['GAG']='E'
code['TTC']='F'
code['TTT']='F'
code['GGA']='G'
code['GGC']='G'
code['GGG']='G'
code['GGT']='G'
code['CAC']='H'
code['CAT']='H'
code['ATA']='I'
code['ATC']='I'
code['ATT']='I'
code['AAA']='K'
code['AAG']='K'
code['CTA']='L'
code['CTC']='L'
code['CTG']='L'
code['CTT']='L'
code['TTA']='L'
code['TTG']='L'
code['ATG']='M'
code['AAC']='N'
code['AAT']='N'
code['CCA']='P'
code['CCC']='P'
code['CCG']='P'
code['CCT']='P'
code['CAA']='Q'
code['CAG']='Q'
code['AGA']='R'
code['AGG']='R'
code['CGA']='R'
code['CGC']='R'
code['CGG']='R'
code['CGT']='R'
code['AGC']='S'
code['AGT']='S'
code['TCA']='S'
code['TCC']='S'
code['TCG']='S'
code['TCT']='S'
code['ACA']='T'
code['ACC']='T'
code['ACG']='T'
code['ACT']='T'
code['GTA']='V'
code['GTC']='V'
code['GTG']='V'
code['GTT']='V'
code['TGG']='W'
code['TAC']='Y'
code['TAT']='Y'
code['TAA']='*'
code['TAG']='*'
code['TGA']='*'
#compteur des codons
codon_nbr={}
codon_nbr['GCA']=0
codon_nbr['GCC']=0
codon_nbr['GCG']=0
codon_nbr['GCT']=0
codon_nbr['TGC']=0
codon_nbr['TGT']=0
codon_nbr['GAC']=0
codon_nbr['GAT']=0
codon_nbr['GAA']=0
codon_nbr['GAG']=0
codon_nbr['TTC']=0
codon_nbr['TTT']=0
codon_nbr['GGA']=0
codon_nbr['GGC']=0
codon_nbr['GGG']=0
codon_nbr['GGT']=0
codon_nbr['CAC']=0
codon_nbr['CAT']=0
codon_nbr['ATA']=0
codon_nbr['ATC']=0
codon_nbr['ATT']=0
codon_nbr['AAA']=0
codon_nbr['AAG']=0
codon_nbr['CTA']=0
codon_nbr['CTC']=0
codon_nbr['CTG']=0
codon_nbr['CTT']=0
codon_nbr['TTA']=0
codon_nbr['TTG']=0
codon_nbr['ATG']=0
codon_nbr['AAC']=0
codon_nbr['AAT']=0
codon_nbr['CCA']=0
codon_nbr['CCC']=0
codon_nbr['CCG']=0
codon_nbr['CCT']=0
codon_nbr['CAA']=0
codon_nbr['CAG']=0
codon_nbr['AGA']=0
codon_nbr['AGG']=0
codon_nbr['CGA']=0
codon_nbr['CGC']=0
codon_nbr['CGG']=0
codon_nbr['CGT']=0
codon_nbr['AGC']=0
codon_nbr['AGT']=0
codon_nbr['TCA']=0
codon_nbr['TCC']=0
codon_nbr['TCG']=0
codon_nbr['TCT']=0
codon_nbr['ACA']=0
codon_nbr['ACC']=0
codon_nbr['ACG']=0
codon_nbr['ACT']=0
codon_nbr['GTA']=0
codon_nbr['GTC']=0
codon_nbr['GTG']=0
codon_nbr['GTT']=0
codon_nbr['TGG']=0
codon_nbr['TAC']=0
codon_nbr['TAT']=0
codon_nbr['TAA']=0
codon_nbr['TAG']=0
codon_nbr['TGA']=0
#compteurs des AA
AA_nbr={}
AA_nbr['A']=0
AA_nbr['C']=0
AA_nbr['D']=0
AA_nbr['E']=0
AA_nbr['F']=0
AA_nbr['G']=0
AA_nbr['H']=0
AA_nbr['I']=0
AA_nbr['K']=0
AA_nbr['L']=0
AA_nbr['M']=0
AA_nbr['N']=0
AA_nbr['P']=0
AA_nbr['Q']=0
AA_nbr['R']=0
AA_nbr['S']=0
AA_nbr['T']=0
AA_nbr['V']=0
AA_nbr['W']=0
AA_nbr['Y']=0
AA_nbr['*']=0
#Pourcentages de GC
pGC=0
pGC1=0
pGC2=0
pGC3=0
#usage du codon par AA
codon_fraction={}
#fréquence du codon dans une séquence de 1000AA
codon_fqc={}

#Définition des fonctions:

def traitementEmbl(fichier):
	f=open(fichier,"r")
	#utilisation d'un booléen pour marquer si l'on est dans une séquence ou non.
	bSQ=False
	#Recherche des différentes lignes portant les informations recherchées
	for line in f:
		if line[0:2] == 'SQ':
			bSQ=True
			SQ=""
		#Si la ligne commence par des espaces et que bSQ en True alors on récupère la séquence.
		elif line[0:2] == '  ' and bSQ==True:
			SQ=SQ+line[5:72]
			SQ=SQ.upper()
			SQ=SQ.replace(' ','')
		#Si aucune des conditions n'est satisfaite on laisse les booléens en False.
		else:
			if (bSQ==True):
				# Lancer ici la recherche d'une séquence codante et le comptage
				SQcod=codante(SQ)
				comptage(SQcod)
			bSQ=False
	f.close()

def traitementFASTA(fichier):
	f=open(fichier,"r")
	bSQ=False
	for line in f:
	#Pour chaque ligne chercher si elle commence par un ">"
		if line[0] == ">" and bSQ==False:
			bSQ=True
			SQ=''
			#Et prendre toutes les lignes suivantes commençant par une base ATCG.
		elif line[0] in 'ATCG' and bSQ == True:
			SQ=SQ+line
			SQ=SQ.strip('\n')
			# Si on enchaîne sur une nouvelle séquence, lancer la recherche de séquence codante et le comptage puis réinitialiser SQ
		elif line[0] == ">" and bSQ == True:
			SQcod=codante(SQ)
			comptage(SQcod)
			SQ=''
			# Si on tombe sur une ligne vide alors traiter la séquence et réinitialiser bSQ
		elif line[0] == "" and bSQ == True:
			SQcod=codante(SQ)
			comptage(SQcod)
			bSQ=False
	#Si on arrive à la fin du fichier et que l'on a une séquence en attente, traiter la séquence.
	if bSQ == True:
		SQcod=codante(SQ)
		comptage(SQcod)
		bSQ=False
	f.close()

def codante(SQ):
	#Recherche d'un ATG
	posATG=SQ.find('ATG')
	SQ=SQ[posATG:]
	#Recherche d'un stop sur le même cadre de lecture
	stops=['TAA','TAG','TGA']
	for ntp in range(0,len(SQ),3):
		codon=SQ[ntp:ntp+3]
		if codon in stops:
			SQ=SQ[:ntp+3]
			break
	return SQ

def comptage(SQ):
	global nb_SQ, taille_SQ, GC, GC1, GC2, GC3
	nb_SQ = nb_SQ + 1
	taille_SQ = taille_SQ + len(SQ)
	for ntp in SQ:
		if ntp=='G' or ntp=='C':
			GC += 1
	for ntp in range(0,len(SQ),3):
		codon=SQ[ntp:ntp+3]
		if codon[0]=='G' or codon[0]=='C':
			GC1 += 1
		if codon[1]=='G' or codon[1]=='C':
			GC2 += 1
		if codon[2]=='G' or codon[2]=='C':
			GC3 += 1
		codon_nbr[codon]=codon_nbr[codon]+1
		AA_nbr[code[codon]]=AA_nbr[code[codon]]+1

def usage():
	# Calcul des pourcentages en GC et de l'usage du code.
	global nb_SQ, taille_SQ, GC, GC1, GC2, GC3, pGC, pGC1, pGC2, pGC3
	pGC=(GC*100)/taille_SQ
	pGC1=(GC1*100)/(taille_SQ/3)
	pGC2=(GC2*100)/(taille_SQ/3)
	pGC3=(GC3*100)/(taille_SQ/3)
	for codon in codon_nbr:
		if AA_nbr[code[codon]] != 0:
			codon_fraction[codon]=codon_nbr[codon]/AA_nbr[code[codon]]
		else:
			codon_fraction[codon]=0
	rapport=1000/(taille_SQ/3)
	for codon in codon_nbr:
		codon_fqc[codon]=codon_nbr[codon]*rapport
	for codon in codon_nbr:
		type(codon_nbr[codon])
		codon_nbr[codon]=int(codon_nbr[codon])
		type(codon_nbr[codon])


if __name__ == '__main__':
	#Traitement des fichiers entrés
	input=sys.argv[:]
	for fichier in input:
		if ".txt" in fichier or ".embl" in fichier:
			traitementEmbl(fichier)
		elif ".fasta" in fichier:
			traitementFASTA(fichier)
	usage()
	#Ecriture des résultats
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
		fo.write( "% 9.3f" %(codon_fqc[codon])+" ")
		fo.write( "% 6i" %(codon_nbr[codon]))
		fo.write("\n")
	fo.close()










