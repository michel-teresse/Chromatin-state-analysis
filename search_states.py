#!/usr/bin/python3

CHROMOSOME_COUNT = 5

import sys
from bisect import bisect

# Vérifie les arguments
if len(sys.argv) -1 != 2 :
    print( "syntax: ", sys.argv[0], "chrom_states_file data_file" )
    sys.exit(1)

chrom_states_file = sys.argv[1]
data_file = sys.argv[2]

# Les 3 tableaux ci-dessous seront indexés par le n° de chromosome (de 1 à n)
# Chaque éléments des tableaux sera une liste de position ou d'état pour le chromosome concerné
# L'indice zéro sera initialisé à None car non utilisé
a_start = []
a_end = []
a_state = []

# Initialise les 3 tableaux
for i in range(0, CHROMOSOME_COUNT +1) :
    a_start.append(i)
    a_start[i] = [] if i > 0 else None
    a_end.append(i)
    a_end[i] = [] if i > 0 else None
    a_state.append(i)
    a_state[i] = [] if i > 0 else None

# Lit le fichier des états et initialise les tableaux a_start, a_end et a_state
with open( chrom_states_file ) as f:  
    line_nb = 0
    for line in f :  
        line_nb += 1
        # saute la 1ère ligne (en-tête)
        if line_nb == 1 :
            continue
        (chrom, start, end, state) = line.split()  
        chrom = int(chrom)
        a_start[chrom].append( int(start) )
        a_end[chrom].append( int(end) )
        a_state[chrom].append( state )

# DEBUG: affiche les tableaux
# for chrom in range(1, CHROMOSOME_COUNT+1) :
#     for i in range(0, len(a_start[chrom])) :
#         print( chrom, a_start[chrom][i], a_end[chrom][i], a_state[chrom][i] )

# Lit le fichier de données
with open( data_file ) as f:  
    for line in f :  
        fields = line.split()  
        chrom     = int(fields[0])
        start     = int(fields[2])
        end       = int(fields[3])
        direction = fields[4]

        # Recherche l'état corespondant à la position de départ du gene courant
        i = bisect(a_start[chrom], start) -1
        state_list = a_state[chrom][i]

        # Recherche les états suivants dans l'intervalle
        for j in range(i +1, len(a_start[chrom])) :
            if a_start[chrom][j] > end :
                break
            # state est une chaîne de caractère !
            state_list = state_list + a_state[chrom][j]

        # Si la direction du gène est négative, il faut inverser la liste des états
        if direction == "-":
            state_list = state_list[::-1]

        # Affiche la ligne avec la liste des états en dernière colonne
        print( line.rstrip(), "\t", state_list )
