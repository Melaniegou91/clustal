from Bio import SeqIO
import pandas as pd
import numpy as np


with open("sequences.fasta", "r") as fasta_file:  # Ouverture du fichier
    # Creation d'un dictionnaire
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    blosum = pd.read_csv("blossum62.txt")  # Lecture de la matrice blosum 62
    id_sequence = list(record_dict.keys())  # Contient les clés du dictionnaire
    def needlman(key_sequence1, key_sequence2, gap, blosum):
        sequence1 = record_dict[key_sequence1].seq
        sequence2 = record_dict[key_sequence2].seq
        # sequence1 = key_sequence1
        # sequence2 = key_sequence2
        value = 0
        matrice = np.zeros((len(sequence1) + 1, len(sequence2) + 1))
        # matrice[0][0] = 0
        for i in range(1, len(sequence1) + 1):
            matrice[i][0] = value + gap  # Affecte les penalites (lignes)
            value = matrice[i][0]

        value = 0

        for j in range(1, len(sequence2) + 1):
            matrice[0][j] = value + gap  # Affecte les penalites (lignes)
            value = matrice[0][j]


        for i in range(1, len(sequence1) + 1):
            for j in range(1, len(sequence2) + 1):
                # print(f"seq1:{sequence1[i]}, seq2:{sequence2[j]}")
                # print(blosum[sequence1[i]][sequence2[j]])
                diagonale = matrice[i - 1][j - 1] + blosum[sequence1[i-1]][sequence2[j-1]]
                gap_haut = matrice[i - 1][j] + gap
                gap_gauche = matrice[i][j - 1] + gap
                matrice[i][j] = max(diagonale, gap_haut, gap_gauche)        
        return matrice[i][j],matrice
    
    def alignement_seq(seq1, seq2,gap,score,blosum):
        seq1 = record_dict[seq1].seq
        seq2 = record_dict[seq2].seq
        i = len(seq1)
        j = len(seq2)
        alignement1 = ""
        alignement2 = ""
        while i > 0 and j > 0:
            if (score[i][j] == score[i - 1][j - 1] + blosum[seq1[i-1]][seq2[j-1]]):
                i = i - 1
                j = j - 1
                alignement1+=seq1[i - 1]
                alignement2+=seq2[j - 1]

            elif (score[i][j] == score[i][j - 1] + gap):
                j = j - 1
                alignement1+="-"
                alignement2+=seq2[j - 1]

            else:
                i = i - 1
                alignement1+=seq1[i - 1]
                alignement2+="-"
        while i > 0:
            i-=1
        while j > 0:
            j-= 1
        return alignement1, alignement2
    
    def matrice_dist(id_sequence, gap, blosum):
        matrice_dist = np.zeros((len(id_sequence), len(id_sequence)))
        for i in range(1, len(id_sequence)):
            for j in range(0, i ):
                print(i,j)
                
                matrice_dist[i][j] = needlman(id_sequence[i], 
                                              id_sequence[j], 
                                              gap = gap, 
                                              blosum = blosum)[0]
                print(matrice_dist[i][j])
                print(matrice_dist)
        print(matrice_dist)
            


#score = needlman(key_sequence1='MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH', key_sequence2='MVHLTAHHFGLWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLNALAHKYH', gap = - 5, blosum=blosum)
score = needlman(id_sequence[1], id_sequence[0], gap = -5, blosum=blosum)
print(score[0])

#alignement = alignement_seq(id_sequence[0], id_sequence[1], gap = -5, score = score[1], blosum=blosum)
#print("L'alignement optimal est : " + "\n","".join(alignement[0]), "".join(alignement[1]))

matrice_dist(id_sequence=id_sequence,gap = -5, blosum=blosum)
 