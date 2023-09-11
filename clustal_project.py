from Bio import SeqIO
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage,dendrogram


with open("sequences.fasta", "r") as fasta_file:  # Ouverture du fichier
    # Dictionary containing the sequences
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    # reading blosum 62 matrix 
    blosum = pd.read_csv("blossum62.txt") 

    # List containing sequence keys 
    id_sequence = list(record_dict.keys())  
    def needlman(key_sequence1, key_sequence2, gap, blosum):

        '''
        Needleman algorithm to provide
        the similarity score for each 
        pair of sequences

        Parameters
        ----------
        key_sequence1 : float
            The keys to the first sequence to study .
        key_sequence2 : float
            The keys to the second sequence to study.
        gap : int
            The penality of a gap
        blosum : matrix
            the blosum 62 matrix


    Returns
    -------
        Returns the alignment score from sequences
        '''

        # Retrive the sequences
        sequence1 = record_dict[key_sequence1].seq
        sequence2 = record_dict[key_sequence2].seq
        value = 0

        # Produce a matrix of size (n + 1) * (m + 1)
        # Where n and m correspond to the lengths of 
        # the sequences 
        matrice = np.zeros((len(sequence1) + 1, len(sequence2) + 1))
        # Add to all the lines in the first column the penalties
        for i in range(1, len(sequence1) + 1):
            matrice[i][0] = value + gap  # Affecte les penalites (lignes)
            value = matrice[i][0]
        value = 0
        # Add to all the columns in the first line the penalties
        for j in range(1, len(sequence2) + 1):
            matrice[0][j] = value + gap  # Affecte les penalites (lignes)
            value = matrice[0][j]
        # Allows to calculate the score according 
        # to the matrix blosum 62 and the penalty due to the gap    
        for i in range(1, len(sequence1) + 1):
            for j in range(1, len(sequence2) + 1):
                diagonale = matrice[i - 1][j - 1] + blosum[sequence1[i-1]][sequence2[j-1]]
                gap_haut = matrice[i - 1][j] + gap
                gap_gauche = matrice[i][j - 1] + gap
                # Allows you to select the maximum between the 3 possibilities
                matrice[i][j] = max(diagonale, gap_haut, gap_gauche)        
        return matrice[i][j],matrice
    
    def alignement_seq(key_sequence1, key_sequence2,gap,score,blosum):
        '''
        Algorithm for pairwise alignment

        Parameters
        ----------
        key_sequence1 : float
            The keys to the first sequence to study .
        key_sequence2 : float
            The keys to the second sequence to study.
        gap : int
            The penality of a gap
        blosum : matrix
            the blosum 62 matrix

    Returns
    -------
        returns the alignment in pairs
        '''

        # Retrive the sequences
        seq1 = record_dict[key_sequence1].seq
        seq2 = record_dict[key_sequence2].seq
        # i and j contain the length of the sequences
        i = len(seq1)
        j = len(seq2)
        # Empty sequence initialization
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
        
    def matrice_score(id_sequence, gap, blosum):
        matrice_dist = np.zeros((len(id_sequence), len(id_sequence)))
        for i in range(1, len(id_sequence)):
            for j in range(0, i):
                matrice_dist[i][j] = needlman(id_sequence[i], 
                                              id_sequence[j], 
                                              gap = gap, 
                                              blosum = blosum)[0]
        print(matrice_dist)
        return matrice_dist
    
    def value_max(matrice):        
        value_max= matrice[1][0]
        max_i = 1
        max_j = 0
        for i in range(1, len(id_sequence)):
            for j in range(0, i):
                if (matrice[i][j] > value_max):
                    value_max = matrice[i][j]
                    max_i = i
                    max_j = j
        return max_i, max_j, value_max
    
    def matrice_dist(matrice_score):
        '''

        '''
        print(matrice_score)
        maximum = value_max(matrice)[2]
        for i in range(1, len(id_sequence)):
            for j in range(0, i):
                matrice_score[i][j] = 1/ matrice_score[i][j]
        print(matrice_score)
        return matrice_score

    def embranchement_sequentiel(dist_matrix) : 
        dist_matrix = matrice_dist(matrice)
        print(f"{dist_matrix =}")
        #distance_matrix_triu = np.triu(dist_matrix)
        #print(distance_matrix_triu)
        #dist_matrix = np.array(dist_matrix)
        print(f"{dist_matrix =}")
        sequential = linkage(dist_matrix,method = "average")
        print(sequential)


matrice = matrice_score(id_sequence = id_sequence,gap = -5, blosum=blosum)
dist_matrix = matrice_dist(matrice)
embranchement_sequentiel(dist_matrix)


    
'''
    def UPGMA(matrice):
        print(matrice)
        clusters = [[i] for i in range(len(matrice))]
        print(clusters)
        while len(clusters) > 2 :
            max_i = value_max(matrice)[0]
            max_j = value_max(matrice)[1]
            new_cluster = clusters[max_i] + clusters[max_j]
            clusters.append(new_cluster)
            print(clusters)
            del clusters[max_i]
            print(f"first del : {clusters}")
            del clusters[max_j]
            print(f"second del : {clusters}")
            score = 0
            for i in range(len(clusters[0])):
                for j in range(len(clusters[1])):
                    seq1 = clusters[0][i]
                    seq2 = clusters[1][j]
                    print(matrice[seq1][seq2])
                    score += matrice[seq1][seq2]
            score = score / (len(clusters[0]) * len(clusters[1]))
            #matrice2 = matrice
            #matrice2 = matrice2.drop(columns = [max_j], inplace = True)
            #matrice2 = matrice2.drop(index = [max_i], inplace = True)
            print(matrice2)
            
            #matrice 
'''
    
       #clusters.append(clusters[max_i]+clusters[max_j])
       # tree = []
       # Permet de calculer la distance entre le groupe et l'arbre 
       # pos_last_ancestor = maximum / 2
       # tree.append()
       #print(matrice, labels,clusters)

    


#score = needlman(key_sequence1='MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH', key_sequence2='MVHLTAHHFGLWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLNALAHKYH', gap = - 5, blosum=blosum)
# score = needlman(id_sequence[1], id_sequence[0], gap = -5, blosum=blosum)
# print(score[0])

#alignement = alignement_seq(id_sequence[0], id_sequence[1], gap = -5, score = score[1], blosum=blosum)
#print("L'alignement optimal est : " + "\n","".join(alignement[0]), "".join(alignement[1]))






# max = UPGMA(matrice)



 