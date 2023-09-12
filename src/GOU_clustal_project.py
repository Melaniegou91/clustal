from Bio import SeqIO
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--file_name", help="give the file name")
parser.add_argument("--file_matrix_blosum", help="blosum 62 matrix")
parser.add_argument("--seq1",
                    help="first sequence",
                    type=int)
parser.add_argument("--seq2",
                    help="second sequence",
                    type=int)
args = parser.parse_args()
file = args.file_name
blosum = args.file_matrix_blosum
choice_seq1 = args.seq1
choice_seq2 = args.seq2


with open(file, "r") as fasta_file:  # Ouverture du fichier
    # Dictionary containing the sequences
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # reading blosum 62 matrix
    blosum = pd.read_csv(blosum)

    # List containing sequence keys
    id_sequence = list(record_dict.keys())

    #
    print(f"Votre fichier comporte {len(id_sequence)} séquences")

    def needleman(key_sequence1, key_sequence2, gap, blosum):
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
        blosum : pd.dataframe
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
                diagonale = matrice[i - 1][j - 1] + \
                    blosum[sequence1[i-1]][sequence2[j-1]]
                gap_haut = matrice[i - 1][j] + gap
                gap_gauche = matrice[i][j - 1] + gap
                # Allows you to select the maximum between the 3 possibilities
                matrice[i][j] = max(diagonale, gap_haut, gap_gauche)
        return matrice[i][j], matrice

    def matrice_score(id_sequence, gap, blosum):
        '''
        This algorithm fills the similarity score matrix
        by filling a square matrix
        Each box is filled using the Needleman algorithm

        Parameters
        ----------
        id_sequence : list
            The list containing the id of the sequences
        gap : int
            The penality of a gap
        blosum : pd.dataframe
            The blosum 62 matrix

    Returns
    -------
        returns the score matrix
        '''
        # Creating a square matrix filled with zeros
        score_matrix = np.zeros((len(id_sequence), len(id_sequence)))
        # Browse the boxes of the matrix to compare
        # the different sequences
        for i in range(1, len(id_sequence)):
            for j in range(0, i):
                # Assigns for each sequence comparison, a similarity score
                score_matrix[i][j] = needleman(id_sequence[i],
                                               id_sequence[j],
                                               gap=gap,
                                               blosum=blosum)[0]
        print(f"La matrice de score est : \n {score_matrix}")
        return score_matrix

    def dist_matrix(score_matrix):
        '''

        This method makes it possible
         to transform a similarity score matrix
         into a distance matrix

        Parameters
        ----------
        score_matrix : np.array
            The score matrix to transform


    Returns
    -------
        returns the distance matrix
        '''
        matrix_dist = np.zeros((len(score_matrix), len(score_matrix)))

        # Allows you to find the maximum value in the matrix
        maximum = np.nanmax(score_matrix)

        # Allows you to find the minimum value in the matrix
        minimum = np.nanmin(score_matrix)
        # Allows you to loop through all score matrix values
        # ​​to transform them into a distance score
        for i in range(1, len(id_sequence)):
            for j in range(0, i):
                matrix_dist[i][j] = maximum - (score_matrix[i][j] + minimum)
        return matrix_dist

    def embranchement_sequentiel(dist_matrix):
        """
        This method makes it possible to create the phylogenetic tree

        Parameters
        ----------
        dist_matrix : np.array
            The distance matrix used to make the tree


    Returns
    -------
        returns the phylogenetic tree

        """
        print(f"La matrice de distance est : \n \n {dist_matrix }")
        sequential = linkage(dist_matrix, method="average")
        print(sequential)
        phylogenetic_tree = dendrogram(sequential)
        plt.title('phylogenetic tree')
        plt.xlabel('Sequence ID')
        plt.ylabel('Distance')
        plt.show()

    #

    def alignement_seq(key_sequence1, key_sequence2, gap, score, blosum):
        '''
        The algorithm allows you to start
        from the lowest point of the matrix
        on the right to go up using the path
        worth the highest score

        Parameters
        ----------
        key_sequence1 : float
            The keys to the first sequence to study .
        key_sequence2 : float
            The keys to the second sequence to study.
        gap : int
            The penality of a gap
        blosum : pd.DataFrame
            the blosum 62 matrix

    Returns
    -------
        returns the alignment in pairs
        '''
        print(score)
        # Retrive the sequences
        seq1 = record_dict[key_sequence1].seq
        print(seq1)
        seq2 = record_dict[key_sequence2].seq
        print(seq2)
        # i and j contain the length of the sequences
        i = len(seq1)
        j = len(seq2)
        # Empty sequence initialization
        alignement1 = ""
        alignement2 = ""

        while i > 1 and j > 1:
            # Lets you know if the score comes from the upper left box
            if (score[i][j] == score[i - 1][j - 1] + blosum[seq1[i - 1]][seq2[j - 1]]): 
                i = i - 1
                j = j - 1
                alignement1 += seq1[i - 1]
                alignement2 += seq2[j - 1]

            # Allows you to know if the score comes from the box above
            elif (score[i][j] == score[i][j - 1] + gap):
                j = j - 1
                alignement1 += "-"
                alignement2 += seq2[j - 1]

            # Allows you to know if the score comes from the box below
            else:
                i = i - 1
                alignement1 += seq1[i - 1]
                alignement2 += "-"
        while i > 0:
            i -= 1
        while j > 0:
            j -= 1
            """if i > 0:
                i -= 1
            if j > 0:
                j -= 1"""
        alignement1 = "".join(reversed(alignement1))
        alignement2 = "".join(reversed(alignement2))
        nb_gap = alignement2.count("-")/len(alignement2)
        print(f"le nombre de gap dans l'alignement 2 ets : {nb_gap}")

        print(f"{alignement1=}")
        print(f"{alignement2=}")
        return alignement1, alignement2

matrice = matrice_score(id_sequence=id_sequence, gap=-5, blosum=blosum)
dist_matrix = dist_matrix(matrice)
embranchement_sequentiel(dist_matrix)
score = needleman(id_sequence[choice_seq1],
                  id_sequence[choice_seq2], gap=-5, blosum=blosum)
alignement_seq(id_sequence[choice_seq1], id_sequence[choice_seq2],
               gap=-5, score=score[1], blosum=blosum)


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
