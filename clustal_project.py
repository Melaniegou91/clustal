from Bio import SeqIO
import pandas as pd
import numpy as np
with open("sequences.fasta" , "r") as fasta_file : # Ouverture du fichier contenant les séquences
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file,"fasta")) # Creation d'un dictionnaire contenant les séquences
    #print(record_dict) # Affiche le contenu du dictionnaire
    blosum = pd.read_csv("blossum62.txt") # Permet de lire la matrice blosum 62
    id_sequence = list(record_dict.keys()) # Contient toutes les clés du dictionnaire record_dict
    #print(id_sequence)
    def needlman(key_sequence1,key_sequence2,gap,blosum):
        sequence1 = record_dict[key_sequence1]
        sequence2 = record_dict[key_sequence2]

        matrix= pd.DataFrame(index = sequence1,columns = sequence2) # Creation de la matrice vide de longueur n*m


        for i in range(1:len(sequence1)) :
            print(i)
        '''
        for i in matrix.index :
            value_init = 0
            matrix.index[i]=value_init + gap
            value_init = matrix.index[i]
        print(matrix)
        return matrix
        '''
        
    
    
    
    matrix = needlman(id_sequence[1],id_sequence[2],gap = 2, blosum=blosum)
    print(matrix)

#networkX
    
    



    # print(blosum["A"]["A"])
    



    '''
    with open("blossum62.txt","r") as blosum : 
        first_am = blosum.read()
       # first_am = first_am[0:-2]

        print(first_am)
   '''




