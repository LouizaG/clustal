import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio.Align import PairwiseAligner
from Bio.SubsMat import MatrixInfo
from itertools import combinations
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt


# Load the subtitution matrix : BLOSUM62
matrix = substitution_matrices.load("BLOSUM62")

def read_fasta_files(directory):
    """
    This function Reads all FASTA files needed for the multiple alignement.

    Parameters:
    directory (str): Path to the directory containing the FASTA files.

    Returns:
    dict: Dictionary with sequence identifiers as keys and sequences as values.
    """
    sequences = {}

    for filename in os.listdir(directory):
        if filename.endswith(".fasta"):
            with open(os.path.join(directory, filename), 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    sequences[record.id] = str(record.seq)

    return sequences

def compute_kmer_similarity(sequences, k):
    """
    Compute the k-mer similarity matrix for a set of sequences using BLOSUM62.

    Parameters:
    sequences (dict): A dictionary where keys are sequence IDs and values are sequences.
    k (int): The k-mer size.

    Returns:
    np.array: A k-mer similarity matrix.
    """
    def get_kmers(sequence, k):
        """Helper function to get all k-mers from a sequence."""
        return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

    # Create a list of all sequence IDs
    sequence_ids = list(sequences.keys())

    # Get the number of sequences
    n = len(sequence_ids)

    # Initialize an nxn matrix with zeros
    similarity_matrix = np.zeros((n, n))

    # Get the BLOSUM62 matrix
    blosum62 = MatrixInfo.blosum62

    # Fill the similarity matrix
    for i, j in combinations(range(n), 2):
        seq1_kmers = get_kmers(sequences[sequence_ids[i]], k)
        seq2_kmers = get_kmers(sequences[sequence_ids[j]], k)

        max_score = 0
        for kmer1 in seq1_kmers:
            for kmer2 in seq2_kmers:
                # Calculate similarity score using BLOSUM62 matrix
                score = sum(blosum62.get((aa1, aa2), blosum62.get((aa2, aa1), 0)) for aa1, aa2 in zip(kmer1, kmer2))
                max_score = max(max_score, score)
        
        similarity_matrix[i, j] = max_score
        similarity_matrix[j, i] = max_score  # Symmetric matrix

    # Fill the diagonal with the maximum possible similarity value (self-similarity)
    np.fill_diagonal(similarity_matrix, k * max(blosum62.values()))

    return pd.DataFrame(similarity_matrix, index=sequence_ids, columns=sequence_ids)

def compute_distance_matrix(similarity_matrix):
    """
    Compute the distance matrix from the similarity matrix.

    Parameters:
    similarity_matrix (pd.DataFrame): A DataFrame containing similarity scores.

    Returns:
    pd.DataFrame: A distance matrix.
    """
    # Normalize the similarity matrix
    max_sim = similarity_matrix.max().max()
    normalized_similarity_matrix = similarity_matrix / max_sim

    # Compute the distance matrix
    distance_matrix = 1 - normalized_similarity_matrix

    return distance_matrix

def compute_UPGMA_tree(distance_matrix, index_to_id_pair):
    """
    Compute the UPGMA tree using hierarchical clustering.

    Parameters:
    distance_matrix (np.ndarray): The upper triangular distance matrix.
    index_to_id_pair (dict): A dictionary mapping indices to sequence ID pairs.

    Returns:
    np.ndarray: The linkage matrix representing the hierarchical clustering.
    tuple: A nested tuple representing the hierarchical clustering structure.
    """
    # We need to convert the upper triangular matrix to a condensed distance matrix
    condensed_dist_matrix = distance_matrix[distance_matrix.nonzero()]
    
    # Perform hierarchical clusering to obtain the linkage matrix
    linkage_matrix = linkage(condensed_dist_matrix, method='average')

    # Create a list of labels based on the index_to_id_pair dictionary
    labels = [index_to_id_pair[(i, i)][0] for i in range(len(distance_matrix))]
    
    # Initialize clusters with individual sequences
    clusters = {i: labels[i] for i in range(len(labels))}
    
    # Build the nested tuple representation
    for i, (idx1, idx2, _, _) in enumerate(linkage_matrix):
        new_cluster_id = len(labels) + i
        clusters[new_cluster_id] = (clusters.pop(idx1), clusters.pop(idx2))
    
    # Get the final nested tuple representation
    nested_tuple_representation = list(clusters.values())[0]

    # Plot the UPGMA dendrogram
    dendro = dendrogram(linkage_matrix, labels=labels)
    plt.title('UPGMA Tree')
    plt.xlabel('Sequence ID')
    plt.ylabel('Distance')
    plt.show()
    
    return linkage_matrix, nested_tuple_representation


#---------------------------------#
#--------------MAIN---------------#
#---------------------------------#

# Read fasta files 
directory = "/Users/galou/M2/progra3/data"
sequences = read_fasta_files(directory)
# simiularity matrix with K-mers 
k = 3
similarity_matrix_df = compute_kmer_similarity(sequences, k)
print(similarity_matrix_df)
# Distance matrix 
distance_matrix_df = compute_distance_matrix(similarity_matrix_df)
print(distance_matrix_df)
# Compute and plot the UPGMA tree
linkage_matrix, nested_tuple_representation = compute_UPGMA_tree(distance_matrix_triu, index_to_id_pair)
print(linkage_matrix)
print(nested_tuple_representation)

