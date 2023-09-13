import os
import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices as sm
from itertools import combinations
from scipy.cluster.hierarchy import linkage, dendrogram

BLOSUM62 = sm.load("BLOSUM62")
KMERS = 3
GAP_PENALITY = -5
MATCHING_GAP_PENALITY = +4


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
    sequences (dict): A dictionary where keys are sequence IDs
    and values are sequences.
    k (int): The k-mer size.

    Returns:
    np.array: A k-mer similarity matrix.
    """
    def get_kmers(sequence, k):
        """Helper function to get all k-mers from a sequence."""
        return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]

    # Create a list of all sequence IDs
    sequence_ids = list(sequences.keys())

    # Get the number of sequences
    n = len(sequence_ids)

    # Initialize an nxn matrix with zeros
    similarity_matrix = np.zeros((n, n))

    # Fill the similarity matrix
    for i, j in combinations(range(n), 2):
        seq1_kmers = get_kmers(sequences[sequence_ids[i]], k)
        seq2_kmers = get_kmers(sequences[sequence_ids[j]], k)

        max_score = 0
        for kmer1 in seq1_kmers:
            for kmer2 in seq2_kmers:
                # Calculate similarity score using BLOSUM62 matrix
                score = sum(
                    BLOSUM62.get((aa1, aa2), BLOSUM62.get((aa2, aa1), 0))
                    for aa1, aa2 in zip(kmer1, kmer2)
                )
                max_score = max(max_score, score)

        similarity_matrix[i, j] = max_score
        similarity_matrix[j, i] = max_score  # Symmetric matrix

    # Fill the diagonal with the maximum possible similarity value
    # (self-similarity)
    np.fill_diagonal(similarity_matrix, k * max(BLOSUM62.values()))

    return similarity_matrix


def similarity_to_distance_matrix(similarity_matrix):
    """
    Convert a similarity matrix to a distance matrix.

    Parameters:
    similarity_matrix (np.array): The similarity matrix.

    Returns:
    np.array: The distance matrix.
    """
    distance_matrix = np.zeros(
        (len(similarity_matrix), len(similarity_matrix)))
    max_value = np.nanmax(similarity_matrix)

    for i in range(similarity_matrix.shape[0]):
        for j in range(i + 1, similarity_matrix.shape[1]):
            distance_matrix[i][j] = max_value - (similarity_matrix[i][j])

    return distance_matrix


def sequential_clustering(distance_matrix, labels):
    """
    Build a tree using sequential clustering.

    Parameters:
    distance_matrix (np.array): The distance matrix.
    labels (list): The labels for the sequences.

    Returns:
    tuple: The hierarchical tree represented as nested tuples.
    """
    # Copy the distance matrix to avoid modifying the original
    dist_mat = distance_matrix.copy()

    # Initialize the cluster with the labels of individual sequences
    clusters = [[label] for label in labels]

    while len(clusters) > 1:
        # Find the indices of the two closest clusters
        min_dist = np.inf
        min_index = (0, 0)
        for i in range(dist_mat.shape[0]):
            for j in range(i + 1, dist_mat.shape[1]):
                if dist_mat[i, j] < min_dist:
                    min_dist = dist_mat[i, j]
                    min_index = (i, j)

        # Merge the two closest clusters
        new_cluster = [clusters[min_index[0]], clusters[min_index[1]]]
        clusters.append(new_cluster)

        # Update the distance matrix with the distances of the new cluster
        new_dist_row = [
            (dist_mat[min_index[0], i] + dist_mat[min_index[1], i]) / 2
            for i in range(dist_mat.shape[0])
        ]
        new_dist_row = np.array(new_dist_row)[:, np.newaxis]
        dist_mat = np.hstack((dist_mat, new_dist_row))

        new_dist_col = np.append(new_dist_row, 0)[np.newaxis, :]
        dist_mat = np.vstack((dist_mat, new_dist_col))

        # Remove the rows and columns corresponding to the clusters that we
        # merged
        dist_mat = np.delete(dist_mat, [min_index[0], min_index[1]], axis=0)
        dist_mat = np.delete(dist_mat, [min_index[0], min_index[1]], axis=1)

        # Remove the merged clusters from our list of clusters
        clusters.pop(min_index[1])
        clusters.pop(min_index[0])

    return tuple(clusters[0])


def pairwise_alignment(seq_m, seq_n, traceback='no'):
    """
    Performs pairwise alignment using the Needleman-Wunsch algorithm.

    Parameters:
    seq_m, seq_n (str): The sequences to be aligned.
    traceback (str): Whether to perform traceback to
    get the alignment ('yes' or 'no', default is 'no').

    Returns:
    If traceback is 'no', returns the alignment score (int).
    If traceback is 'yes', returns the aligned sequences (list of str).
    """
    # The length of the sequences
    m = len(seq_m)
    n = len(seq_n)

    # Initialize the score matrix with zeros
    scores = np.zeros((m + 1, n + 1))

    # Fill the first column and the first row with gap penalties
    scores[:, 0] = [GAP_PENALITY * i for i in range(m + 1)]
    scores[0, :] = [GAP_PENALITY * i for i in range(n + 1)]

    # Fill in the score matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = scores[i - 1,
                           j - 1] + BLOSUM62.get((seq_m[i - 1],
                                                  seq_n[j - 1]),
                                                 BLOSUM62.get((seq_n[j - 1],
                                                               seq_m[i - 1]),
                                                              None))
            delete = scores[i - 1, j] + GAP_PENALITY
            insert = scores[i, j - 1] + GAP_PENALITY
            scores[i, j] = max(match, delete, insert)

    alignment_score = scores[m, n]
    if traceback == 'no':
        return alignment_score
    else:

        # traceback to find the optimal alignment

        align_m = ""
        align_n = ""

        # starting point in the bottom right cell in the matrix
        i = m
        j = n

        while i > 0 and j > 0:
            current_score = scores[i, j]
            diagonal_score = scores[i - 1, j - 1]
            upwards_score = scores[i, j - 1]
            left_score = scores[i - 1, j]

            # check which cell the current score comes from, update i and j
            # accordingly
            if current_score == (diagonal_score +
                                 BLOSUM62.get((seq_m[i - 1], seq_n[j - 1]),
                                              BLOSUM62.get((seq_n[j - 1], seq_m[i - 1]), None))
                                 ):

                align_m += seq_m[i - 1]
                align_n += seq_n[j - 1]
                i -= 1
                j -= 1
            elif current_score == upwards_score + GAP_PENALITY:
                align_n += seq_n[j - 1]
                align_m += '-'
                j -= 1
            elif current_score == left_score + GAP_PENALITY:
                align_n += '-'
                align_m += seq_m[i - 1]
                i -= 1

        while j > 0:
            align_n += seq_n[j - 1]
            align_m += '-'
            j -= 1
        while i > 0:
            align_n += '-'
            align_m += seq_m[i - 1]
            i -= 1

        align_m = align_m[::-1]
        align_n = align_n[::-1]

        aligned_seq = [align_m, align_n]

        return aligned_seq


def flatten_tree(tree):
    """
    Recursively flattens a nested list (tree) into a single list.

    The function takes a nested list (tree) where each element can
    either be a value or another nested list. It flattens this structure
    into a single list of values, maintaining the order of the elements.

    Parameters:
    tree (list): A nested list structure to be flattened.
    The nested list can contain other lists or non-list values.

    Returns:
    list: A single list containing all the values from the input nested list,
    in the order they appeared in the nested structure.
    """

    flat_list = []
    for item in tree:
        if isinstance(item, list):
            flat_list.extend(flatten_tree(item))
        else:
            flat_list.append(item)
    return flat_list


def compute_multiple_alignment(sequence_dict, tree):
    """
    Perform a multiple sequence alignment using progressive alignment.

    This function takes a dictionary of sequences and a guiding tree
    representing the alignment order. It progressively aligns the sequences
    according to the tree structure, utilizing pairwise alignment
    and a scoring system defined by gap penalties and a BLOSUM62 matrix.

    Parameters:
    - sequence_dict (dict): A dictionary where keys are sequence identifiers
    and values are sequences to be aligned.
    - tree (list): A nested list representing the hierarchical
    structure in which sequences should be aligned.

    Returns:
    - list: A list of strings representing the aligned sequences.

    """
    def alignment_score(algt_list, seq_2, position_i, position_j):
        """
        Compute the alignment score at a specific position based
        on a list of aligned sequences and a new sequence.

        The function iterates over all pairs of amino acids at the given
        positions in the sequences, and accumulates a score based on a
        scoring matrix and penalties for matching gaps and other gaps.

        Parameters:
        - algt_list (list): A list of strings representing the aligned sequences up until this point.
        - seq_2 (str): A string representing the new sequence to be aligned.
        - position_i (int): The index representing the current position in the aligned sequences.
        - position_j (int): The index representing the current position in the new sequence to be aligned.

        Returns:
        - int: The alignment score for the current positions based on the scoring system.
        """
        score = 0

        for seq in algt_list:
            if seq[position_i] == '-' and seq_2[position_j] == '-':
                score += MATCHING_GAP_PENALITY
            elif seq[position_i] == '-' or seq_2[position_j] == '-':
                score += GAP_PENALITY
            else:
                pair = (seq[position_i], seq_2[position_j])
                score += BLOSUM62.get(pair, BLOSUM62.get(pair[::-1], 0))

        for i in range(len(algt_list)):
            for j in range(i + 1, len(algt_list)):
                pair = (algt_list[i][position_i], algt_list[j][position_i])
                if pair[0] == '-' and pair[1] == '-':
                    score += MATCHING_GAP_PENALITY
                elif pair[0] == '-' or pair[1] == '-':
                    score += GAP_PENALITY
                else:
                    score += BLOSUM62.get(pair, BLOSUM62.get(pair[::-1], 0))

        return score

    def get_optimal_alignment(scores, alignment_list, seq_n):
        """
        Finds the optimal alignment using a scoring matrix and a list of alignments.
        This function backtracks through a scoring matrix to find the optimal alignment
        of a list of aligned sequences with a new sequence.

        Parameters:
        - scores (numpy.ndarray): A 2D array containing the alignment
        scores at each position.
        - alignment_list (list of str): A list of strings representing the
        aligned sequences up until this point.
        - seq_n (str): A string representing the new sequence to be aligned.

        Returns:
        - list of str: A list of strings representing the optimal alignment of
        all sequences including the new sequence.
        """
        m, n = scores.shape[0] - 1, scores.shape[1] - 1
        align_m = [''] * len(alignment_list)
        align_n = ''

        while m > 0 and n > 0:
            score = scores[m, n]
            score_diag = scores[m - 1, n - 1]
            score_up = scores[m, n - 1]
            score_left = scores[m - 1, n]

            if score == score_diag + \
                    alignment_score(alignment_list, seq_n, m - 1, n - 1):
                for k, seq in enumerate(alignment_list):
                    align_m[k] += seq[m - 1]
                align_n += seq_n[n - 1]
                m -= 1
                n -= 1
            elif score == score_up + GAP_PENALITY:
                for k in range(len(align_m)):
                    align_m[k] += '-'
                align_n += seq_n[n - 1]
                n -= 1
            else:
                for k, seq in enumerate(alignment_list):
                    align_m[k] += seq[m - 1]
                align_n += '-'
                m -= 1

        # Finalizing the alignment with remaining gaps if necessary
        while n > 0:
            align_n += seq_n[n - 1]
            for seq in align_m:
                seq += '-'
            n -= 1
        while m > 0:
            align_n += '-'
            for k, seq in enumerate(alignment_list):
                align_m[k] += seq[m - 1]
            m -= 1

        # Reversing the alignments to get the final sequences
        align_m = [seq[::-1] for seq in align_m]
        align_n = align_n[::-1]

        return align_m + [align_n]

    flat_tree = list(flatten_tree(tree))
    alignment_list = pairwise_alignment(
        sequence_dict[flat_tree[0]],
        sequence_dict[flat_tree[1]],
        traceback='yes'
    )

    for a in range(2, len(flat_tree)):
        seq_n = sequence_dict[flat_tree[a]]
        m, n = len(max(alignment_list, key=len)), len(seq_n)
        scores = np.zeros((m + 1, n + 1))
        scores[:, 0] = np.arange(m + 1) * GAP_PENALITY
        scores[0, :] = np.arange(n + 1) * GAP_PENALITY

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = scores[i - 1, j - 1] + \
                    alignment_score(alignment_list, seq_n, i - 1, j - 1)
                gap1 = scores[i - 1, j] + GAP_PENALITY
                gap2 = scores[i, j - 1] + GAP_PENALITY
                scores[i, j] = max(match, gap1, gap2)

        alignment_list = get_optimal_alignment(scores, alignment_list, seq_n)

    return alignment_list


def print_multiple_alignement(aligned_sequences, sequence_ids, filename):
    """
    Print the multiple sequence alignment and save it to a txt file.

    Parameters:
    - aligned_sequences (list): A list of strings representing the aligned sequences.
    - sequence_ids (list): A list of strings representing the identifiers for each sequence.
    - filename (str): The name of the file where the aligned sequences should be saved.

    Returns:
    - None
    """
    # Get the number of sequences
    num_sequences = len(aligned_sequences)

    max_id_length = max(len(id) for id in sequence_ids)

    with open(filename, 'w') as file:
        for i in range(num_sequences):
            line = f"{sequence_ids[i].rjust(max_id_length)}: {aligned_sequences[i]}"

            file.write(line + "\n")

            print(line)


def main():
    directory = "/home/CLUSTAL/data"
    sequences = read_fasta_files(directory)
    similarity_matrix = compute_kmer_similarity(sequences, KMERS)

    distance_matrix = similarity_to_distance_matrix(similarity_matrix)

    # Create a list of labels
    labels = list(sequences.keys())
    # get the guiding tree from sequential clustering
    tree = sequential_clustering(distance_matrix, labels)

    flat_list = flatten_tree(tree)
    # print(flat_list)
    # get the sequences IDs
    sequence_ids = list(sequences.keys())

    # Get the aligned sequences
    aligned_sequences = compute_multiple_alignment(sequences, tree)

    # Save the aligned sequences with their ids to a file and print them to
    # the console
    print_multiple_alignement(
        aligned_sequences,
        sequence_ids,
        "aligned_sequences.txt")


if __name__ == "__main__":
    main()
