import os
import shutil
import numpy as np

# create folder
def create_folder(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Folder '{folder_path}' created.")
    else:
        print(f"Folder '{folder_path}' already exists.")

        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            if os.path.isdir(file_path):
                shutil.rmtree(file_path)
            else:
                os.remove(file_path)
        print(f"All contents of '{folder_path}' have been removed.")

# Needleman_Wunsch
def needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-2):

    n = len(seq1) + 1
    m = len(seq2) + 1
    score_matrix = [[0] * m for _ in range(n)]
    traceback_matrix = [[None] * m for _ in range(n)]

    for i in range(n):
        score_matrix[i][0] = i * gap_penalty
        traceback_matrix[i][0] = 'up'
    for j in range(m):
        score_matrix[0][j] = j * gap_penalty
        traceback_matrix[0][j] = 'left'
    traceback_matrix[0][0] = 'done'

    for i in range(1, n):
        for j in range(1, m):
            if seq1[i-1] == seq2[j-1]:
                diag = score_matrix[i-1][j-1] + match_score
            else:
                diag = score_matrix[i-1][j-1] + mismatch_penalty
            up = score_matrix[i-1][j] + gap_penalty
            left = score_matrix[i][j-1] + gap_penalty

            max_score = max(diag, up, left)
            score_matrix[i][j] = max_score

            if max_score == diag:
                traceback_matrix[i][j] = 'diag'
            elif max_score == up:
                traceback_matrix[i][j] = 'up'
            else:
                traceback_matrix[i][j] = 'left'

    align1, align2 = '', ''
    i, j = len(seq1), len(seq2)

    while traceback_matrix[i][j] != 'done':
        if traceback_matrix[i][j] == 'diag':
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 'up':
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        elif traceback_matrix[i][j] == 'left':
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1

    return align1, align2, score_matrix[-1][-1]


# Smith-Waterman
def smith_waterman(seq1, seq2, match_score=2, mismatch_penalty=-1, gap_penalty=-1):

    m, n = len(seq1), len(seq2)
    score_matrix = np.zeros((m+1, n+1), dtype=int)
    max_score = 0
    max_pos = None

    for i in range(1, m+1):
        for j in range(1, n+1):
            if seq1[i-1] == seq2[j-1]:
                diag = score_matrix[i-1][j-1] + match_score
            else:
                diag = score_matrix[i-1][j-1] + mismatch_penalty
            up = score_matrix[i-1][j] + gap_penalty
            left = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(0, diag, up, left)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    align1, align2 = '', ''

    i, j = max_pos
    while score_matrix[i][j] != 0:
        current_score = score_matrix[i][j]
        diagonal_score = score_matrix[i-1][j-1]
        up_score = score_matrix[i-1][j]
        left_score = score_matrix[i][j-1]

        if current_score == diagonal_score + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif current_score == up_score + gap_penalty:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        elif current_score == left_score + gap_penalty:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1

    return align1, align2, max_score


# similarity
def calculate_identity(align1, align2, method="alignment_length"):
    matches = 0
    aligned_length = len(align1)
    non_gap_positions = 0

    for a, b in zip(align1, align2):
        if a != '-' and b != '-':
            non_gap_positions += 1
            if a == b:
                matches += 1

    if method == "alignment_length":
        return matches / aligned_length * 100
    elif method == "non_gap":
        return matches / non_gap_positions * 100
    else:
        raise ValueError("method must be 'alignment_length' or 'non_gap'")


def needleman_wunsch_similarity(seq1, seq2):
    align1, align2, score = needleman_wunsch(seq1, seq2)
    similarity = calculate_identity(align1, align2)
    return similarity


def smith_waterman_similarity(seq1, seq2):
    align1, align2, score = smith_waterman(seq1, seq2)
    similarity = calculate_identity(align1, align2)
    return similarity

