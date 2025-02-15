# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    # initialize the NeedlemanWunsch object
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    # performing alignment between test seq1 and seq2 
    nw.align(seq1, seq2)

    # check that the alignment matrix is filled correctly
    assert nw._align_matrix is not None, "alignment matrix not filled"
    assert isinstance(nw._align_matrix, np.ndarray), "alignment matrix not a numpy array"

    # check that the gap matrices are filled correctly
    assert nw._gapA_matrix is not None, "gap matrix for seqA is not filled"
    assert nw._gapB_matrix is not None, "gap matrix for seqB is not filled"
    assert isinstance(nw._gapA_matrix, np.ndarray), "Gap matrix for seqA is not a numpy array"
    assert isinstance(nw._gapB_matrix, np.ndarray), "Gap matrix for seqB is not a numpy array"

    # check that the backtrace matrices are initialized
    assert nw._back is not None, "backtrace matrix is not initialized"
    assert nw._back_A is not None, "backtrace matrix for seqA is not initialized"
    assert nw._back_B is not None, "backtrace matrix for seqB is not initialized"
    assert isinstance(nw._back, np.ndarray), "backtrace matrix is not a numpy array"
    assert isinstance(nw._back_A, np.ndarray), "backtrace matrix for seqA is not a numpy array"
    assert isinstance(nw._back_B, np.ndarray), "backtrace matrix for seqB is not a numpy array"


def test_nw_empty_sequence():
    """
    ensures that aligning an empty sequence raises a ValueError
    """
    seq1 = "ACTG"  # example non-empty sequence
    empty_seq = ""  

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

  
    with pytest.raises(ValueError, match="sequences cannot be empty"):
        nw.align(seq1, empty_seq)



def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    score, aligned_seqA, aligned_seqB = nw.align(seq3, seq4)
    
    assert score == 17, f"Expected alignment score of 17, but got {score}"
    expected_seqA = "MAVHQLIRRP"
    expected_seqB = "M---QLIRHP"
    assert aligned_seqA == expected_seqA, f"Expected aligned sequence A to be {expected_seqA}, but got {aligned_seqA}"
    assert aligned_seqB == expected_seqB, f"Expected aligned sequence B to be {expected_seqB}, but got {aligned_seqB}"
    print(aligned_seqB)




