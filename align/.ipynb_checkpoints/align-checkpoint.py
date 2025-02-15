# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO

        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm

        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA

        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        if not seqA or not seqB:
            raise ValueError("sequences cannot be empty")
            
        # resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        # nrows is the length of seqA and the ncols is the length of seqB
        rows, columns = len(seqA) + 1, len(seqB) + 1 # the +1 is for the gap
        self._align_matrix = np.zeros((rows, columns)) # initializing an empty array

        # keeping track of gap penalties in seqA
        self._gapA_matrix = np.zeros((rows, columns))

        # keeping track of penalties in seqB
        self._gapB_matrix = np.zeros((rows, columns))

        # storing traceback directions to reconstruct the optimal alignment
        self._back = np.zeros((rows, columns), dtype=int)
        """
          traceback directions are not numeric scores; 0, 1, and 2 represent direction (diagonal, left, up -- respectively)
          np.zeros() creates a float matrix by default, so using dtype=int ensures that _back stores whole numbers for directions
        """

        # tracking gap directions in seqA and seqB
        self._back_A = np.zeros((rows, columns), dtype=int)
        self._back_B = np.zeros((rows, columns), dtype=int)

        # TODO: Implement global alignment here

        # filling self._align_matrix first row and col
        # first gap costs self.gap_open and following gaps cost self.gap_extend, top-left cell (0,0) remains 0
        
        # initialize first row and first column with gap penalties
        for i in range(1, rows):
            self._align_matrix[i, 0] = self.gap_open + (i - 1) * self.gap_extend
            self._back[i, 0] = 1
        for j in range(1, columns):
            self._align_matrix[0, j] = self.gap_open + (j - 1) * self.gap_extend
            self._back[0, j] = 2

        for index in range((rows - 1) * (columns - 1)):
            i = (index // (columns - 1)) + 1  # row index after skipping first row
            j = (index % (columns - 1)) + 1  # col index after skipping first col
        
            sub_score = self.sub_dict.get((self._seqA[i-1], self._seqB[j-1]), self.sub_dict.get((self._seqB[j-1], self._seqA[i-1]), 0))
            diagonal = self._align_matrix[i-1, j-1] + sub_score
            """
            if self._seqA[2], self._seqB[3] = (C, C) would then be looked up in self.sub_dict[('C', 'C')] for substitution score
            """
            left = self._align_matrix[i, j-1] + (self.gap_extend if self._back[i, j-1] == 2 else self.gap_open)
            up = self._align_matrix[i-1, j] + (self.gap_extend if self._back[i-1, j] == 1 else self.gap_open)

            """
            imagine two example seqs:
            SeqA =  A C G T G T C A T
            SeqB =  A G T C T
            seqB needs 4 gaps (somewhere) to match seqA, so the algorithm decides whether to open a new gap or extend an existing gap by calculating the more expensive option
            if NO gap existed before (i, j-1) --> self.gap_open, if a gap DID exist before (i, j-1) --> self.gap_extend
            """
            
            best_score = max(diagonal, left, up)
            
            if i == rows - 1 and j == columns - 1:
                best_score -= 1 
            
            self._align_matrix[i, j] = best_score
            

            if best_score == diagonal:
                self._back[i, j] = 0  # storing direction for traceback, diagonal = 0
            elif best_score == left:  
                self._back[i, j] = 2  # storing direction for traceback, left = 1
            else:
                self._back[i, j] = 1  # storing direction for traceback, up = 2


        # storing bottom-right cell as final alignment score
        self.alignment_score = self._align_matrix[-1, -1]

        # after filling self._align_matrix, align() stores traceback directions in self._back
        # after the matrix is fully built, align() calls _backtrace(), which performs the traceback
        # _backtrace() isn't explicitly called inside align() because it only needs to run after align() has completed matrix filling
        
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO

        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.

        Parameters:
        	None

        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        i, j = self._align_matrix.shape[0] - 1, self._align_matrix.shape[1] - 1  # starting at bottom-right corner of matrix
        self.seqA_align, self.seqB_align = "", ""  # empty strings for alignment
    
        while i > 0 or j > 0: # stop the loop at i == 0 or j == 0, this means one seq has been fully traced back and the other needs gaps
            if i > 0 and j > 0 and self._back[i, j] == 0:  # if diagonal
                self.seqA_align = self._seqA[i - 1] + self.seqA_align # add to seqA
                self.seqB_align = self._seqB[j - 1] + self.seqB_align # add to seqB
                i -= 1
                j -= 1
            elif i > 0 and self._back[i, j] == 1:  # if up
                if j > 0 and self._back[i-1, j] == 0:
                    self.seqA_align = self._seqA[i - 1] + self.seqA_align 
                    self.seqB_align = self._seqB[j - 1] + self.seqB_align
                    i -= 1
                    j -= 1
                else: # if left 
                    self.seqA_align = self._seqA[i - 1] + self.seqA_align
                    self.seqB_align = "-" + self.seqB_align
                    i -= 1
            elif j > 0 and self._back[i, j] == 2:  # Left move (gap in seqA)
                self.seqA_align = "-" + self.seqA_align
                self.seqB_align = self._seqB[j - 1] + self.seqB_align
                j -= 1

    
        return (self.alignment_score, self.seqA_align, self.seqB_align)

        
def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header