# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    alignment_scores = {
        "Gallus gallus": nw.align(gg_seq, hs_seq)[0], # nw.align() returns (self.alignment_score, self.seqA_align, self.seqB_align)
        "Mus musculus": nw.align(mm_seq, hs_seq)[0],  # [0] extracts only the alignment score from the tuple
        "Balaeniceps rex": nw.align(br_seq, hs_seq)[0],
        "Tursiops truncatus": nw.align(tt_seq, hs_seq)[0]
    }
    ranked_species = sorted(alignment_scores.items(), key=lambda x: x[1], reverse=True) # sort based on alignment score, then flip (default ascending --> descending; high score first)
    """
    sorted() usually ranks by x[0], the first element in the list, so in alignment_scores.items() = (species, alignment_score), x[1] = alignment_score
    lambda x: x[1] is shorthand for:
    def get_second_element(x):
        return x[1]
    sorted(alignment_scores, key=get_second_item, reverse=True) is equivalent to sorted(alignment_scores, key=lambda x: x[1], reverse=True)
    lambda is just an inline way to define a function without using def
    """

    print("\nSpecies ranked by similarity to human BRD2:")
    for species, score in ranked_species:
        print(f"{species}: {score}")
    

if __name__ == "__main__":
    main()