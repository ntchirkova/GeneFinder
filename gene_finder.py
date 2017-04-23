# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Nina Tchirkova

"""
from datetime import datetime
import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
dna = load_seq("./data/X73525.fa")

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        I think the tests are sufficient because the method is so simple.

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'



def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        The test is sufficient because all the method does is reverses the string and
        calls on get_complement(), there are no edge cases.

        Edited to change efficiency

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    r_list = []
    length = len(dna)
    for i in range(length):
        r_list.append(get_complement(dna[length - 1 - i]))
    return ''.join(r_list)



def rest_of_ORF(start, dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        These tests are sufficient because they test the two cases possible (no stop codon and stop codon)

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF(0,"ATGTGAA")
    'ATG'
    >>> rest_of_ORF(0,"ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF(0,"AAAAAAAT")
    'AAAAAAAT'
    >>> rest_of_ORF(0,"ATGATG")
    'ATGATG'
    """
    length = len(dna)
    i = start;
    stopindex = start;
    while i < length:
        #case when less than 3 nucleotides left
        if length - i < 3:
            return(dna[start:])
        #Checks if stop codon
        if dna[i] == 'T' and ((dna[i+1] == 'A' and (dna[i+2] == 'A' or dna[i+2] == 'G')) or (dna[i+1] == 'G' and dna[i+2] == 'A')):
            break
        else:
            i = i+3
            stopindex = stopindex+3
    return dna[start:stopindex]




def find_start(start, dna):
    """ returns index of first start codon (if there is no start codon returns -1)
        this is my own function for me

        dna: string of DNA
        return: index of the start codon

        >>> find_start(0,'ATGGGGGGGGGG')
        0
        >>> find_start(0,'GGGGGGATGAAA')
        6
        >>> find_start(0,'GGGGGG')
        -1
    """
    for i in range(start, len(dna), 3):
        if len(dna) - i < 3:
            break
        if dna[i] == 'A' and dna[i+1] == 'T' and dna[i+2] == 'G':
            return i
    return -1

def find_all_ORFs_oneframe(start, dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe(0,"ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe(0,'GGGGGGGG')
    []
    >>> find_all_ORFs_oneframe(0,'AATGGGGGGTGAAA')
    []
    """
    orfList = [] #list of orfs
    startIndex = start
    while True:
        startIndex = find_start(startIndex, dna)
        if startIndex == -1:
            break
        else:
            orfList.append(rest_of_ORF(startIndex, dna)) #cleaned orf is added to orf list
            startIndex = startIndex + len(orfList[-1])
    return orfList



def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("AATGCATGAATGTAG")
    ['ATG', 'ATGCATGAATGTAG', 'ATGAATGTAG']
    >>> find_all_ORFs("AAATGCATGAATGTAG")
    ['ATGAATGTAG', 'ATG', 'ATGCATGAATGTAG']
    """
    orfList = find_all_ORFs_oneframe(0, dna)
    orfList.extend(find_all_ORFs_oneframe(1, dna))
    orfList.extend(find_all_ORFs_oneframe(2, dna))
    return(orfList)



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("AATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCATT']
    >>> find_all_ORFs_both_strands("AAATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCATTT']
    """
    orfList = find_all_ORFs(dna)
    reverse = get_reverse_complement(dna)
    orfList.extend(find_all_ORFs(reverse))
    return orfList


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    orf_list = find_all_ORFs_both_strands(dna)
    longest = ''
    for orf in orf_list:
        if len(orf) > len(longest):
            longest = orf
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    length = 0
    for i in range(num_trials):
        dna = shuffle_string(dna)
        longest = longest_ORF(dna)
        if length < len(longest):
            length = len(longest)
    return length



def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    acids = ''
    #ensures that only codons are counted and no stragglers at the end
    for i in range(0,len(dna)-len(dna)%3,3):
        codon = dna[i:i+3]
        acids += aa_table[codon]
    return acids




def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """

    threshold = longest_ORF_noncoding(dna, 10)
    aa_list = []
    orf_list = find_all_ORFs_both_strands(dna)
    for orf in orf_list:
        if len(orf) > threshold:
            acid = coding_strand_to_AA(orf)
            aa_list.append(acid)
    return aa_list



if __name__ == "__main__":
    import doctest
    doctest.testmod()
    startTime = datetime.now()
    bb =gene_finder(dna)
    print(str(datetime.now() - startTime))
