# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Nina Tchirkova

"""

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

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse = ''
    length = len(dna)
    for i in range(length):
        reverse = reverse + get_complement(dna[length - 1 - i])
    return reverse



def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        These tests are sufficient because they test the two cases possible (no stop codon and stop codon)

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    length = len(dna)
    cleanDNA = ''
    for i in range (0, length, 3):
        codon = dna[i:i+3]
        if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
            break
        else:
            cleanDNA = cleanDNA + codon
    return cleanDNA


def find_start(dna):
    """ returns index of first start codon (if there is no start codon returns -1)
        this is my own function for me

        dna: string of DNA
        return: index of the start codon

        >>> find_start('ATGGGGGGGGGG')
        0
        >>> find_start('GGGGGGATGAAA')
        6
        >>> find_start('GGGGGG')
        -1
    """
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        if codon == 'ATG':
            return i
    return -1

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe('GGGGGGGG')
    []
    >>> find_all_ORFs_oneframe('AATGGGGGGTGAAA')
    []
    """
    orfList = [] #list of orfs
    while True:
        startIndex = find_start(dna)
        if startIndex == -1:
            break
        else:
            dna = dna[startIndex:] #cuts off all DNA before start codon
            orf = rest_of_ORF(dna)
            orfList.append(orf) #cleaned orf is added to orf list
            newStart = len(orf)
            dna = dna[newStart:] #prepares sequence to be tested for more start cod
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
    orfList = find_all_ORFs_oneframe(dna)
    orfList = orfList + find_all_ORFs_oneframe(dna[1:])
    orfList = orfList + find_all_ORFs_oneframe(dna[2:])
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
    orfList = orfList + find_all_ORFs(reverse)
    return orfList


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    orfList = find_all_ORFs_both_strands(dna)
    longest = ''
    for orf in orfList:
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
    for i in range(0,len(dna)-len(dna)%3,3):
        codon = dna[i:i+3]
        acids += aa_table[codon]
    return acids




def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    aa_list = []
    orfList = find_all_ORFs_both_strands(dna)
    for orf in orfList:
        if len(orf) > threshold:
            acid = coding_strand_to_AA(orf)
            aa_list.append(acid)

    return aa_list



if __name__ == "__main__":
    import doctest
    doctest.testmod()
    #print(gene_finder(dna))
    #print("done thanks Nina")
