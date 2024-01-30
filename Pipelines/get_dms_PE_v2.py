# Determine number of mutations and read coverage at every nucleotide in an aligned .bam file
# optimized for parallel processing
# Leonard Schaerfen
# 8/2/22

import numpy as np
import pickle
import time
from Bio import SeqIO
from multiprocessing import Pool
import pysam
import re
import sys
from collections import defaultdict

def coverage_mutations_PE(file, genome_file, n_cpu=4, clip_filter=99999, clip_3=0, return_read_stats=False):
    # parallel wrapper for coverage_mutations_chrom().

    start = time.time()

    # load the genome    
    chroms = []
    for i in SeqIO.parse(genome_file, 'fasta'):
        chroms.append(i)
    
    vec_chroms = []
    
    # using istarmap, a patch for starmap to allow for tqdm to work
    # get coverage and mutation count one chromosome per thread
    vars_iterable = [(file, chrom_nr, chroms[chrom_nr], clip_filter, clip_3, return_read_stats) for chrom_nr in range(len(chroms))]
    with Pool(processes=n_cpu) as pool:
        for vec in pool.starmap(coverage_mutations_chrom_PE, vars_iterable):
            vec_chroms.append(vec)

    end = time.time()
    print('time elapsed: %.2f min'%((end-start)/60))
    # unpack the individual vectors
    if not return_read_stats:
        return([i[0] for i in vec_chroms], [i[1] for i in vec_chroms], [i[2] for i in vec_chroms], [i[3] for i in vec_chroms])
    else:
        read_stats = []
        for i in vec_chroms:
            read_stats += i[4]
        return([i[0] for i in vec_chroms], [i[1] for i in vec_chroms], [i[2] for i in vec_chroms], [i[3] for i in vec_chroms], read_stats)

def coverage_mutations_chrom_PE(file, chrom_nr, chrom, clip_filter=9999999, clip_3=0, return_read_stats=False):
    # same as coverage_mutations(), for the parallel loop, runs on a single chromosome and it doesn't return the arrays because that would
    # cost too much RAM. Bitvectors are added up on-the-fly.

    # open bam file and get chromosome names
    bamfile = pysam.AlignmentFile(file, 'rb')

    chrom_names = bamfile.header.references
    
    # initialize the bitvector for the specified chromosome
    vec_mut_p = np.zeros(len(chrom), dtype=int)
    vec_cov_p = np.zeros(len(chrom), dtype=int)
    vec_mut_m = np.zeros(len(chrom), dtype=int)
    vec_cov_m = np.zeros(len(chrom), dtype=int)
    read_stats = []

    dms_dict = {'A':1, 'G':1, 'T':1, 'C':1, '0':0, '?':0, '.':0, '1':1, 'N': 0} # specify the value of bit entries         
    
    for read1, read2 in read_pair_generator(bamfile, reference=chrom_names[chrom_nr]):

        if not read1.is_proper_pair or read1.is_secondary or read1.is_supplementary:
            continue

        bv1, clipped_3_1, clipped_5_1 = Convert_Read(read1, with_clipped=True)
        bv2, clipped_3_2, clipped_5_2 = Convert_Read(read2, with_clipped=True)
        
        # skip reads that have a specified amount of clipped bases elsewhere
        if (len(clipped_3_1) > clip_filter) or (len(clipped_5_1) > clip_filter) or (len(clipped_5_2) > clip_filter) or (len(clipped_3_2) > clip_filter):
            continue 

        bitvector = join_reads(bv1, bv2)
        
        n_muts = np.sum(np.isin(np.array(list(bitvector.values())), ['A', 'G', 'C', 'T', '1']))
        read_stats.append(n_muts)

        # max. 10% of the read can be mutated
        if n_muts > 0.1 * len(bitvector):
            continue

        # PLUS strand
        if not read1.is_reverse:
            # remove specified number of bases from the 3'-end (because Pol protected)
            if clip_3 > 0:
                # find max bp value and remove specified number of previous bases
                end3 = max(bitvector.keys())
                [bitvector.pop(i, None) for i in np.arange(end3-(clip_3-len(clipped_5_2)), end3+1)]
            for pos in bitvector:
                if not bitvector[pos] == '.': # don't need to care about clipped bases for now
                    vec_mut_p[pos] += dms_dict[bitvector[pos]]
                    vec_cov_p[pos] += 1
        # MINUS strand
        else:
            # remove specified number of bases from the 3'-end (because Pol protected)
            if clip_3 > 0:
                # find max bp value and remove specified number of previous bases
                end3 = min(bitvector.keys())
                [bitvector.pop(i, None) for i in np.arange(end3, end3+(clip_3-len(clipped_5_2))+1)]
            for pos in bitvector:
                if not bitvector[pos] == '.': # don't need to care about clipped bases for now
                    vec_mut_m[pos] += dms_dict[bitvector[pos]]
                    vec_cov_m[pos] += 1
        
    bamfile.close()

    if not return_read_stats:
        return(vec_mut_p, vec_cov_p, vec_mut_m, vec_cov_m)
    else:
        return(vec_mut_p, vec_cov_p, vec_mut_m, vec_cov_m, read_stats)


def join_reads(read1, read2):

    #############################################################
    # FORBID ILLEGAL OVERLAP DUE TO TRIMMING OF READ1 AND READ2 #
    #############################################################
    
    miss_info, ambig_info = '.', '?'
    nomut_bit, del_bit = '0', '1'
    bases = ['A', 'T', 'G', 'C']
    
    for pos in read2.keys():
        # non-overlapping nts
        if pos not in read1: # add new read2 nts
            read1[pos] = read2[pos]
        # overlapping nts
        else:
            if read1[pos] == read2[pos]:
                continue
            elif read1[pos] not in [miss_info, ambig_info] and read2[pos] in [miss_info, ambig_info]:
                continue
            elif read1[pos] in [miss_info, ambig_info] and read2[pos] not in [miss_info, ambig_info]:
                read1[pos] = read2[pos]
            elif read1[pos] in [miss_info, ambig_info] and read2[pos] in [miss_info, ambig_info]:
                read1[pos] = ambig_info
            else: # disagreement
                read1[pos] = ambig_info
    return(read1)



def Convert_Read(mate, qscore_cutoff=20, sur_bases=10, with_clipped=False):
    """
    Convert a read's sequence to a bit vector of 0s & 1s and substituted bases
    Args:
        mate (Mate): Read (pysam.AlignedSegment() object)
        refs_seq (dict): Sequences of the ref genomes in the file
        phred_qscore (dict): Qual score - ASCII symbol mapping
    Returns:
        bitvector_mate (dict): Bitvector. Format: d[pos] = bit

    Original code from https://codeocean.com/capsule/6175523/tree/v1

    Refactored completely to use pysam and consider additional scenarios
    """

    # 3'-end clipped bases are marked in the bit vector, 5'-end bases are not

    if mate.is_unmapped:
        return('')

    bitvector_mate = {}  # Mapping of read to 0s and 1s
    clipped_3 = ''
    clipped_5 = ''

    # if the read is on the - strand we have to get the complementary base
    if mate.is_reverse:
        reverse = True
        read_seq = reverse_complement(mate.get_forward_sequence()).upper()
    else:
        reverse = False
        read_seq = mate.get_forward_sequence().upper()  # Sequence of the read
    ref_seq = mate.get_reference_sequence().upper()  # Sequence of the ref genome
    q_scores = mate.get_forward_qualities()  # Qual scores of the bases in the read
    ref_pos = fix_ref_pos(mate.get_reference_positions(full_length=True))  # Pos in the ref sequence
    
    # print(ref_pos)

    miss_info, ambig_info = '.', '?'
    del_bit = '1'
    nomut_bit = '0'

    i = 0  # Pos in the ref sequence
    j = 0  # Pos in the read sequence
    l = 0  # Pos in the ref position list
    CIGAR_Ops = CIGAR_Ops =re.findall(r'(\d+)([A-Z]{1})', mate.cigarstring)
    # print(CIGAR_Ops)
    op_index = 0
    while op_index < len(CIGAR_Ops):  # Each CIGAR operation
        op = CIGAR_Ops[op_index]
        desc, length = op[1], int(op[0])

        if desc == 'M':  # Match or mismatch
            for _ in range(length):  # Each base
                if q_scores[j] >= qscore_cutoff:
                    if not reverse: # register the base that was mutated, depending on strand
                        bitvector_mate[ref_pos[l]] = ref_seq[i] if read_seq[j] != ref_seq[i] else nomut_bit
                    else:
                        bitvector_mate[ref_pos[l]] = ref_seq[i] if read_seq[j] != ref_seq[i] else nomut_bit
                else:  # < Qscore cutoff
                    bitvector_mate[ref_pos[l]] = ambig_info
                i += 1  # Update ref index
                j += 1  # Update read index
                l += 1  # Update position index

        elif desc == 'D': # Deletion
            if ref_pos[l-1] is None:  # if insertion is followed directly by deletion
                break

            for k in range(length - 1):  # All bases except the 3' end
                bitvector_mate[ref_pos[l-1]+k+1] = ambig_info
                i += 1  # Update ref index
            # 3' end of deletion
            ambig = Calc_Ambig_Reads(ref_seq, i, length, sur_bases)
            bitvector_mate[ref_pos[l-1]+length] = ambig_info if ambig else del_bit
            i += 1  # Update ref index

        elif desc == 'I':  # Insertion
            j += length  # Update read index
            l += length  

        elif desc == 'S':  # Soft clipping
            if (not reverse and op_index == len(CIGAR_Ops) - 1) or (reverse and op_index == 0):  # Soft clipped at the 3'-end
                for _ in range(length):
                    bitvector_mate[ref_pos[l]] = miss_info
                    l += 1  # Update position index (soft-clipped bases are part of the position list)
                clipped_3 = read_seq[j:j+length] # clipped 3'-end bases for studying poly(A) tail properties
            else: # Soft clipped at 5'-end
                for _ in range(length):
                    bitvector_mate[ref_pos[l]] = miss_info
                    l += 1  # Update position index (soft-clipped bases are part of the position list)
                clipped_5 = read_seq[j:j+length] # clipped 5'-end bases for studying RT-stop tail properties
            j += length  # Update read index

        elif desc == 'N': # Intron, already accounted for by pysam in the ref seq
            pass

        elif desc == 'H': # hard clip, already accounted for by minimap2
            pass

        else:
            print('Unknown CIGAR op encountered: %s'%(desc))
            break

        op_index += 1

    # return the vector and clipped bases if desired
    if with_clipped:
        if reverse:
            return(bitvector_mate, reverse_complement(clipped_3), reverse_complement(clipped_5))
        else:
            return(bitvector_mate, clipped_3, clipped_5)
    else:
        return(bitvector_mate)

def Calc_Ambig_Reads(ref_seq, i, length, num_surBases):
    """
    Determines whether a deletion is ambiguous or not by looking at the
    sequence surrounding the deletion. Edge cases not handled right now.
    Args:
        ref_seq (string): Reference sequence
        i (int): 3' index of del at ref sequence
        length (int): Length of deletion
        num_surBases (int): Number of surrounding bases to consider
    Returns:
        boolean: Whether deletion is ambiguous or not
    """
    orig_del_start = i - length + 1
    orig_sur_start = orig_del_start - num_surBases
    orig_sur_end = i + num_surBases
    orig_sur_seq = ref_seq[orig_sur_start - 1: orig_del_start - 1] + ref_seq[i:orig_sur_end]
    for new_del_end in range(i - length, i + length + 1):  # Alt del end points
        if new_del_end == i:  # Orig end point
            continue
        new_del_start = new_del_end - length + 1
        sur_seq = ref_seq[orig_sur_start - 1: new_del_start - 1] + ref_seq[new_del_end:orig_sur_end]
        if sur_seq == orig_sur_seq:
            return True
    return False


def fix_ref_pos(ref_pos):
    """
    Fills in None values at the soft clipped positions that pysam produces
    """
    if None not in ref_pos:
        return(ref_pos)

    i = 0
    clip_counter = 0
    while ref_pos[i] is None:
        i += 1
        clip_counter += 1
    
    if i > 0:
        for j in range(clip_counter):
            ref_pos[j] = ref_pos[i] - clip_counter + j
   
    if i == len(ref_pos)-1:
        return(ref_pos)
    
    ref_pos = ref_pos[::-1]
    i = 0
    clip_counter = 0
    while ref_pos[i] is None:
        i += 1
        clip_counter += 1

    if i > 0:
        for j in range(clip_counter):
            ref_pos[j] = ref_pos[i] + clip_counter - j
    
    return(ref_pos[::-1])

def complement(dna):
    # this implementation is speedy
    def _complement(x):
        if x == 'A':
            return('T')
        elif x == 'T':
            return('A')
        elif x == 'C':
            return('G')
        elif x == 'G':
            return('C')
        elif x == 'U':
            return('A')
        elif x == 'a':
            return('t')
        elif x == 't':
            return('a')
        elif x == 'c':
            return('g')
        elif x == 'g':
            return('c')
        elif x == 'u':
            return('a')
        elif x == 'N':
            return('N')
    return(''.join([_complement(x) for x in dna]))

def reverse_complement(dna):
    return(complement(dna)[::-1])

def read_pair_generator(bam, reference=None, start=None, end=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    from https://www.biostars.org/p/306041/
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(reference=reference, start=start, end=end):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


# main function, adjust inputs and parameters here!

filename = sys.argv[1]
genome = sys.argv[2]

def run():
	mut_p, cov_p, mut_m, cov_m = coverage_mutations_PE(filename, genome, n_cpu=32, clip_filter=8, return_read_stats=False)

	with open(filename.split('.')[0] + '_mut_cov_PE.pkl', 'wb') as f:
   	    pickle.dump((mut_m, cov_m), f)

if __name__ == '__main__':
	run()
