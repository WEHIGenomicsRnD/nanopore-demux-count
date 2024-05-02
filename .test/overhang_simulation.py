# Written by Layla Wang, WEHI Genomics Lab

import numpy as np
from random import choice, randint

# Constants for the simulation
NUM_READS = 1000
NUCLEOTIDES = ['A', 'T', 'C', 'G']
QUALITY_SCORE = 'I'  # Placeholder for a high-quality score

# Primer Sequences
FORWARD_PRIMER = 'GGGGGATAACATTGAACTTC'
REVERSE_PRIMER = 'CCTAATATACGGACGCAATC'

# Index Sequences and Identifiers
FWD_INDEX_1 = 'TAGATCGC'
FWD_INDEX_2 = 'CTCTCTAT'
REV_INDEX_1 = 'TCGCCTTA'
REV_INDEX_2 = 'CTAGTACG'
IDENTIFIER_FWD_1 = 'TACG'
IDENTIFIER_FWD_2 = 'GTAC'
IDENTIFIER_REV = 'TCAG'

# Full Index Combinations
INDEX_COMBINATIONS = [
    (IDENTIFIER_FWD_1 + FWD_INDEX_1, REV_INDEX_1 + IDENTIFIER_REV),
    (IDENTIFIER_FWD_1 + FWD_INDEX_2, REV_INDEX_1 + IDENTIFIER_REV),
    (IDENTIFIER_FWD_1 + FWD_INDEX_1, REV_INDEX_2 + IDENTIFIER_REV),
    (IDENTIFIER_FWD_1 + FWD_INDEX_2, REV_INDEX_2 + IDENTIFIER_REV),
    (IDENTIFIER_FWD_2 + FWD_INDEX_1, REV_INDEX_1 + IDENTIFIER_REV),
    (IDENTIFIER_FWD_2 + FWD_INDEX_2, REV_INDEX_1 + IDENTIFIER_REV),
    (IDENTIFIER_FWD_2 + FWD_INDEX_1, REV_INDEX_2 + IDENTIFIER_REV),
    (IDENTIFIER_FWD_2 + FWD_INDEX_2, REV_INDEX_2 + IDENTIFIER_REV)
]

# Guide parameters
NUM_GUIDES = 20
GUIDE_LENGTH = 1000


# Functions for introducing variations
def introduce_mismatch(sequence):
    if sequence:
        position = randint(0, len(sequence) - 1)
        original_base = sequence[position]
        bases = {'A', 'T', 'C', 'G'}
        bases.remove(original_base)
        sequence = list(sequence)
        sequence[position] = choice(list(bases))
    return ''.join(sequence)


def introduce_deletion(sequence):
    if sequence:
        position = randint(0, len(sequence) - 1)
        sequence = sequence[:position] + sequence[position+1:]
    return sequence


def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))


# Function to generate a single read with variations
def generate_varied_read(read_id, combination_index, guide_sequences):
    forward_index, reverse_index = INDEX_COMBINATIONS[combination_index]

    # Apply variations to primers
    fwd_primer = FORWARD_PRIMER if np.random.rand() >= 0.05 else ''
    rev_primer = REVERSE_PRIMER if np.random.rand() >= 0.05 else ''
    fwd_primer = introduce_mismatch(fwd_primer) if np.random.rand() < 0.05 else fwd_primer
    fwd_primer = introduce_deletion(fwd_primer) if np.random.rand() < 0.05 else fwd_primer

    # Create the random sequences
    start_random_sequence = ''.join(choice(NUCLEOTIDES) for _ in range(randint(0, 100)))
    end_random_sequence = ''.join(choice(NUCLEOTIDES) for _ in range(randint(0, 100)))

    # pick a random guide sequence from the list
    guide_seq = guide_sequences[randint(0, len(guide_sequences) - 1)]

    # Assemble the read
    sequence = f"{start_random_sequence}{forward_index}{fwd_primer}{guide_seq}{rev_primer}{reverse_index}{end_random_sequence}"

    # Reverse complement 50% of the reads
    if read_id % 2 == 0:
        sequence = reverse_complement(sequence)

    quality = QUALITY_SCORE * len(sequence)
    read_header = f"@varied_read_{read_id}"

    return f"{read_header}\n{sequence}\n+\n{quality}\n"


# Generate guide sequence fasta
guide_sequences = []
guides_filename = "guides_simulated.fasta"
with open(guides_filename, 'w') as guides_file:
    for guide_id in range(1, NUM_GUIDES + 1):
        guide_sequence = ''.join(choice(NUCLEOTIDES) for _ in range(GUIDE_LENGTH))
        guide_sequences.append(guide_sequence)
        guides_file.write(f">guide_{guide_id}\n{guide_sequence}\n")

# Write the FastQ file
fastq_filename = "overhang_simulated_reads.fastq"
with open(fastq_filename, 'w') as fastq_file:
    for read_id in range(1, NUM_READS + 1):
        combination_index = (read_id - 1) % len(INDEX_COMBINATIONS)
        fastq_file.write(generate_varied_read(read_id, combination_index, guide_sequences))

print(f"FastQ file generated: {fastq_filename}")
