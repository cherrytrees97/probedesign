"""Consensus sequence of alignment

This module facilitates generating a consensus sequence from a multiple sequence
alignment. The multiple sequence alignment must be in fasta format. 

This module requires that 'Biopython' and 'pandas' is installed within the Python
environment that you are using this module in. 

This module can also be run as a script to generate a consensus sequence.
"""
#Standard libraries
import argparse
import pathlib
#Third-party libraries
from Bio import AlignIO, SeqIO
import numpy as np
import pandas as pd

def parse_args(): 
    parser = argparse.ArgumentParser(__doc__)

    parser.add_argument(
        "target_path",
        action='store', 
        type=pathlib.Path, 
        help='Path to multiple sequence alignment file, in fasta format.'
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_path", 
        action='store',
        default=None, 
        type=pathlib.Path,
        help='Path to directory that the consensus fasta will be output to. '
    )
    parser.add_argument(
        '--min_cons',
        '-c',
        action='store',
        dest='min_cons',
        type=float,
        default=0.9,
        help='Minimum percentage (decimal format) of identical bases for consensus base to be called (0.0-1.0). Default=0.9'
    )
    parser.add_argument(
        '--min_rep',
        '-r', 
        action='store',
        dest='min_rep',
        type=float,
        default=0.5,
        help='Minimum percentage of sequences represented for consensus to be generated (0.0-1.0). Default=0.5',
    )

    args = parser.parse_args()

    #Assign output path to the path of the input file if no output path
    #was specified. 
    if not args.output_path: 
        args.output_path = args.target_path.parent

    #Arg checker
    if not(
        args.output_path.exists()
        and args.target_path.exists()
    ):
        parser.error('Invalid paths.')
    
    return args

class Alignment: 
    """
    Represents a DNA sequence alignment.

    Attributes
    ----------
    alignment : Bio.Align.MultipleSeqAlignment
        The MultipleSeqAlignment object returned by AlignIO.read()
    seq_position_data : pandas.DataFrame
        The sequence alignment represented as a dataframe.
        Representing the alignment as a dataframe simplifies creation of a consensus sequence.
        Each row corrresponds to a position on the alignment. 
        Each column corresponds to an accession.
    seq_regions : dict
        Dictionary containing tuples indexed by accession ID. 
        Tuples contain (start_index, end_index) of first and last base of a sequence 
        in the alignment. 
    consensus : str
        String representing the consensus nucleotide sequence of the alignment. 
    
    Methods
    -------
    get_consensus(threshold: float=0.9) -> None
        Determines the consensus sequence of the alignment. 
    get_accessions() -> list
        Get a list of accession IDs for sequences in the alignment.
    _get_sequence_regions() -> dict
        Get a dictionary with tuples containing start and end indices of the first
        and last base of sequences in the alignment. 
    _get_sequence_position_data() -> pandas.Dataframe
        Get a dataframe representation of the sequence alignment.
    """

    def __init__(self, alignment_path: pathlib.Path): 
        """
        Parameters
        ----------
        alignment_path : pathlib.Path
            The path to the multiple sequence alignment file, in fasta format. 
        """
        self.alignment = AlignIO.read(alignment_path, 'fasta')
        self.sequences = self._get_sequences()
        self.sequence_regions = self._get_sequence_regions()
        self.consensus = None

    def __repr__(self): 
        return self.alignment

    def _get_sequences(self) -> dict:
        """Get ungapped sequences, put into a nice dictionary keyed by accession"""
        sequence_dict = {}
        for sequence in self.alignment: 
            sequence_dict[sequence.id] = str(sequence.seq.ungap())
        return sequence_dict

    def get_consensus(self, min_con: float = 0.9, min_rep: float = 0.5) -> str:
        """
        Get the consensus sequence
        """
        cons_sequence = []
        num_sequence = len(self.alignment)
        
        #Determine sequence representation
        sequence_regions = []
        for sequence in self.alignment: 
            start_base = sequence.seq.ungap()[0]
            end_base = sequence.seq.ungap()[-1]
            start_index = sequence.seq.find(start_base) #inclusive
            end_index = sequence.seq.rfind(end_base)+1 #exclusive
            sequence_regions.append((start_index, end_index))
        
        seq_rep = []

        for position in range(self.alignment.get_alignment_length()): 
            seq_rep.append(0)
            for region in sequence_regions:
                if (position >= region[0]) and (position < region[1]): 
                    seq_rep[position] = seq_rep[position] + 1
        
        seq_rep_ratio = [position/num_sequence for position in seq_rep]

        #Determine which sequences pass the min_rep score
        #Start with the left index
        #Finish with the right index
        left_index = 0
        reverse_right_index = -1
        while (seq_rep_ratio[left_index]) < min_rep: 
            left_index = left_index + 1
        while (seq_rep_ratio[reverse_right_index]) < min_rep: 
            reverse_right_index = reverse_right_index - 1
        right_index = self.alignment.get_alignment_length() + reverse_right_index + 1

        #Determine the dominant base at each position
        for position in range(left_index, right_index): 
            bases = self.alignment[:,position].upper()
            
            num_bases = seq_rep[position]

            counts = {
                'A':bases.count('A'),
                'C':bases.count('C'),
                'T':bases.count('T'),
                'G':bases.count('G'),
            }
            cons_base = max(counts, key=counts.get)
            if (counts[cons_base]/num_bases) >= min_con: 
                cons_sequence.append(cons_base)
            else: 
                cons_sequence.append('N')

        self.consensus = ''.join(cons_sequence)

        #return self.consensus

    def _get_sequence_regions(self) -> dict: 
        """
        Determine the start and end of each sequence in an alignment.

        Goes through each sequence in the alignment, and takes the ungapped sequence.
        The first and last nucleotide of the ungapped sequence is then searched from either
        end of the gapped sequence. 
        The indices are recorded in the sequence_regions dictionary, keyed by the sequence
        accession. 

        Parameters
        ----------
        None

        Return
        ------
        sequence_regions : dict
            Key = sequence accession. Each entry is a tuple of start and end index of the sequence
            in the alignment. 

        """  
        sequence_regions = dict()

        for sequence in self.alignment: 
            ungap_sequence = sequence.seq.ungap()
            start_base = ungap_sequence[0]
            end_base = ungap_sequence[-1]
            start_index = sequence.seq.find(start_base)
            end_index = sequence.seq.rfind(end_base)
            sequence_regions[sequence.id] = (start_index, end_index+1)

        return sequence_regions

def main():
    args = parse_args()

    if args.target_path.is_file(): 
        target_alignment = Alignment(args.target_path)
        target_alignment.get_consensus(args.min_cons, args.min_rep)

        #Write consensus fasta
        consensus_fasta_path = args.output_path.joinpath(f'{args.target_path.stem}_consensus.fasta')
        with open(consensus_fasta_path, 'w') as output_file: 
            output_file.write(f">{args.target_path.stem}\n")
            output_file.write(target_alignment.consensus)

    elif args.target_path.is_dir(): 
        for path in args.target_path.glob('*.fasta'): 
            target_alignment = Alignment(path)
            target_alignment.get_consensus(args.min_cons, args.min_rep)

            #Write consensus fasta
            consensus_fasta_path = args.output_path.joinpath(f'{path.stem}_consensus.fasta')
            with open(consensus_fasta_path, 'w') as output_file: 
                output_file.write(f">{path.stem}\n")
                output_file.write(target_alignment.consensus)

    else: 
        print('Invalid target path provided.')

if __name__ == "__main__":
    main()