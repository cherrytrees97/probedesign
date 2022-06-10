from Bio import AlignIO
import numpy as np
import pandas as pd
import argparse
import pathlib

def parse_args(): 
    parser = argparse.ArgumentParser("Generate a consensus sequence.")
    parser.add_argument(
        "target_path",
        action='store', 
        type=pathlib.Path, 
        help='Path to alignment'
    )
    parser.add_argument(
        "--output",
        dest="output_path", 
        action='store',
        default=None, 
        type=pathlib.Path,
        help='Output path'
    )
    args = parser.parse_args()
    if not args.output_path: 
        args.output_path = args.target_path.parent
    return (args.target_path, args.output_path)

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

    def __init__(self, alignment_path): 
        """
        Parameters
        ----------
        alignment_path : pathlib.Path
            The path to the multiple sequence alignment file, in fasta format. 
        """
        self.alignment = AlignIO.read(alignment_path, 'fasta')
        self.seq_position_data = None
        self.sequence_regions = dict()
        self._get_sequence_regions()
        self._get_sequence_position_data()
        self.consensus = None
    def __repr__(self): 
        return self.alignment
    def get_consensus(self, threshold: float=0.9) -> None: 
        """
        Get the consensus sequence of the alignment. 

        Uses seq_position_data (dataframe representation) of the alignment to generate
        the consensus sequence. 
        TODO: Figure out what the algorithm is and describe it here. 

        Parameters
        ----------
        threshold : float=0.9
            Base consensus is only called if the most frequent basecall exceeds the
            threshold percentage. 
        """

        consensus=[]
        for base_position in self.seq_position_data.iterrows(): 
            nucleotide_counts = base_position[1].value_counts(normalize=True)
            if nucleotide_counts[0] >= threshold: 
                consensus.append(nucleotide_counts.index[0])
            else: 
                consensus.append("n")
        self.consensus = "".join(consensus).replace("-","")
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
        """  

        for sequence in self.alignment: 
            ungap_sequence = sequence.seq.ungap()
            start_base = ungap_sequence[0]
            end_base = ungap_sequence[-1]
            start_index = sequence.seq.find(start_base)
            end_index = sequence.seq.rfind(end_base)
            self.sequence_regions[sequence.id] = (start_index, end_index+1)
        return self.sequence_regions
    def _get_sequence_position_data(self) -> pd.DataFrame:
        """
        Converts MultipleSeqAlignment object into a dataframe representation. 

        The dataframe representation of the alignment is easier to use, especially when 
        implementing the threshold calculation method of the consensus sequence. 

        Parameters
        ----------
        None
        """
        dict_series = dict()
        for sequence in self.alignment: 
            #Create the series containing all of the data
            series = pd.Series(list(sequence.seq), index=range(len(sequence.seq)), name=sequence.id)
            #Using sequence_regions, convert all non-sequence gap characters to NaN
            sequence_region = self.sequence_regions[sequence.id]
            #Start of the alignment to first base of sequence
            #Base past end of sequence to end of alignment
            series.iloc[0:sequence_region[0]] = np.NaN
            series.iloc[sequence_region[1]:len(sequence.seq)] = np.NaN
            dict_series[sequence.id] = series
        self.seq_position_data = pd.DataFrame(dict_series)
        return self.seq_position_data
    def get_accessions(self) -> list: 
        """Get the accessions for the sequences in the alignment. """
        list_id = []
        for seq in self.alignment: 
            list_id.append(seq.id.split('.')[0])
        return list_id

if __name__ == "__main__":
    target_path, output_path = parse_args()
    target_alignment = Alignment(target_path)
    target_alignment.get_consensus()
    output_file = open(output_path.joinpath(f'{target_path.stem}_consensus.fasta'), 'w')
    output_file.write(f">{target_path.stem}\n")
    output_file.write(target_alignment.consensus)
    output_file.close()
