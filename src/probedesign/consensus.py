from Bio import AlignIO
import numpy as np
import pandas as pd

class alignment: 
    def __init__(self, alignment_path): 
        self.alignment = AlignIO.read(alignment_path, 'fasta')
        self.seq_position_data = None
        self.sequence_regions = dict()
        self._get_sequence_regions()
        self._get_sequence_position_data()
        self.consensus = None
    def __repr__(self): 
        return self.alignment
    def get_consensus(self, threshold=0.9): 
        """
        
        """
        consensus=[]
        for base_position in self.seq_position_data.iterrows(): 
            nucleotide_counts = base_position[1].value_counts(normalize=True)
            if nucleotide_counts[0] >= threshold: 
                consensus.append(nucleotide_counts.index[0])
            else: 
                consensus.append("n")
        self.consensus = "".join(consensus).replace("-","")
    def _get_sequence_regions(self): 
        """
        Function to determine the start and end of sequences in an alignment
        alignment - Bio.Align.MultipleSeqAlignment object
        sequence region - list of tuples(start_index, end_index)
        0-based indices, half open
        """
        for sequence in self.alignment: 
            ungap_sequence = sequence.seq.ungap()
            start_base = ungap_sequence[0]
            end_base = ungap_sequence[-1]
            start_index = sequence.seq.find(start_base)
            end_index = sequence.seq.rfind(end_base)
            self.sequence_regions[sequence.id] = (start_index, end_index+1)
        return self.sequence_regions
    def _get_sequence_position_data(self):
        """
        Function to store information in a dataframe
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
    def get_accessions(self): 
        list_id = []
        for seq in self.alignment: 
            list_id.append(seq.id.split('.')[0])
        return list_id