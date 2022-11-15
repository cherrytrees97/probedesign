"""

"""
#Standard libraries
import csv
import tempfile
import multiprocessing
import subprocess
import io
import pathlib
from math import floor
#Third-party libraries
from Bio import SeqIO, Seq
import pandas as pd
from primer3 import calcTm
#Local modules
from consensus import Alignment   
from tmcalc import CalcProbeTm

class Oligo: 
    """
    Oligo - representation of a oligonucleotide. Either primer or probe. 

    Attributes
    ----------
    seq : str 
        oligonucleotide sequence
    root_pos : int 
        position of first nucleotide on the template sequence used to generate the oligo.
        0-based.
    len : int
        length of the oligonucleotide sequence
    tm : float 
        calculated melting temperature of the oligo
    id : str 
        sequence ID 
    sensitivity : float 
        0-1; proportion of target sequences that contain the oligonucleotide sequence
    specificity : float 
        0-1; proportion of database sequences that do not match oligonucleotide sequence
    score : float 
        0-2; sum of the sensitivity and specificity score
    target_accessions : list 
        list of strings representing target accession IDs
    amplified_accessions : list 
        list of strings representing target accessions this oligo binds to

    Methods
    -------
    calculate_sensitivity : calculate the sensitivity of the oligo on the target accessions.
    calculate_specificity : calculate the specificity of the oligo using BLAST results against.
                            database
    calculate_score       : return the sum of the sensitivity and specificity value.
    """

    def __init__(self, root_pos: int, seq: str, tm: float):
        """
        """
        #Input parameters
        self.seq = seq
        self.root_pos = root_pos
        self.len = len(self.seq)
        self.tm = tm
        #Calculated
        self.id =  f"{str(self.root_pos)}-{str(self.len)}"
        self.amplified_accessions = []
        #Defined after calling methods
        self.sensitivity = 0
        self.specificity = 0
        self.score = 0
        self.target_accessions = []

    def _calculate_target_accessions(self, seq_regions: dict):
        """
        Determines which target accessions contributed to the creation of the oligo. This 
        is defined as target accessions that have start and end coordinates that contain the
        oligo. 

        This method is only called if sensitivty/specificty check is being done on the oligo. 

        Parameters
        ----------
        seq_regions : dict
            List of tuples containing start and end coordinates of the sequence in alignment.
            Refer to the Alignment class in consensus.py for more information on how
            this is generated. 

        Return
        ------
        target_accessions : list
            List of accession IDs stored as strings
        """

        #Determine the range the probe spans
        oligo_start = self.root_pos - 1
        oligo_end = oligo_start + self.len

        #Determine target_accessions.
        #Target accessions are defined as accessions where it contributes 
        #to the alignment that the oligo sequence is derived from. 
        target_accessions = []
        for target in seq_regions: 
            if (
                oligo_start >= seq_regions[target][0]
                and oligo_end <= seq_regions[target][1]
            ):
                target_accessions.append(target)
        
        self.target_accessions = target_accessions

        return target_accessions

    def calculate_sensitivity(self, alignment: Alignment, reverse: bool = False) -> float:
        """
        Calculates sensitivity using solely the alignment. Essentially, a string search
        is done for the oligo sequence on each sequence in the alignment

        Parameters
        ----------
        alignment : Alignment
            Alignment object. Refer to Alignment class in consensus.py 
        
        Return
        ------
        sensitivity : float
            Sensitivity score. Ranges between 0 to 1. 1 is most sensitive. 
        """

        #Identify target accessions if it is not there
        amplified_accessions = []

        if not self.target_accessions: 
            seq_regions = alignment.sequence_regions
            self._calculate_target_accessions(seq_regions)

        #TODO: MAKE THIS NOT GARBAGE BUT LOGICALLY THIS SHOULD MAKE SENSE
        i = 0
        for accession in self.target_accessions: 
            if not reverse:
                sequence = alignment.sequences[accession].upper()
            else: 
                sequence = str(Seq.Seq(alignment.sequences[accession]).reverse_complement()).upper()

            if self.seq in sequence: 
                i += 1
                amplified_accessions.append(accession)

        sensitivity = i/len(self.target_accessions)
        self.amplified_accessions = amplified_accessions
        self.sensitivity = sensitivity

        return sensitivity

    def calculate_specificity(self, alignment: Alignment, blast_results: pd.DataFrame, blastdb_len: int): 
        """
        Calculates specificity of the oligo (binding to non-target sequences). 
        Binding is defined as simply appearing in the BLAST results --> this will be a overestimation of the specificity,
        making it the 'worst case' scenario. 
        Specificity is calculated by the following formula: 
            TN / (TN + FP)
            where:
            TN = non-target accessions that were not amplified
            TN = (total non-target) - (amplified non-target)
            FP = amplified non-target
            resulting in the following formula: 
            ((Total non-target) - (amplified non-target)) / total non-target

        Parameters
        ----------
        alignment : Alignment
            Alignment object. Refer to Alignment class in consensus.py 
        blast_results : DataFrame
            Blast results
        blastdb_len : int
            Length of the BLAST database. 

        Return
        ------
        specificity : float
            Specificity score. Ranges between 0 to 1. 1 is most specific. 
        """

        seq_regions = alignment.sequence_regions
        #Identify target accessions if it is not there
        if not self.target_accessions: 
            self._calculate_target_accessions(seq_regions)

        #Get all BLAST accession hits, even partial hits
        blast_match_accessions = set(blast_results.loc[:,'sacc'])

        #Remove every target_accession from the blast_match_accessions list.
        #Store the number of every target_accession removed. 
        num_removed_target_accession = 0
        for accession in self.target_accessions:
            if accession in blast_match_accessions:  
                blast_match_accessions.remove(accession)
                num_removed_target_accession += 1

        #Calculate specificity
        #Total non-target = all_blast_sequences - target_accessions
        total_non_target = blastdb_len - num_removed_target_accession
        specificity = (blastdb_len - len(blast_match_accessions))/total_non_target

        self.specificity = specificity

        return specificity

    def calculate_score(self): 
        """Return sum of sensitivity and specificity"""
        self.score = self.sensitivity + self.specificity
        return self.score

class PrimerPair: 
    """
    PrimerPair - representation of a primer pair. 

    Attributes
    ----------
    fw_primer : Oligo
        fw primer object
    rev_primer : Oligo
        rev primer object
    sensitivity : float
        0-1; proportion of target sequences that contain the oligonucleotide sequence
    specificity : float
        0-1; proportion of database sequences that do not match oligonucleotide sequence
    score : float
        0-2; sum of sensitivity and specificity scores. 
    target_accesions : list 
        list of strings representing target accession IDs
    amplified_accessions : list
        list of strings representing amplified accessions

    Methods
    -------
    _calculate_combined_target_accessions : 
    calc_tm_diff : 
    calculate_sensitivity : calculate the sensitivity of the oligo on the target accessions.
    calculate_specificity : calculate the specificity of the oligo using BLAST results against.
                            database
    calculate_score       : return the sum of the sensitivity and specificity value.
    """

    def __init__(self, fw_primer, rev_primer): 
        self.fw_primer = fw_primer
        self.rev_primer = rev_primer
        self.sensitivity = 0
        self.specificity = 0
        self.score = 0
        self.target_accessions = self._calculate_combined_target_accessions()
        self.amplified_accessions = None

    def _calculate_combined_target_accessions(self):
        """Return list of accessions where they're found in both fw and rev primer"""

        target_accessions = []

        for accession in self.fw_primer.target_accessions: 
           if accession in self.rev_primer.target_accessions: 
                target_accessions.append(accession)

        return target_accessions

    def calc_tm_diff(self):
        """Calculate the difference between the forward and reverse primer tm """ 
        return abs(self.fw_primer.tm - self.rev_primer.tm)

    def calculate_sensitivity(self): 
        """
        Function calculates the sensitivity of a primer pair given the BLAST results for both. 
        Amplification is defined as an accession where both the forward and reverse primers have 
        100% query coverage and 100% percent identity to the target sequence. 
        Sensitivity is calculated by the following formula: 
            TP / (TP + FN)
            where: 
            TP = succesfully amplified accessions
            FN = target accessions that were not amplified
            TP + FN = total number of target accessions

        Parameters
        ----------
        self.fw_primer.amplified_accessions : list
        self.rev_primer.amplified_accessions : list
        self.target_accessions : list

        Return
        ------
        self.amplified_accessions : list
        sensitivity : float
            Sensitivity score. Ranges from 0 to 1. 
        """

        fw_primer_amplified = self.fw_primer.amplified_accessions
        rev_primer_amplified = self.rev_primer.amplified_accessions
        target_accessions = self.target_accessions

        num_amplified = 0
        amplified_accessions = []

        for accession_fw in fw_primer_amplified: 
            if accession_fw in rev_primer_amplified: 
                num_amplified += 1
                amplified_accessions.append(accession_fw)
        
        sensitivity = num_amplified/len(target_accessions)
        self.amplified_accessions = amplified_accessions

        self.sensitivity = sensitivity

        return sensitivity

    def calculate_specificity(self, fw_blast_results, rev_blast_results, blastdb_len):
        """
        Calculates specificity of a primer pair using BLAST results. 

        Parameters
        ----------
        fw_blast_results : DataFrame
            Blast results for the forward primer.
        rev_blast_results : DataFrame
            Blast results for the reverse primer.
        blastdb_len : int
            Length of the BLAST database used to generate the BLAST results.

        Return
        ------
        specificity : float
            Specificity score. Ranges between 0 to 1. 1 is most specific. 
        """

        #Take the set of the accessions
        fw_match_accessions = set(fw_blast_results.loc[:,'sacc'])
        rev_match_accessions = set(rev_blast_results.loc[:,'sacc'])
        amplified_accessions = set(fw_match_accessions&rev_match_accessions)
        num_removed_target_accession = 0

        #Remove every target_accession from the blast_match_accessions list
        #Store the number of every target_accession removed. 
        for accession in self.amplified_accessions:
            if accession in amplified_accessions:  
                amplified_accessions.remove(accession)
                num_removed_target_accession += 1

        #Calculate specificity
        #Total non-target = all_blast_sequences - target_accessions
        total_non_target = blastdb_len - num_removed_target_accession
        specificity = (blastdb_len - len(amplified_accessions))/total_non_target

        self.specificity = specificity

        return specificity

    def calculate_score(self):
        """Return sum of sensitivity and specificity""" 
        self.score = self.sensitivity + self.specificity

class ProbeGenerator: 
    """
    Class to handle the generation of probes from a target sequence.

    Attributes
    ----------
    template : str
        DNA sequence that is used as the template for the generation of probes
    start : int
        The start index on the template sequence that is used as the first root position
        for the generation of probes. 0-based. 
    end : int
        The end index on the template sequence that defines the sequence used for probe design.
        0-based. 
    min_length : int
        Minimum length of any generated probe. 
    max_length : int
        Maximum length of any generated probe. 
    probes : list
        List of probes that are legal. 

    Methods
    -------
    get_probes : finds all viable probes in given template sequence in given range. 
    output : output all probes and their associated information as a .csv file. 
    """

    def __init__(
        self, 
        template_sequence, 
        start, 
        end, 
        min_length, 
        max_length, 
        ):
        #TODO: figure out if comment below is true..? I don't think it is.
        #Note that the coordinates are 1-based, so the template start has to be -1
        #Set template sequence to uppercase. 
        self.template_sequence = template_sequence[start:end].upper()
        self.start = start
        self.end = end
        self.min_length = min_length
        self.max_length = max_length
        self.probes = []

    def get_probes(self) -> list: 
        """
        Function finds all viable probes. 
        Probe criteria: 
        1) 5' end of the probe is not a G
        2) More C's than G's
        3) No more than 4 nucleotides in a row
        4) 30 - 80% CG content

        Parameters
        ----------
        self.template_sequence : str
            DNA sequence that is used as the template for the generation of probes.
        self.start : int
            The start index on the template sequence that is used as the first root position
            for the generation of probes. 0-based. 
        self.end : int
            The end index on the template sequence that defines the sequence used for probe design.
            0-based. 
        self.min_length : int
            Minimum length of any generated probe. 
        self.max_length : int
            Maximum length of any generated probe. 

        Return
        ------
        probes : list
            List of Oligo objects, where each Oligo object is a viable probe. 
            Note that the Oligo object's indices are stored as 1-based, fully closed coordinates. 
        """

        def check_probe(input_seq): 
            """Check probe sequence for probe conditions."""

            def count_c_g(seq):
                """Count the number of CGs in the sequence."""
                return seq.count('C') + seq.count('G')

            def chk_run(seq): 
                """
                Ensure that there are no runs of four or more identical nucleotides 
                in the probe.
                """

                #Store current nucleotide identity
                current_nucleotide = seq[0]

                #Length of the run of current nucleotide
                len_run = 1
                detected = False

                #Go through each nucleotide in the primer.
                #Increment len_run if nucleotide position matches current_nucleotide.
                #If nucleotide changes, reset len_run to 1 and change current_nucleotide. 
                #If len_run exceeds 4, terminate
                for nucl in seq: 
                    if nucl == current_nucleotide: 
                        len_run = len_run + 1
                        if len_run > 3: 
                            detected = True
                            break
                    else: 
                        current_nucleotide = nucl
                        len_run = 1

                return detected

            def chk_last_5(seq): 
                """Check if the last five nucleotides have more than 2 C or Gs"""
                return count_c_g(seq[-5:]) > 2
            
            seq = input_seq.upper()

            #Check the three conditions
            percent_c_g = count_c_g(seq)/len(seq)
            run_flag = chk_run(seq)
            last5_flag = chk_last_5(seq)
            if(
                0.3 <= percent_c_g <= 0.8 
                and run_flag is False
                and last5_flag is False
                and seq[0] != 'G'
                and not('N' in seq)
            ):
                return True
            else: 
                return False
        
        #Trim template to specified region
        template = self.template_sequence[self.start:self.end]

        probes = []

        #For each root position in the target sequence
        for i in range(len(template) - self.min_length + 1):
            #Iterate over all probe lengths
            for probe_len in range(self.min_length, self.max_length + 1): 
                probe_seq = template[i:i+probe_len]
                if check_probe(probe_seq) is True: 
                    #Note that the coordinates are converted back to 1-based
                    probe_tm = CalcProbeTm(probe_seq).Tm
                    probes.append(Oligo(self.start+i+1, probe_seq, probe_tm))

        self.probes = probes

        return probes

    def output(
        self, 
        path: pathlib.Path, 
        f_seq_rep: bool, 
        min_seq_rep: float,
        num_seq: int,
        f_tm: bool,
        min_tm: float,
        max_tm: float,
        ): 
        """
        Output a .csv file containing all of the probes and associated information.

        Parameters
        ----------
        path : pathlib.Path
            Path to output directory. 
        
        Returns
        -------
        None

        """

        probe_data = []

        #Create list of tuples from list of Oligo objects
        for probe in self.probes: 
            if probe.target_accessions: 
                accessions_covered = probe.target_accessions
            else:
                accessions_covered = None
            
            #seq rep filter
            if f_seq_rep and ((len(accessions_covered)/num_seq) < min_seq_rep):
                continue
            #tm filter
            if f_tm and not(min_tm < probe.tm <= max_tm):
                continue

            probe_data.append(
                (
                    probe.root_pos,
                    probe.len, 
                    probe.seq,
                    len(accessions_covered),
                    probe.tm, 
                    probe.sensitivity, 
                    probe.specificity,
                    probe.score,
                )
            )

        probe_data.sort(key=lambda x: x[7], reverse=True)

        #Output csv
        output_path = path.joinpath('probe_candidates.csv')
        with open(output_path, 'w', newline='') as csv_file: 
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(
                (
                    'probe_root',
                    'probe_len',
                    'probe_seq',
                    'seq_rep',
                    'tm',
                    'sens',
                    'spec',
                    'score',
                )
            )
            csv_writer.writerows(probe_data)

class PrimerGenerator: 
    """
    Class to handle the generation of primers from a target sequence, targetting a specific probe. 

    Attributes
    ----------
    template : str
        DNA sequence that is used as the template for the generation of primers.
    pb_start : int
        The start index of the probe on the template sequence to design primers for. 
    pb_len : int
        The length of the probe to design primers for.  
    min_length : int
        Minimum length of any generated primer. 
    max_length : int
        Maximum length of any generated primer. 
    max_tm_diff : int
        Maximum annealing temperature difference allowed between any primer pair. 
    fw_primers : list
        List of fw primers that are legal. 
    rev_primers : list
        List of rev primers that are legal. 

    Methods
    -------
    get_probes : finds all viable probes in given template sequence in given range. 
    output : output all probes and their associated information as a .csv file. 
    """

    def __init__(
        self,
        template, 
        pb_start, 
        pb_len,
        min_length, 
        max_length, 
        max_tm_diff
        ):
        #Input
        #Set template sequence to uppercase. 
        self.template = template.upper()
        self.pb_start = pb_start
        self.pb_len = pb_len
        self.min_length = min_length
        self.max_length = max_length
        self.max_tm_diff = max_tm_diff
        #Calculated
        self.pb_end = pb_start-1+pb_len
        #Output
        self.fw_primers = []
        self.rev_primers = []
        self.primer_pairs = []

    def find_fw_primers(self):
        """
        Function finds all viable FW primers.
        FW primer criteria: 
        1) 3'-end of the FW needs to be within 50 bp of the 5'-end of the probe.
        2) Between min and max length 
        3) % GC is 30% to 80%
        4) Last five nucleotides at the 3' end contain no more than two G + C residues
        5) No more than 4 consecutive nucleotides within the primer 

        Input data: 
        1) target_seq - str - target sequence - ATCTGATCATGATCATGACTAGTCATGGC
        2) pb_start - int - start index of 5'-end of the probe - 607

        Output data: 
        [
            {root_pos:0, len:17, seq:''}
        ]

        Algorithm: 
        1st root position is the nucleotide right before the 5'-end of the probe. 
        1st_root_pos_index = pb_start - 1
        However, for list-slicing purposes, take the pb_start index. 
        For each root position, slice the sequence from the root position to fw_length for each fw_length allowed. 
        Ex. 
            MIN_LENGTH = 17
            MAX_LENGTH = 22
            pb_start = 607
            1st_root_position = pb_start - 1
            Given the above: 
            primer1 = [590:607] #17 bp primer
            primer2 = [589:607]
            primer3 = [588:607]
            primer4 = [587:607]
            primer5 = [586:607]
            primer6 = [585:607] #22 bp primer

            fw_end --> pb_start - i
            fw_start --> fw_end - fw_len

            {root_pos:pb_start - 1 - i, len:fw_len, seq:[pb_start - i - fw_len :pb_start - i]}
        
        Parameters
        ----------
        self.template : str
            DNA sequence that is used as the template for the generation of probes.
        self.pb_start : int
            The start index of the probe. 0-based. 
        self.min_length : int
            Minimum length of any generated primer. 
        self.max_length : int
            Maximum length of any generated primer.
        
        Return
        ------
        fw_primers : list
            List of Oligo objects, where each Oligo object is a viable fw primer. 
            Note that the Oligo object's indices are stored as 1-based, fully-closed coordinates. 
        """

        def check_primer(input_seq): 
            """
            Check primer sequence for primer conditions.
            
            1) Between specified min and max primer length.
            2) % GC is 30% to 80%.
            3) Last five nucleotides at the 3' end contain no more than two G + C residues.
            4) No more than 4 consecutive nucleotides within the primer.

            Parameters
            ----------
            seq : str
                Sequence of the primer. 

            Return
            ------
            True if primer passes, false if primer fails. 
            """
            
            def count_c_g(seq): 
                """Count the number of CGs in the sequence."""
                return seq.count('C') + seq.count('G')

            def chk_run(seq): 
                """
                Ensure that there are no runs of four or more identical nucleotides 
                in the probe.
                """

                #Store current nucleotide identity
                current_nucleotide = seq[0]

                #Length of the run of current nucleotide
                len_run = 1
                detected = False

                #Go through each nucleotide in the primer.
                #Increment len_run if nucleotide position matches current_nucleotide.
                #If nucleotide changes, reset len_run to 1 and change current_nucleotide. 
                #If len_run exceeds 4, terminate
                for nucl in seq: 
                    if nucl == current_nucleotide: 
                        len_run = len_run + 1
                        if len_run > 3: 
                            detected = True
                            break
                    else: 
                        current_nucleotide = nucl
                        len_run = 1

                return detected

            def chk_last_5(seq):
                """Count the number of CGs in the sequence.""" 
                return count_c_g(seq[-5:]) > 2
            
            seq = input_seq.upper()
            
            #Check the three conditions
            percent_c_g = count_c_g(seq)/len(seq)
            run_flag = chk_run(seq)
            last5_flag = chk_last_5(seq)

            #Return true if all of the conditions are passed, otherwise return false
            if(
                0.3 <= percent_c_g <= 0.8 
                and run_flag is False
                and last5_flag is False
                and not('N' in seq)
            ):
                return True
            else: 
                return False
        
        fw_primer_indices = []
        fw_primers = []

        #Determine if from the current root position of the probe, all of the indices before are valid. 
        #Otherwise, use a truncated range to search for forward primers
        if self.pb_start - self.max_length - 50 < 0: 
            range_len = self.pb_start - self.max_length + 1
            print("Start too close to beginning of seq. Truncated cases")

            #This range length only applies longest case; there's another few cases before it
            #Calculate edge cases
            edge_range = 50 - (self.pb_start - self.max_length) + 1

            for edge_index in range(1, edge_range): 
                length = self.max_length - edge_index
                if length >= self.min_length: 
                    fw_primer_indices.append((0, length))

        else: 
            range_len = 51

        for gap_size in range(range_len): 
            for fw_len in range(self.min_length, self.max_length + 1):
                fw_start = self.pb_start - gap_size - fw_len
                fw_end = self.pb_start - gap_size
                fw_primer_indices.append((fw_start, fw_end))

        #Get the sequence and check
        for index in fw_primer_indices:
            sequence = self.template[index[0]:index[1]]
            if check_primer(sequence): 
                fw_primers.append(
                    Oligo(
                        index[0] + 1,
                        sequence,
                        float(calcTm(sequence, dv_conc=1.5)),
                    )
                )

        self.fw_primers = fw_primers
        
        return fw_primers

    def find_rev_primers(self):
        """
        Function finds all viable REV primers.

        Parameters
        ----------
        self.template : str
            DNA sequence that is used as the template for the generation of reverse primers.
        self.pb_len : int
            Length of the probe the primers are being designed for.
        self.pb_end : int
            Index of 3' end of the probe that the primers are being designed for. 0-based.
        min_length : int
            Minimum length of any generated primer. 
        max_length : int
            Maximum length of any generated primer. 

        Return
        ------
        rev_primers : list
            List of Oligo objects, where each Oligo object is a viable reverse primer. 
            Note that the Oligo object's indices are stored as 1-based, fully closed coordinates.
        """
        
        def check_primer(input_seq): 
            """
            Check primer sequence for primer conditions.
            
            1) Between specified min and max primer length.
            2) % GC is 30% to 80%.
            3) Last five nucleotides at the 3' end contain no more than two G + C residues.
            4) No more than 4 consecutive nucleotides within the primer.

            Parameters
            ----------
            seq : str
                Sequence of the primer. 

            Return
            ------
            True if primer passes, false if primer fails. 
            """
            
            def count_c_g(seq): 
                """Count the number of CGs in the sequence."""
                return seq.count('C') + seq.count('G')

            def chk_run(seq): 
                """
                Ensure that there are no runs of four or more identical nucleotides 
                in the probe.
                """

                #Store current nucleotide identity
                current_nucleotide = seq[0]

                #Length of the run of current nucleotide
                len_run = 1
                detected = False

                #Go through each nucleotide in the primer.
                #Increment len_run if nucleotide position matches current_nucleotide.
                #If nucleotide changes, reset len_run to 1 and change current_nucleotide. 
                #If len_run exceeds 4, terminate
                for nucl in seq: 
                    if nucl == current_nucleotide: 
                        len_run = len_run + 1
                        if len_run > 3: 
                            detected = True
                            break
                    else: 
                        current_nucleotide = nucl
                        len_run = 1

                return detected

            def chk_last_5(seq):
                """Count the number of CGs in the sequence.""" 
                return count_c_g(seq[-5:]) > 2
            
            seq = input_seq.upper()
            
            #Check the three conditions
            percent_c_g = count_c_g(seq)/len(seq)
            run_flag = chk_run(seq)
            last5_flag = chk_last_5(seq)

            #Return true if all of the conditions are passed, otherwise return false
            if(
                0.3 <= percent_c_g <= 0.8 
                and run_flag is False
                and last5_flag is False
                and not('N' in seq)
            ):
                return True
            else: 
                return False

        template = str(self.template)

        rev_primer_indices = []
        rev_primer_indices_real = []
        rev_primers = []

        #last_root_pos --> assumes fw primer of min length directly adjacent to probe
        #indicates boundary of last viable rev primer root given the 150 bp limit
        last_root_pos = 150 - 2*self.min_length - self.pb_len

        #Iterate through all of the root positions and get primer indices
        for i in range(last_root_pos + 1):
            for rev_len in range(self.min_length, self.max_length + 1):
                rev_start = self.pb_end + i
                rev_end = rev_start + rev_len
                rev_primer_indices.append((rev_start, rev_end))

        #Remove any indices where they are past end of the sequence
        for index in rev_primer_indices: 
            if index[1] <= len(template): 
                rev_primer_indices_real.append(index)

        #Get rev primer sequence and check for conditions
        for index in rev_primer_indices_real: 
            rev_primer_seq = str(Seq.Seq(template[index[0]:index[1]]).reverse_complement())
            if check_primer(rev_primer_seq) is True: 
                rev_primers.append(
                    Oligo(
                        index[0], 
                        rev_primer_seq, 
                        float(calcTm(rev_primer_seq, dv_conc=1.5))
                    )
                )
        
        self.rev_primers = rev_primers

        return rev_primers

    def find_primer_pairs(self): 
        """
        Function will pair the fw and rev primers together. 
        The function iterates over the list of FW primers, and pairs them with all viable reverse primers. 
        Viable reverse primers --> a reverse primer such that the amplicon length is equal to or less than 
        150 bp. 

        Algorithm: 
        Amplicon length = (rev_root_pos + rev_len) - fw_root_pos
        Root_pos_limit = FW_root_pos + 150 - MIN_PRIMER_LEN

        Parameters
        ----------
        self.fw_primers : list
            List of all viable fw primers. 
        self.rev_primers : list
            List of all viable rev primers. 
        Return
        ------

        """
        primer_pairs = []

        for fw_primer in self.fw_primers: 
            
            root_pos_limit = fw_primer.root_pos + 150 - self.min_length
            
            #Check two conditions:
            #1) Amplicon size limit
            #2) Max tm difference
            #Using continue instead of break just in case the list is or some reason
            #not sorted by its root position.
            for rev_primer in self.rev_primers:
                if (
                    rev_primer.root_pos < root_pos_limit
                    and (rev_primer.root_pos + rev_primer.len - fw_primer.root_pos) <= 150
                ):
                    primer_pairs.append(
                        PrimerPair(
                            fw_primer, 
                            rev_primer
                        )
                    )
                else: 
                    continue
        
        self.primer_pairs = primer_pairs

        return primer_pairs

    def output(
        self,
        path: pathlib.Path,
        f_tm: bool,
        max_tm_diff: float,
        max_tm: float,
        min_tm: float,
        f_seq_rep: bool,
        min_seq_rep: float,
        num_seq: int,
        ): 
        """
        Output a .csv file containing all of the primer pairs and associated information.

        Parameters
        ----------
        path : pathlib.Path
            Path to output directory. 
        
        Returns
        -------
        None

        """

        primer_data = []

        #Create a list of tuples from list of PrimerPair objects. 
        for primer_pair in self.primer_pairs:
            #Apply temperature filter, if specified 
            tm_diff = abs(primer_pair.fw_primer.tm - primer_pair.rev_primer.tm)
            if (
                f_tm
                and
                (
                    tm_diff > max_tm_diff
                    or primer_pair.fw_primer.tm > max_tm
                    or primer_pair.fw_primer.tm < min_tm
                    or primer_pair.rev_primer.tm > max_tm
                    or primer_pair.rev_primer.tm < min_tm
                )
            ): 
                continue
            
            #Apply sequence representation filter, if specified
            p_fw_accession = len(primer_pair.fw_primer.target_accessions)/num_seq
            p_rev_accession = len(primer_pair.rev_primer.target_accessions)/num_seq
            if (
                f_seq_rep
                and (
                    p_fw_accession < min_seq_rep
                    or p_rev_accession < min_seq_rep
                )
            ):
                continue
            
            primer_data.append(
                (
                    primer_pair.fw_primer.root_pos,
                    primer_pair.fw_primer.len,
                    primer_pair.fw_primer.seq,
                    len(primer_pair.fw_primer.target_accessions),
                    primer_pair.fw_primer.tm,
                    primer_pair.fw_primer.sensitivity,
                    primer_pair.fw_primer.specificity,
                    primer_pair.rev_primer.root_pos,
                    primer_pair.rev_primer.len,
                    primer_pair.rev_primer.seq,
                    len(primer_pair.rev_primer.target_accessions), 
                    primer_pair.rev_primer.tm,
                    primer_pair.rev_primer.sensitivity,
                    primer_pair.rev_primer.specificity,
                    primer_pair.sensitivity,
                    primer_pair.specificity,
                    primer_pair.score,
                    len(primer_pair.target_accessions),
                )
            )
        
        output_path = path.joinpath('primer_pairs.csv')

        #Output csv
        with open(output_path, 'w', newline='') as csv_file: 
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(
                (
                    'fw_root',
                    'fw_len',
                    'fw_seq',
                    'fw_seq_rep', 
                    'fw_tm',
                    'fw_sens',
                    'fw_spec',
                    'rev_root',
                    'rev_len',
                    'rev_seq',
                    'rev_seq_rep',
                    'rev_tm',
                    'rev_sens',
                    'rev_spec',
                    'sens',
                    'spec',
                    'score',
                    'pp_seq_rep',
                )
            )
            csv_writer.writerows(primer_data)

class Blast: 
    """
    Class to handle the BLAST jobs for generated oligos.
    
    #TODO: blast() can only be called through multi_blast() if output is desired. 
    #      This behaviour needs to be fixed. 

    Attributes
    ----------
    blastdb : pathlib.Path
        Path to the BLASTdb that is being used for the blast runs. 
    blastdb_len : int
        Length of the BLASTdb. (# of sequences)

    Methods
    -------
    _get_blastdb_len : see name
    blast : handles the BLAST job
    multi_blast : runs multiprocessing BLAST jobs
    output : outputs BLAST results into independent CSV fiels
     
    """
    def __init__(self, blastdb: pathlib.Path) -> None:
        self.blastdb = blastdb
        self.blastdb_len = self._get_blastdb_len()

    def _get_blastdb_len(self) -> int: 
        """Get length of the BLASTdb (# of sequences)"""
        
        with open(self.blastdb) as blastdb_file: 
            data = blastdb_file.read()
            return data.count('>')

    def blast(self, oligos: list) -> pd.DataFrame: 
        """
        Execute BLAST job on all oligos.

        Parameters
        ----------
        oligos : list
            List of Oligo objects that are the queries for the BLAST jobs. 
        
        Return
        ------
        blast_output : pandas.DataFrame
            DataFrame holding BLAST results for all oligo sequences. 
        """

        #Generate the oligo temporary file.
        #Each oligo is its own FASTA entry in one FASTA temp file
        fasta = tempfile.NamedTemporaryFile(delete=True)
        for oligo in oligos:
            fasta.write(f">{str(oligo.id)}\n{str(oligo.seq)}\n".encode())
        fasta.seek(0)

        #Run the BLAST job
        args = [
            "blastn",
            "-task",
            "blastn-short",
            "-db",
            str(self.blastdb),
            "-num_alignments",
            str(self.blastdb_len),
            "-outfmt",
            "10 qacc sacc ssciname pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-query",
            fasta.name,
        ]

        #Handle BLAST output
        result = subprocess.run(args, capture_output=True)
        decoded = result.stdout.decode('utf-8')
        output = io.StringIO(decoded)
        
        #Output formatting into dataframe
        headers=[
            'qacc',
            'sacc',
            'ssciname',
            'pident',
            'qlen',
            'length',
            'mismatch', 
            'gapopen', 
            'qstart', 
            'qend', 
            'sstart', 
            'send', 
            'evalue', 
            'bitscore',
        ]

        fasta.close()

        blast_output = pd.read_csv(output, sep=',', header=None, names=headers)

        return blast_output

    def multi_blast(self, oligos: list, num_jobs: int) -> dict:
        """
        Function for calling BLAST in several independent processes. 

        Parameters
        ----------
        oligos : list
            List of Oligo objects that are the queries for the BLAST job. 
        num_jobs : int
            Number of processes to separate the BLAST job into. 

        Returns
        -------
        blast_results : dict
            Each key is the oligo.id, and corresponds to that specific oligo's BLAST results. 
        
        """ 
        def job_allocator(oligos: list, num_jobs: int) -> list:
            """Create job list based on number of processes to run."""
            
            list_size = floor(len(oligos)/num_jobs)
            remainder = len(oligos)%num_jobs
            job_list = []

            for group in range(num_jobs-1): 
                job_list.append(oligos[0+(list_size*group):list_size+(list_size*group)])

            #Remainder group is smaller than the other groups.
            #Accounts for case where # of oligos is not divisible by num_jobs
            job_list.append(oligos[(list_size*(num_jobs-1)):(list_size*num_jobs+remainder)])

            return job_list

        #Allocate the jobs
        job_list = job_allocator(oligos, num_jobs)

        #Run the BLAST
        pool = multiprocessing.Pool(num_jobs)
        results = pool.map(self.blast, job_list)

        #Combine and return
        blast_results = dict()
        for job in results: 
            list_oligo_ids = set(job['qacc'])
            for oligo_id in list_oligo_ids: 
                blast_results[str(oligo_id)]=job.loc[job['qacc']==oligo_id]

        return blast_results

    def output(self, blast_results: dict, path: pathlib.Path, tag: str) -> None: 
        """Output all of the BLAST results for each oligo into independent CSV files"""
        #Make path to store all of the blast results
        blast_folder_path = path.joinpath(f'{tag}_blast')
        blast_folder_path.mkdir(exist_ok=True)

        #Go through blast result dictionary and output all of the data
        for blast_result_key in blast_results: 
            blast_output_path = blast_folder_path.joinpath(f"{blast_result_key}_{tag}.csv")
            blast_output_file = open(blast_output_path, 'w')
            blast_results[blast_result_key].to_csv(blast_output_file)
            blast_output_file.close()
