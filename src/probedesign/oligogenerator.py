from Bio import SeqIO
import csv
import tempfile
import multiprocessing
import subprocess
import io
import pandas as pd
import primer3
from math import floor
from functools import partial

class oligo: 
    def __init__(self, root_pos, seq):
        self.seq = seq
        self.root_pos = root_pos
        self.len = len(self.seq)
        #self.tm = tm
        self.id =  f"{str(self.root_pos)}-{str(self.len)}"
        self.sensitivity = 0
        self.specificity = 0
        self.score = 0
    def calculate_sensitivity(self, blast_results, target_accessions):
        """
        Function calculates sensitivity of the oligo (binding to target sequences). 
        Binding is defined as 100% query coverage and 100% percent identity of the oligo to the target sequence.
        Sensitivity is calculated by the following formula: 
            TP / (TP + FN)
            where: 
            TP = succesfully amplified accessions
            FN = target accessions that were not amplified
            TP + FN = total number of target accessions
        """
        #Take only accessions where there was a perfect match --> full query coverage, 100% identity
        perfect_match_results = blast_results.loc[(blast_results['qlen']==blast_results['length'])&(blast_results['pident']==100.0)]
        #Retrieve only the accessions list
        amplified_accessions = set(perfect_match_results.loc[:,'sacc'])
        target_match = 0
        #Count number of target accessions in the amplified accession list 
        for accession in amplified_accessions: 
                if accession in target_accessions: 
                    target_match = target_match + 1
        #Calculate sensitivity and return
        self.sensitivity = target_match/len(target_accessions)
    def calculate_specificity(self, blast_results, target_accessions, blastdb_len): 
        """
        Function calculates specificity of the oligo (binding to non-target sequences). 
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
        """
        blast_match_accessions = set(blast_results.loc[:,'sacc'])
        #Remove every target_accession from the blast_match_accessions list
        for accession in target_accessions:
            if accession in blast_match_accessions:  
                blast_match_accessions.remove(accession)
        #Calculate specificity
        #Total non-target = all_blast_sequences - target_accessions
        self.specificity = (blastdb_len - len(blast_match_accessions))/(blastdb_len - len(target_accessions))
    def calculate_score(self): 
        self.score = self.sensitivity + self.specificity

class primerpair: 
    def __init__(self, fw_primer, rev_primer): 
        self.fw_primer = fw_primer
        self.rev_primer = rev_primer
        self.sensitivity = 0
        self.specificity = 0
        self.score = 0 
    def calc_tm_diff(self): 
        return abs(self.fw_primer.tm - self.rev_primer.tm)
    def calculate_sensitivity(self, fw_blast_results, rev_blast_results, target_accessions): 
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
        """
        #Take only accessions where there was a perfect match --> full query coverage, 100% identity
        fw_perfect_match_results = fw_blast_results.loc[(fw_blast_results['qlen']==fw_blast_results['length'])&(fw_blast_results['pident']==100.0)]
        rev_perfect_match_results = rev_blast_results.loc[(rev_blast_results['qlen']==rev_blast_results['length'])&(rev_blast_results['pident']==100.0)]
        #Retrieve only the accessions list
        fw_perfect_match_accessions = set(fw_perfect_match_results.loc[:,'sacc'])
        rev_perfect_match_accessions = set(rev_perfect_match_results.loc[:,'sacc'])
        amplified_accessions = set(fw_perfect_match_accessions&rev_perfect_match_accessions)
        target_match = 0
        #Count number of target accessions in the amplified accession list 
        for accession in amplified_accessions: 
                if accession in target_accessions: 
                    target_match = target_match + 1
        #Calculate sensitivity and return
        self.sensitivity = target_match/len(target_accessions)
    def calculate_specificity(self, fw_blast_results, rev_blast_results, target_accessions, blastdb_len):
        fw_match_accessions = set(fw_blast_results.loc[:,'sacc'])
        rev_match_accessions = set(rev_blast_results.loc[:,'sacc'])
        amplified_accessions = set(fw_match_accessions&rev_match_accessions)
        #Remove every target_accession from the blast_match_accessions list
        for accession in target_accessions:
            if accession in amplified_accessions:  
                amplified_accessions.remove(accession)
        #Calculate specificity
        #Total non-target = all_blast_sequences - target_accessions
        self.specificity = (blastdb_len - len(amplified_accessions))/(blastdb_len - len(target_accessions))
    def calculate_score(self): 
        self.score = self.sensitivity + self.specificity

class probeGenerator: 
    def __init__(
        self, 
        template, 
        start, 
        end, 
        min_length, 
        max_length, 
        ):
        #Note that the coordinates are 1-based, so the template start has to be -1
        self.template = template[start:end]
        self.start = start
        self.end = end
        self.min_length = min_length
        self.max_length = max_length
        self.probes = []
    def get_probes(self): 
        """
        Function finds all viable probes. 
        Probe criteria: 
        1) 5' end of the probe is not a G
        2) More C's than G's
        3) No more than 4 nucleotides in a row
        4) 30 - 80% CG content
        """
        def check_probe(seq): 
            def count_c_g(seq): 
                return seq.count('c') + seq.count('g') + seq.count('C') + seq.count('G')
            def chk_run(seq): 
                """
                Ensure that there are no runs of four or more identical nucleotides in the probe
                """
                current = seq[0]
                identical_len = 1
                detected = False
                #Go through each nucleotide in the primer
                for nucl in range(1, len(seq)): 
                    if seq[nucl] == current: 
                        identical_len = identical_len + 1
                        if identical_len > 3: 
                            detected = True
                            break
                    else: 
                        current = seq[nucl]
                        identical_len = 1
                return detected
            def chk_last_5(seq): 
                return count_c_g(seq[-5:]) > 2
            #Check the three conditions
            percent_c_g = count_c_g(seq)/len(seq)
            run_flag = chk_run(seq)
            last5_flag = chk_last_5(seq)
            if(
                0.3 <= percent_c_g <= 0.8 
                and run_flag is False
                and last5_flag is False
                and seq[0] != 'G'
                and seq[0] != 'g'
                and not('n' in seq)
                and not('N' in seq)
            ):
                return True
            else: 
                return False
        
        #For each root position in the target sequence
        for i in range(len(self.template) - self.min_length + 1):
            #Iterate over all probe lengths
            for probe_len in range(self.min_length, self.max_length + 1): 
                probe_seq = self.template[i:i+probe_len]
                if check_probe(probe_seq) is True: 
                    #Note that the coordinates are converted back to 1-based
                    self.probes.append(oligo(self.start+i+1, probe_seq))
    def output(self, path): 
        probe_data = []
        for probe in self.probes: 
            probe_data.append(
                (
                    probe.root_pos,
                    probe.len, 
                    probe.seq, 
                    probe.sensitivity, 
                    probe.specificity,
                    probe.score,
                )
            )
        output_path = path.joinpath('probe_candidates.csv')
        csv_file = open(output_path, 'w', newline='')
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(
            (
                'probe_root',
                'probe_len',
                'probe_seq',
                'sens',
                'spec',
                'score',
            )
        )
        csv_writer.writerows(probe_data)
        csv_file.close()

class primerGenerator: 
    def __init__(
        self,
        template_seq_path, 
        pb_start, 
        pb_len,
        min_length, 
        max_length, 
        max_tm_diff
        ):
        #Input properties
        template_seq_file = SeqIO.read(template_seq_path, 'fasta')
        self.template = template_seq_file.seq
        self.pb_start = pb_start
        self.pb_len = pb_len
        self.pb_end = pb_start-1+pb_len
        self.min_length = min_length
        self.max_length = max_length
        self.max_tm_diff = max_tm_diff
        #Output properties
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
        """
        def check_primer(seq): 
            def count_c_g(seq): 
                return seq.count('c') + seq.count('g') + seq.count('C') + seq.count('G')
            def chk_run(seq): 
                """
                Ensure that there are no runs of four or more identical nucleotides in the probe
                """
                current = seq[0]
                identical_len = 1
                detected = False
                #Go through each nucleotide in the primer
                for nucl in range(1, len(seq)): 
                    if seq[nucl] == current: 
                        identical_len = identical_len + 1
                        if identical_len > 3: 
                            detected = True
                            break
                    else: 
                        current = seq[nucl]
                        identical_len = 1
                return detected
            def chk_last_5(seq): 
                return count_c_g(seq[-5:]) > 2
            #Check the three conditions
            percent_c_g = count_c_g(seq)/len(seq)
            run_flag = chk_run(seq)
            last5_flag = chk_last_5(seq)
            #print(f"{str(percent_c_g)}{str(run_flag)} {str(last5_flag)}")
            #Return true if all of the conditions are passed, otherwise return false
            if(
                0.3 <= percent_c_g <= 0.8 
                and run_flag is False
                and last5_flag is False
            ):
                #print('Passed!')
                return True
            else: 
                return False
        
        #Search for each root position and each legal primer length
        for i in range(50): 
            for fw_len in range(self.min_length, self.max_length+1):
                fw_start = self.pb_start - i - fw_len
                fw_end = self.pb_start - i
                fw_primer_seq = str(self.template[fw_start:fw_end])
                if check_primer(fw_primer_seq) is True: 
                    self.fw_primers.append(
                        oligo(
                            self.pb_start-i, 
                            fw_primer_seq, 
                            float(primer3.calcTm(fw_primer_seq, dv_conc=1.5))
                        )
                    )
    def find_rev_primers(self):
        """
        Function finds all viable REV primers.
        REV primer criteria: 
        1) Between min and max length 
        2) % GC is 30% to 80%
        3) Last five nucleotides at the 3' end contain no more than two G + C residues
        4) No more than 4 consecutive nucleotides within the primer 
        
        Input data: 
        Input data: 
        1) target_seq - str - target sequence - ATCTGATCATGATCATGACTAGTCATGGC
        2) pb_start - int - start index of 5'-end of the probe - 607
        3) pb_end - int - index of 3'-end of the probe (1 past index of last probe nucleotide) - 623 

        Output data: 
        [
            {root_pos:0, len:17, seq:''}
        ]

        Algorithm: 
        The 1st root position is the pb_end index. 
        The last legal root position is the case where: 
            1) fw_primer is at the MIN_PRIMER_LEN
            2) rev_primer is at the MIN_PRIMER_LEN
            3) fw_primer is directly adjacent to probe (no gap)
        To calculate the last legal root position: 
            150 - 2(MIN_PRIMER_LEN) - pb_len
        Let the last legal position be x, min_primer_len be 17, and max_primer len be 22. 
        From 1 to X - (max_primer_len - min_primer_len + 1), there are max_primer_len - min_primer_len + 1 primers. 
        For each subsequent position, the number of primers decreases by 1 until the last position, 
        where the only legal primer length is min_primer_len. 

        {root_pos:pb_end + i, len:rev_len, seq:[(pb_end+i):(pb_end+i+pb_len)]}
        """
        def check_primer(seq): 
            def count_c_g(seq): 
                return seq.count('c') + seq.count('g') + seq.count('C') + seq.count('G')
            def chk_run(seq): 
                """
                Ensure that there are no runs of four or more identical nucleotides in the probe
                """
                current = seq[0]
                identical_len = 1
                detected = False
                #Go through each nucleotide in the primer
                for nucl in range(1, len(seq)): 
                    if seq[nucl] == current: 
                        identical_len = identical_len + 1
                        if identical_len > 3: 
                            detected = True
                            break
                    else: 
                        current = seq[nucl]
                        identical_len = 1
                return detected
            def chk_last_5(seq): 
                return count_c_g(seq[-5:]) > 2
            #Check the three conditions
            percent_c_g = count_c_g(seq)/len(seq)
            run_flag = chk_run(seq)
            last5_flag = chk_last_5(seq)
            #print(f"{str(percent_c_g)}{str(run_flag)} {str(last5_flag)}")
            #Return true if all of the conditions are passed, otherwise return false
            if(
                0.3 <= percent_c_g <= 0.8 
                and run_flag is False
                and last5_flag is False
            ):
                #print('Passed!')
                return True
            else: 
                return False
          
        last_root_pos = 150 - 2*self.min_length - self.pb_len
        #Case 1: all primer lengths are available at these root positions
        for i in range(last_root_pos - (self.max_length) - self.min_length + 1): 
            for rev_len in range(self.min_length, self.max_length+1): 
                rev_start = self.pb_end + i
                rev_end = rev_start + rev_len
                rev_primer_seq = str(self.template[rev_start:rev_end].reverse_complement())
                if check_primer(rev_primer_seq) is True: 
                    self.rev_primers.append(
                        oligo(
                            self.pb_end + i, 
                            rev_primer_seq, 
                            float(primer3.calcTm(rev_primer_seq, dv_conc=1.5))
                        )
                    )
        #Case 2: end of root positions
        for i in range(self.max_length - self.min_length + 1):
            rev_start = self.pb_end + last_root_pos - i
            for rev_len in range(self.min_length, self.min_length + 1 + i):
                rev_end = rev_start + rev_end
                rev_primer_seq = str(self.template[rev_start:rev_end].reverse_complement())
                if check_primer(rev_primer_seq) is True: 
                    self.rev_primers.append(
                        oligo(
                            self.pb_end + i, 
                            rev_primer_seq, 
                            float(primer3.calcTm(rev_primer_seq, dv_conc=1.5))
                        )
                    )
    def find_primer_pairs(self): 
        """
        Function will pair the fw and rev primers together. 
        The function iterates over the list of FW primers, and pairs them with all viable reverse primers. 
        Viable reverse primers --> a reverse primer such that the amplicon length is equal to or less than 
        150 bp. 

        Algorithm: 
        Amplicon length = (rev_root_pos + rev_len) - fw_root_pos
        Root_pos_limit = FW_root_pos + 150 - MIN_PRIMER_LEN

        """
        def check_tm_diff(fw_tm, rev_tm):
            return abs(fw_tm-rev_tm)

        for fw_primer in self.fw_primers: 
            root_pos_limit = fw_primer.root_pos + 150 - self.min_length
            for rev_primer in self.rev_primers:
                if (
                    rev_primer.root_pos < root_pos_limit
                    and (rev_primer.root_pos + rev_primer.len - fw_primer.root_pos) <= 150
                ):  
                    if check_tm_diff(fw_primer.tm, rev_primer.tm) <= self.max_tm_diff:
                        self.primer_pairs.append(
                            primerpair(
                                fw_primer, 
                                rev_primer
                            )
                    )
                else: 
                    break
    def output(self, path): 
        primer_data = []
        for primer_pair in self.primer_pairs: 
            primer_data.append(
                (
                    primer_pair.fw_primer.root_pos,
                    primer_pair.fw_primer.len,
                    primer_pair.fw_primer.seq, 
                    primer_pair.fw_primer.tm,
                    primer_pair.fw_primer.sensitivity,
                    primer_pair.fw_primer.specificity,
                    primer_pair.rev_primer.root_pos,
                    primer_pair.rev_primer.len,
                    primer_pair.rev_primer.seq, 
                    primer_pair.rev_primer.tm,
                    primer_pair.rev_primer.sensitivity,
                    primer_pair.rev_primer.specificity,
                    primer_pair.sensitivity,
                    primer_pair.specificity,
                    primer_pair.score,
                )
            )
        
        output_path = path.with_name('primer_pairs.csv')
        csv_file = open(output_path, 'w', newline='')
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(
            (
                'fw_root',
                'fw_len',
                'fw_seq',
                'fw_tm',
                'fw_sens',
                'fw_spec',
                'rev_root',
                'rev_len',
                'rev_seq',
                'rev_tm',
                'rev_sens',
                'rev_spec',
                'sens',
                'spec',
                'score'
            )
        )
        csv_writer.writerows(primer_data)
        csv_file.close()

class blast: 
    def __init__(self, blastdb, blastdb_len):
        self.blastdb = blastdb
        self.blastdb_len = blastdb_len
        self.data = None
        self.NUM_POOL = 8
    def blast_all(self, oligos): 
        #Generate the oligo temporary file
        fasta = tempfile.NamedTemporaryFile(delete=True)
        print('Oligos in the temp file: ')
        for oligo in oligos:
            fasta.write(f">{str(oligo.id)}\n{str(oligo.seq)}\n".encode())
            print(oligo.id)
            print(oligo.seq)
        fasta.seek(0)
        #cpu_count = multiprocessing.cpu_count() - 2
        cpu_count = 4
        print(f"NUMBER OF THREADS: {str(cpu_count)}")
        #Run the BLAST job
        args = [
            "blastn",
            "-task",
            "blastn-short",
            "-db",
            self.blastdb,
            "-num_alignments",
            str(self.blastdb_len),
            "-outfmt",
            "10 qacc sacc ssciname pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-query",
            fasta.name,
            "-num_threads",
            str(cpu_count),
            "-mt_mode",
            str(1)
        ]
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
        self.data = pd.read_csv(output, sep=',', header=None, names=headers)
        fasta.close()
        #Split the data
        list_oligo_ids = set(self.data['qacc'])
        #print(f"qacc set:\n{list_oligo_ids}")
        blast_results = dict()
        for oligo_id in list_oligo_ids: 
            blast_results[str(oligo_id)]=self.data.loc[self.data['qacc']==oligo_id]
        return blast_results
    def blast(self, oligos): 
        #Generate the oligo temporary file
        fasta = tempfile.NamedTemporaryFile(delete=True)
        for oligo in oligos:
            fasta.write(f">{str(oligo.id)}\n{str(oligo.seq)}\n".encode())
        fasta.seek(0)
        #cpu_count = multiprocessing.cpu_count() - 2
        cpu_count = 16/self.NUM_POOL
        #Run the BLAST job
        args = [
            "blastn",
            "-task",
            "blastn-short",
            "-db",
            self.blastdb,
            "-num_alignments",
            str(self.blastdb_len),
            "-outfmt",
            "10 qacc sacc ssciname pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-query",
            fasta.name,
            "-num_threads",
            str(cpu_count),
            "-mt_mode",
            str(1)
        ]
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
        return(pd.read_csv(output, sep=',', header=None, names=headers))
    def multi_blast(self, oligos): 
        def job_allocator(oligos, NUM_GROUPS):
            list_size = floor(len(oligos)/NUM_GROUPS)
            remainder = len(oligos)%NUM_GROUPS
            job_list = []
            for group in range(NUM_GROUPS-1): 
                job_list.append(oligos[0+(list_size*group):list_size+(list_size*group)])
            job_list.append(oligos[(list_size*(NUM_GROUPS-1)):(list_size*NUM_GROUPS+remainder)])
            return job_list
        NUM_POOL = 8
        #Allocate the jobs
        job_list = job_allocator(oligos, NUM_POOL)
        print(len(job_list))
        for item in job_list: 
            print(len(item))
            print(item[0].id)
        #Run the BLAST
        pool = multiprocessing.Pool(NUM_POOL)
        results = pool.map(self.blast, job_list)
        #Combine and return
        blast_results = dict()
        for job in results: 
            list_oligo_ids = set(job['qacc'])
            for oligo_id in list_oligo_ids: 
                blast_results[str(oligo_id)]=job.loc[job['qacc']==oligo_id]
        return blast_results
    def output(self, blast_results, path): 
        #Make path to store all of the blast results
        blast_folder_path = path.joinpath('probe_blast')
        blast_folder_path.mkdir(exist_ok=True)
        #Go through blast result dictionary and output all of the data
        for blast_result_key in blast_results: 
            blast_output_path = blast_folder_path.joinpath(f"{blast_result_key}_probe.csv")
            blast_output_file = open(blast_output_path, 'w')
            blast_results[blast_result_key].to_csv(blast_output_file)
            blast_output_file.close()
