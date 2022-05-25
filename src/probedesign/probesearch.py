from oligogenerator import probeGenerator, oligo, nemaBlast
from consensus import alignment
import pathlib
import argparse
import pandas as pd
import numpy as np
import time

def print_runtime(action) -> None:
    """ Print the time and some defined action. """
    print(f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}] {action}')

def parse_args(): 
    parser = argparse.ArgumentParser(description='probesearch.py - identify viable probes in an alignment for given target sequences')
    parser.add_argument(
        'target_alignment_path', 
        action='store', 
        type=pathlib.Path,
        help = 'Path to target alignment file, fasta format'
    )
    parser.add_argument(
        '--output', 
        metavar='output_directory',
        dest='output_path',
        default=None,
        action='store',
        type=pathlib.Path, 
        help='Output path',
    )
    parser.add_argument(
        '--target_start',
        metavar='target_start', 
        dest='target_start',
        default=1,
        action='store',
        type=int,
        help='Start coordinate of target region, 1-based coordinates'
    )
    parser.add_argument(
        '--target_end',
        metavar='target_end', 
        dest='target_end', 
        default=None,
        action='store',
        type=int,
        help='End coordinate of target region, 1-based coordinates'
    )
    parser.add_argument(
        '--min_primer_len',
        action='store',
        type=int, 
        default=17,
        dest='min_primer_len',
        help='Minimum primer length (default=17)'
    )
    parser.add_argument(
        '--max_primer_len',
        action='store',
        type=int, 
        default=22,
        dest='max_primer_len',
        help='Maximum primer length (default=22)'
    )
    parser.add_argument(
        '--no_sens_spec_check',
        action='store_true',
        dest='sens_spec_flag',
        help='Flag to not check the putative probes for their specificity and sensitivity'
    )
    #Arguments for specificity checking
    parser.add_argument(
        '--blastdb',
        action='store',
        type=str,
        dest='blastdb',
        default='',
        help='Name of blastdb'
    )
    parser.add_argument(
        '--blastdb_len',
        action='store',
        type=int,
        dest='blastdb_len',
        help='Length of blastdb'
    )
    args = parser.parse_args()
    if not args.output_path:
        args.output_path = args.target_alignment_path.parent
    #Note that coordinates are converted to 0-based half-open coordinates
    return (
        args.target_alignment_path, 
        args.output_path, 
        args.target_start-1, 
        args.target_end, 
        args.min_primer_len, 
        args.max_primer_len, 
        args.sens_spec_flag, 
        args.blastdb, 
        args.blastdb_len, 
    )

def main():
    #Arguments
    target_alignment_path, output_path, target_start, target_end, min_primer_len, max_primer_len, check_flag, blastdb, blastdb_len = parse_args()
    #Process the alignment
    target_alignment = alignment(target_alignment_path)
    target_consensus = target_alignment.get_consensus()
    target_accessions = target_alignment.get_accessions()
    #Generate Probes
    print_runtime("Start")
    pb_gen = probeGenerator(target_consensus, target_start, target_end, min_primer_len, max_primer_len)
    print("Generating probes...")
    pb_gen.get_probes()
    print("Probes finished!")
    #Do the specificity check
    if check_flag is False: 
        #Read target accessions
        #target_accessions = get_target_accessions(target_accession_path)
        #Generate BLAST results
        pb_blast = nemaBlast(blastdb, blastdb_len)
        blast_results = pb_blast.blast_all_proper(pb_gen.probes)
        #Output BLAST results
        pb_blast.output(blast_results, output_path)
        for probe in pb_gen.probes: 
            probe.calculate_sensitivity(blast_results[probe.id], target_accessions)
            probe.calculate_specificity(blast_results[probe.id], target_accessions, blastdb_len)
            probe.calculate_score()
        #Output probe list
        pb_gen.output(output_path)
    else: 
        pb_gen.output(output_path)
    print_runtime("End")
if __name__ == '__main__': 
    main()