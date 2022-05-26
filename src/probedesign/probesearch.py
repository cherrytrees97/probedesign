from oligogenerator import ProbeGenerator, Blast
from consensus import Alignment
import pathlib
import argparse
import pandas as pd
import numpy as np
import time

def print_runtime(action) -> None:
    """ Print the time and some defined action. """
    print(f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}] {action}')

def output_timers(timers): 
    print(f"Total runtime: {timers['end'] - timers['start']}")
    print(f"Probe generation runtime: {timers['pb_gen-end'] - timers['pb_gen-start']}")
    print(f"BLAST runtime: {timers['blast-end'] - timers['blast-start']}")
    print(f"Calculation runtime: {timers['calc-end'] - timers['calc-start']}")

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
        type=pathlib.Path,
        dest='blastdb',
        default='',
        help='Name of blastdb'
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
    )

def main():
    #Timer
    timers ={
        "start":time.monotonic(),
    }
    #Arguments
    target_alignment_path, output_path, target_start, target_end, min_primer_len, max_primer_len, check_flag, blastdb, = parse_args()
    #Process the alignment
    target_alignment = Alignment(target_alignment_path)
    target_alignment.get_consensus()
    target_accessions = target_alignment.get_accessions()

    #Generate Probes
    print_runtime("Start")
    pb_gen = ProbeGenerator(target_alignment.consensus, target_start, target_end, min_primer_len, max_primer_len)
    print("Generating probes...")
    timers['pb_gen-start'] = time.monotonic()
    pb_gen.get_probes()
    timers['pb_gen-end'] = time.monotonic()
    print("Probes finished!")

    print(f"Total number of probes generated: {len(pb_gen.probes)}")

    #Do the specificity check
    if check_flag is False: 
        #Read target accession        
        #Generate BLAST results
        timers['blast-start'] = time.monotonic()
        pb_blast = Blast(blastdb)
        #blast_results = pb_blast.blast_all(pb_gen.probes)
        blast_results = pb_blast.multi_blast(pb_gen.probes)
        timers['blast-end'] = time.monotonic()
        print("Blast complete.")
        #Output BLAST results
        print("Outputting BLAST results...")
        pb_blast.output(blast_results, output_path, 'probe')
        print("Output complete...")
        print("Calculating sensitivity and specificity...")
        timers['calc-start'] = time.monotonic()
        for probe in pb_gen.probes: 
            probe.calculate_sensitivity(blast_results[probe.id], target_alignment.sequence_regions)
            probe.calculate_specificity(blast_results[probe.id], target_alignment.sequence_regions, pb_blast.blastdb_len)
            probe.calculate_score()
        timers['calc-end'] = time.monotonic()
        print("Calculation complete.")
        #Output probe list
        pb_gen.output(output_path)
    else: 
        pb_gen.output(output_path)
    print_runtime("End")
    timers['end'] = time.monotonic()
    output_timers(timers)
if __name__ == '__main__': 
    main()
