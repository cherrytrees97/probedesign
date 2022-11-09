#Standard libraries
import pathlib
import argparse
import logging
#Third-party libraries
import pandas as pd
import numpy as np
#Local modules
from oligogenerator import ProbeGenerator, Blast
from consensus import Alignment

def parse_args(): 
    parser = argparse.ArgumentParser(
        description='probesearch.py - identify viable probes in an alignment for given target sequences'
    )
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
    #Arguments for specificity checking
    parser.add_argument(
        '--no_sens_spec_check',
        action='store_true',
        dest='sens_spec_flag',
        help='Flag to not check the putative probes for their specificity and sensitivity'
    )
    parser.add_argument(
        '--blastdb',
        action='store',
        type=pathlib.Path,
        dest='blastdb',
        default='',
        help='Name of blastdb'
    )
    parser.add_argument(
        '--mp_job',
        '-m',
        action='store',
        type=int,
        default=1,
        dest='num_jobs',
        help='Number of processes to spawn to handle BLAST jobs. (Default=1)'
    )

    args = parser.parse_args()
    if not args.output_path:
        args.output_path = args.target_alignment_path.parent
    
    #Note that coordinates are converted to 0-based half-open coordinates
    args.target_start = args.target_start - 1
    
    return args

def main():
    #Arguments
    args = parse_args()
    
    #Set up logger
    logging.basicConfig(
        encoding='utf-8',
        level=logging.INFO,
        handlers=[
            logging.FileHandler(args.output_path.joinpath('probesearch.log')),
            logging.StreamHandler()
        ],
        format='%(levelname)s:%(message)s'
    )

    #Process the alignment
    logging.info(f'Probesearch - Designing probes for {args.target_alignment_path.stem}.')

    target_alignment = Alignment(args.target_alignment_path)
    target_alignment.get_consensus()

    #Generate Probes
    logging.info(f'Generating probes...')

    pb_gen = ProbeGenerator(target_alignment.consensus, args.target_start, args.target_end, args.min_primer_len, args.max_primer_len)
    pb_gen.get_probes()

    logging.info(f'Generated {len(pb_gen.probes)} probes.')

    #Do the specificity check
    if args.sens_spec_flag is False:     
        #Generate BLAST results
        logging.info(f'BLASTing probes...')

        pb_blast = Blast(args.blastdb)
        blast_results = pb_blast.multi_blast(pb_gen.probes, args.num_jobs)

        #Output BLAST results
        pb_blast.output(blast_results, args.output_path, 'probe')

        #Calculate sensitivity and specificity
        logging.info(f'Calculating sensitivity and specificity...')

        for probe in pb_gen.probes: 
            probe.calculate_sensitivity(target_alignment)
            probe.calculate_specificity(target_alignment, blast_results[probe.id], pb_blast.blastdb_len)
            probe.calculate_score()

    pb_gen.output(args.output_path)
    logging.info(f'Finished!')

if __name__ == '__main__': 
    main()
