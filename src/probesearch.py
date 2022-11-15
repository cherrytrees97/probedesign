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
    probe_param = parser.add_argument_group('Probe parameters')
    probe_param.add_argument(
        '--target_start',
        metavar='target_start', 
        dest='target_start',
        default=1,
        action='store',
        type=int,
        help='Start coordinate of target region, 1-based coordinates'
    )
    probe_param.add_argument(
        '--target_end',
        metavar='target_end',
        dest='target_end',
        default=None,
        action='store',
        type=int,
        help='End coordinate of target region, 1-based coordinates'
    )
    probe_param.add_argument(
        '--min_probe_len',
        action='store',
        type=int,
        default=17,
        dest='min_probe_len',
        help='Minimum primer length (default=17)'
    )
    probe_param.add_argument(
        '--max_probe_len',
        action='store',
        type=int,
        default=22,
        dest='max_probe_len',
        help='Maximum primer length (default=22)'
    )
    blast_param = parser.add_argument_group('BLAST parameters')
    blast_param.add_argument(
        '--no_sens_spec_check',
        action='store_true',
        dest='sens_spec_flag',
        help='Flag to not check the putative probes for their specificity and sensitivity'
    )
    blast_param.add_argument(
        '--blastdb',
        action='store',
        type=pathlib.Path,
        dest='blastdb',
        default='',
        help='Name of blastdb'
    )
    blast_param.add_argument(
        '--mp_job',
        '-m',
        action='store',
        type=int,
        default=1,
        dest='num_jobs',
        help='Number of processes to spawn to handle BLAST jobs. (Default=1)'
    )
    output_param = parser.add_argument_group('Output parameters')
    output_param.add_argument(
        '--output',
        '-o',
        metavar='output_directory',
        dest='output_path',
        default=None,
        action='store',
        type=pathlib.Path, 
        help='Output path',
    )
    filter_param = parser.add_argument_group('Filter parameters')
    filter_param.add_argument(
        '--filter_seq_rep',
        '-fs',
        action='store_true',
        dest='f_seq_rep',
        help='Filter by probes returned by sequence representation'
    )
    filter_param.add_argument(
        '--filter_min',
        action='store',
        dest='min_seq_rep',
        default=0.5,
        type=float,
        help='Minimum percentage of sequences that need to be represented for probes to be returned. Default = 0.5'
    )
    filter_param.add_argument(
        '--filter_tm',
        '-ft',
        action='store_true',
        dest='f_tm',
        help='Filter probes by tm range'
    )
    filter_param.add_argument(
        '--min_tm',
        action='store',
        dest='min_tm',
        default=68.0,
        type=float,
        help='Minimum tm of probe for it to be returned. Default = 68.0'
    )
    filter_param.add_argument(
        '--max_tm',
        action='store',
        dest='max_tm',
        default=70.0,
        type=float,
        help='Maximum tm of probe for it to be returned. Default = 70.0'
    )

    args = parser.parse_args()

    if not args.output_path:
        args.output_path = args.target_alignment_path.parent
    
    #Note that coordinates are converted to 0-based half-open coordinates
    args.target_start = args.target_start - 1
    
    return args

def get_param_string(args: argparse.Namespace) -> str: 
    """Returns formatted string that describes program parameters."""
    param_string = (
        "\nParameters for this run: \n"
        f"\tInput alignment: {args.target_alignment_path}\n"
        f"Probe parameters\n"
        f"\tProbe search start: {args.target_start}\n"
        f"\tProbe search end: {args.target_end}\n"
        f"\tMin probe len: {args.min_probe_len}\n"
        f"\tMax probe len: {args.max_probe_len}\n"
        f"BLAST parameters\n"
        f"\tCheck sens/spec: {args.sens_spec_flag}\n"
        f"\tBLASTdb: {args.blastdb}\n"
        f"\tAllocated cores: {args.num_jobs}\n"
        f"Filter parameters\n"
        f"\tFilter by sequence representation: {args.f_seq_rep}\n"
        f"\tMin sequence rep: {args.min_seq_rep}\n"
        f"\tFilter by Tm: {args.f_tm}\n"
        f"\tMin probe Tm: {args.min_tm}\n"
        f"\tMax probe Tm: {args.max_tm}\n"
        f"Output parameters\n"
        f"\tOutput path: {args.output_path}"        
    )

    return param_string

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
        format='%(asctime)s:%(levelname)s: %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S',
    )

    #Process the alignment
    logging.info(f'Probesearch - Designing probes for {args.target_alignment_path.stem}.')
    logging.info(f"{get_param_string(args)}")

    target_alignment = Alignment(args.target_alignment_path)
    target_alignment.get_consensus()
    num_seq = len(target_alignment.alignment)

    #Generate Probes
    logging.info(f'Generating probes...')

    pb_gen = ProbeGenerator(
        target_alignment.consensus, 
        args.target_start, 
        args.target_end, 
        args.min_probe_len, 
        args.max_probe_len)
    pb_gen.get_probes()

    logging.info(f'Generated {len(pb_gen.probes)} probes.')

    #Do the specificity check
    if args.sens_spec_flag is False:     
        #Generate BLAST results
        logging.info(f'BLASTing probes...')

        pb_blast = Blast(args.blastdb)
        blast_results = pb_blast.multi_blast(pb_gen.probes, args.num_jobs)

        logging.info(f'BLAST jobs complete!')
        #Output BLAST results
        logging.info(f'Outputting BLAST result .csv to {str(args.output_path)}')

        pb_blast.output(blast_results, args.output_path, 'probe')

        logging.info(f'Output complete.')

        #Calculate sensitivity and specificity
        logging.info(f'Calculating sensitivity and specificity...')

        for probe in pb_gen.probes: 
            probe.calculate_sensitivity(target_alignment)
            probe.calculate_specificity(target_alignment, blast_results[probe.id], pb_blast.blastdb_len)
            probe.calculate_score()

    pb_gen.output(
        args.output_path, 
        args.f_seq_rep, 
        args.min_seq_rep, 
        num_seq,
        args.f_tm,
        args.min_tm,
        args.max_tm,
        )

    logging.info(f'Finished!')

if __name__ == '__main__': 
    main()
