from consensus import Alignment
from oligogenerator import PrimerGenerator, Blast
import pathlib
import argparse
import time
import logging

def parse_args(): 
    parser = argparse.ArgumentParser(description='Search for primers')

    parser.add_argument(
        'target_alignment_path', 
        action='store', 
        type=pathlib.Path,
        help = 'Path to target alignment file, fasta format'
    )
    parser.add_argument(
        'pb_start',
        metavar='pb_start', 
        action='store', 
        type=int, 
        help = 'Start coordinate of probe, 1-based coordinates'
    )
    parser.add_argument(
        'pb_len',
        metavar='pb_len', 
        action='store', 
        type=int, 
        help = 'Length of the probe'
    )

    primer_param = parser.add_argument_group('Primer parameters')
    primer_param.add_argument(
        '--min_primer_len',
        action='store',
        type=int, 
        default=17,
        dest='min_primer_len',
        help='Minimum primer length'
    )
    primer_param.add_argument(
        '--max_primer_len',
        action='store',
        type=int, 
        default=22,
        dest='max_primer_len',
        help='Maximum primer length'
    )
    
    blast_param = parser.add_argument_group('BLAST parameters')
    blast_param.add_argument(
        '--no_sens_spec',
        action='store_true',
        dest='sens_spec_flag',
        help='Flag to not check the putative probes for their specificity and sensitivity'
    )
    blast_param.add_argument('--blastdb',
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

    output_param = parser.add_argument_group('Output parmaeters')
    output_param.add_argument(
        '--output_path',
        '-o',
        action='store',
        type=pathlib.Path, 
        default=None,
        dest='output_path', 
        help='Output path',
    )

    filter_param = parser.add_argument_group('Output Filter parameters')
    filter_param.add_argument(
        '--filter_tm',
        '-ft',
        action='store_true',
        dest='f_tm',
        help='Filter primers by tm diff, as well as min and max tm'
    )
    filter_param.add_argument(
        '--max_tm_diff'
        '-d',
        dest='tm_diff',
        metavar='maximum_tm_diff',
        action='store',
        type=float,
        default=5.0,
        help='Maximum temperature difference between forward and reverse'
    )
    filter_param.add_argument(
        '--max_primer_tm',
        action='store',
        dest='max_tm',
        type=float,
        default=63.0,
        help='Maximum tm for primer.'
    )
    filter_param.add_argument(
        '--min_primer_tm',
        action='store',
        dest='min_tm',
        type=float,
        default=55.0,
        help='Minimum tm for primer.'
    )
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

    args = parser.parse_args()

    if not args.output_path: 
        args.output_path = args.target_alignment_path.parent

    #Arguments
    #Note conversion of pb_start to 0-based coordinate system
    args.pb_start = args.pb_start - 1

    return args

def get_param_string(args: argparse.Namespace) -> str: 
    """Returns formatted string that describes program parameters."""
    param_string = (
        "\nParameters for this run: \n"
        f"\tInput alignment: {args.target_alignment_path}\n"
        f"\tInput probe start: {args.pb_start}\n"
        f"\tInput probe length: {args.pb_len}\n"
        f"Primer parameters\n"
        f"\tMin primer len: {args.min_primer_len}\n"
        f"\tMax primer len: {args.max_primer_len}\n"
        f"BLAST parameters\n"
        f"\tCheck sens/spec: {args.sens_spec_flag}\n"
        f"\tBLASTdb: {args.blastdb}\n"
        f"\tAllocated cores: {args.num_jobs}\n"
        f"Filter parameters\n"
        f"\tFilter by sequence representation: {args.f_seq_rep}\n"
        f"\tMin sequence rep: {args.min_seq_rep}\n"
        f"\tFilter by Tm: {args.f_tm}\n"
        f"\tMin primer Tm: {args.min_tm}\n"
        f"\tMax primer Tm: {args.max_tm}\n"
        f"\tMax Tm diff: {args.tm_diff}\n"
        f"Output parameters\n"
        f"\tOutput path: {args.output_path}"        
    )

    return param_string

def main(): 
    #Get arguments
    args = parse_args()

    #Set up logger
    logging.basicConfig(
        encoding='utf-8',
        level=logging.INFO,
        handlers=[
            logging.FileHandler(args.output_path.joinpath('primersearch.log')),
            logging.StreamHandler()
        ],
        format='%(asctime)s:%(levelname)s: %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S',
    )

    #Process the alignment
    logging.info(f'Primersearch.py - designing primer pairs for probe {args.pb_start}-{args.pb_len}')
    logging.info(f"{get_param_string(args)}")

    target_alignment = Alignment(args.target_alignment_path)
    target_alignment.get_consensus()
    num_seq = len(target_alignment.alignment)
    
    logging.info("Generating primers...")
    
    primer_gen = PrimerGenerator(
        target_alignment.consensus, 
        args.pb_start, 
        args.pb_len, 
        args.min_primer_len, 
        args.max_primer_len, 
        args.tm_diff
    )

    #Generate primer pairs
    primer_gen.find_fw_primers()
    primer_gen.find_rev_primers()

    logging.info("Primers finished!")
    logging.info(f"Total number of forward primers: {len(primer_gen.fw_primers)}")
    logging.info(f"Total number of reverse primers: {len(primer_gen.rev_primers)}")

    #Generate BLAST results
    if args.sens_spec_flag is False: 
        primer_blast = Blast(args.blastdb)

        logging.info(f"Generating BLAST results...")
        fw_blast_results = primer_blast.multi_blast(primer_gen.fw_primers, args.num_jobs)
        rev_blast_results = primer_blast.multi_blast(primer_gen.rev_primers, args.num_jobs)
        logging.info(f"Done!")

        logging.info(f"Calculating sensitivity and specificity scores...")
        for fw_primer in primer_gen.fw_primers: 
            fw_primer.calculate_sensitivity(target_alignment)
            fw_primer.calculate_specificity(
                target_alignment, 
                fw_blast_results[fw_primer.id], 
                primer_blast.blastdb_len,
                )
            fw_primer.calculate_score()

        for rev_primer in primer_gen.rev_primers: 
            rev_primer.calculate_sensitivity(target_alignment, reverse=True)
            rev_primer.calculate_specificity(
                target_alignment,
                rev_blast_results[rev_primer.id], 
                primer_blast.blastdb_len
                )
            rev_primer.calculate_score()
        logging.info(f"Done!")

        logging.info(f"Outputting BLAST results...")
        primer_blast.output(fw_blast_results, args.output_path, 'fw')
        primer_blast.output(rev_blast_results, args.output_path, 'rev')
        logging.info(f"Done!")

    #Generate primer pairs
    logging.info(f"Finding primer pairs...")
    primer_gen.find_primer_pairs()
    logging.info(f"Generated {len(primer_gen.primer_pairs)} primer pairs!")

    if args.sens_spec_flag is False:
        #Generate primer pair data
        logging.info(f"Calculating sensitivity and specificity scores for primer pairs...")
        for primer_pair in primer_gen.primer_pairs: 
            fw_blast = fw_blast_results[primer_pair.fw_primer.id]
            rev_blast = rev_blast_results[primer_pair.rev_primer.id]
            primer_pair.calculate_sensitivity()
            primer_pair.calculate_specificity(fw_blast, rev_blast, primer_blast.blastdb_len)
            primer_pair.calculate_score()
        logging.info(f"Done!")

    #Output
    primer_gen.output(
        args.output_path,
        args.f_tm,
        args.tm_diff,
        args.max_tm,
        args.min_tm,
        args.f_seq_rep,
        args.min_seq_rep,
        num_seq,
    )
    logging.info(f"Program done!")

if __name__ == '__main__':
    main()