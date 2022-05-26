from consensus import alignment
from oligogenerator import primerGenerator, blast
import pathlib
import argparse

def print_runtime(action) -> None:
    """ Print the time and some defined action. """
    print(f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}] {action}')

def output_timers(timers): 
    print(f"Total runtime: {timers['end'] - timers['start']}")
    print(f"Probe generation runtime: {timers['pb_gen-end'] - timers['pb_gen-start']}")
    print(f"BLAST runtime: {timers['blast-end'] - timers['blast-start']}")
    print(f"Calculation runtime: {timers['calc-end'] - timers['calc-start']}")

def parse_args(): 
    parser = argparse.ArgumentParser(description='Search for primers')
    parser.add_argument('target_alignment_path', 
        action='store', 
        type=pathlib.Path,
        help = 'Path to target alignment file, fasta format'
    )
    parser.add_argument('pb_start',
        metavar='pb_start', 
        action='store', 
        type=int, 
        help = 'Start coordinate of probe, 1-based coordinates'
    )
    parser.add_argument('pb_len',
        metavar='pb_len', 
        action='store', 
        type=int, 
        help = 'Length of the probe'
    )
    parser.add_argument('--output_path', 
        action='store',
        type=pathlib.Path, 
        default=None,
        dest='output_path', 
        help='Output path',
    )
    parser.add_argument('--min_primer_len',
        action='store',
        type=int, 
        default=17,
        dest='min_primer_len',
        help='Minimum primer length'
    )
    parser.add_argument('--max_primer_len',
        action='store',
        type=int, 
        default=22,
        dest='max_primer_len',
        help='Maximum primer length'
    )
    parser.add_argument('-d',
        dest='tm_diff',
        metavar='maximum_tm_diff',
        action='store',
        type=float,
        default=5.0,
        help='Maximum temperature difference between forward and reverse'
    )
    #Arguments for specificity checking
    parser.add_argument('--no_sens_spec',
        action='store_true',
        dest='sens_spec_flag',
        help='Flag to not check the putative probes for their specificity and sensitivity'
    )
    parser.add_argument('--blastdb',
        action='store',
        type=str,
        dest='blastdb',
        default='',
        help='Name of blastdb'
    )
    parser.add_argument('--blastdb_len',
        action='store',
        type=int,
        dest='blastdb_len',
        help='Length of blastdb'
    )
    args = parser.parse_args()

    if not args.output_path: 
        args.output_path = args.target_alignment_path.parent

    #Arguments
    #Target seq path, probe start, probe length, minimum primer length, max primer length, max allowable tm difference, sens_spec_flag, blastdb, blastdb length, target accessions path
    #Note conversion of pb_start to 0-based coordinate system
    return (
        args.target_alignment_path,
        args.output_path, 
        args.pb_start-1, 
        args.pb_len, 
        args.min_primer_len, 
        args.max_primer_len, 
        args.tm_diff, 
        args.sens_spec_flag, 
        args.blastdb, 
        args.blastdb_len, 
    )

def main(): 
    #Get arguments
    target_alignment_path, output_path, pb_start, pb_len, min_primer_len, max_primer_len, max_tm_diff, skip_check_flag, blastdb, blastdb_len = parse_args()
    #Process the alignment
    target_alignment = alignment(target_alignment_path)
    target_alignment.get_consensus()
    target_accessions = target_alignment.get_accessions()

    primer_gen = primerGenerator(
        target_alignment.consensus, 
        pb_start, pb_len, 
        min_primer_len, 
        max_primer_len, 
        max_tm_diff
    )

    #Generate primer pairs
    primer_gen.find_fw_primers()
    primer_gen.find_rev_primers()

    #Generate BLAST results
    if skip_check_flag is False: 
        primer_blast = blast(blastdb, blastdb_len)
        fw_blast_results = primer_blast.multi_blast(primer_gen.fw_primers)
        rev_blast_results = primer_blast.multi_blast(primer_gen.rev_primers)
        for fw_primer in primer_gen.fw_primers: 
            fw_primer.calculate_sensitivity(fw_blast_results[fw_primer.id], target_accessions)
            fw_primer.calculate_specificity(fw_blast_results[fw_primer.id], target_accessions, blastdb_len)
            fw_primer.calculate_score()
        for rev_primer in primer_gen.rev_primers: 
            rev_primer.calculate_sensitivity(rev_blast_results[rev_primer.id], target_accessions)
            rev_primer.calculate_specificity(rev_blast_results[rev_primer.id], target_accessions, blastdb_len)
            rev_primer.calculate_score()
        primer_blast.output(fw_blast_results, output_path, 'fw')
        primer_blast.output(rev_blast_results, output_path, 'rev')
    
    #Generate primer pairs
    primer_gen.find_primer_pairs()
    if skip_check_flag is False:
        #Generate primer pair data
        for primer_pair in primer_gen.primer_pairs: 
            fw_blast = fw_blast_results[primer_pair.fw_primer.id]
            rev_blast = rev_blast_results[primer_pair.rev_primer.id]
            primer_pair.calculate_sensitivity(fw_blast, rev_blast, target_accessions)
            primer_pair.calculate_specificity(fw_blast, rev_blast, target_accessions, blastdb_len)
            primer_pair.calculate_score() 

    #Output
    primer_gen.output(output_path)

if __name__ == '__main__':
    main()