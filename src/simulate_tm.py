#Standard libraries
import pathlib
import argparse
import logging
import random
import csv
#Third-party libraries
import pandas as pd
import numpy as np
import primer3
#Local modules
from oligogenerator import ProbeGenerator
from consensus import Alignment

def parse_args(): 
    args = argparse.ArgumentParser('Script to compare Thermofisher Tm calculation to Primer3 calculation')
    args.add_argument(
        '-o',
        '--output',
        dest='output_path',
        type=pathlib.Path,
        default=None, 
        help='Output path for the simulation file'
    )
    args.add_argument(
        '--min_len',
        type=float,
        default=17.0,
        help='Minimum length for oligo'
    )
    args.add_argument(
        '--max_len',
        type=float,
        default=23.0,
        help='Maximum length for simulated oligo',
    )

    args = parse_args()

    if args.output_path is None: 
        args.output_path = pathlib.Path.cwd()
    
    return args

def main(): 
    args = parse_args()
    
    log_path = args.output_path.joinpath('probesearch.log')
    sim_path = args.output_path.joinpath(f'sim_path.csv')    

    #Set up logger
    logging.basicConfig(
        encoding='utf-8',
        level=logging.INFO,
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler()
        ],
        format='%(asctime)s:%(levelname)s: %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S',
    )

    #Process the alignment
    logging.info(f'simulate_tm.py - compare Thermofisher Tm calculation to Primer3.')
    #Generate Probes
    logging.info(f'Generating probes...')

    SIM_LENGTH = 10000
    NUCL = ['A', 'T', 'C', 'G']

    simulated_seq = ''.join(random.choice(NUCL) for i in range(SIM_LENGTH))

    print(simulated_seq)

    pb_gen = ProbeGenerator(
        simulated_seq,
        0,
        None,
        args.min_len,
        args.max_len)
    pb_gen.get_probes()

    logging.info(f'{len(pb_gen.probes)} probes generated...')

    data = []

    logging.info(f'Calculating primer3 tms...')

    for probe in pb_gen.probes:
        data.append(
            (
                probe.seq,
                probe.tm, 
                primer3.calcTm(probe.seq), 
                probe.tm - primer3.calcTm(probe.seq)
            )
        )

    logging.info(f'Outputting data...')

    with open(sim_path, 'w', newline='') as csv_file: 
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(
            (
                'seq',
                'thermo_tm',
                'primer3_tm',
                'dTm'
            )
        )
        csv_writer.writerows(data)

    logging.info(f'Done!')

if __name__ == '__main__': 
    main()