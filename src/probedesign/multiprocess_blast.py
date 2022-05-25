import multiprocessing
import subprocess
import io
import tempfile
import pandas as pd

blastdb = '/data/db/SSU-nematodes.fasta'
blastdb_len = 16094
oligos = {
    '1630-17':'CAGAAGAGATCAACCTCGTTA',
    '1631-18':'AGAAGAGATCAACCTCGTT',
    '1640-19':'GAAGAGATCAACCTCGTT',
    '1641-20':'AAGAGATCAACCTCGTTA',
}


def blast_all(oligos): 
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
        blastdb,
        "-num_alignments",
        str(blastdb_len),
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
    data = pd.read_csv(output, sep=',', header=None, names=headers)
    fasta.close()
    #Split the data
    list_oligo_ids = set(data['qacc'])
    #print(f"qacc set:\n{list_oligo_ids}")
    blast_results = dict()
    for oligo_id in list_oligo_ids: 
        blast_results[str(oligo_id)]=data.loc[data['qacc']==oligo_id]
    return blast_results

def multiblast(oligos): 
    pool = multiprocessing.pool(4)
    result = pool.map(blast_all, oligos)
    print(len(result))

if __name__ == "__main__":
    multiblast(oligos)