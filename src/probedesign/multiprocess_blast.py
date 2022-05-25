import multiprocessing
import subprocess
import io
import tempfile
import pandas as pd

blastdb = '/data/db/SSU-nematodes.fasta'
blastdb_len = 16094
test_oligos = {
    '1630-17':'CAGAAGAGATCAACCTCGTTA',
    '1631-18':'AGAAGAGATCAACCTCGTT',
    '1640-19':'GAAGAGATCAACCTCGTT',
    '1641-20':'AAGAGATCAACCTCGTTA',
}
test_oligos2 = [
    [('1630-17','CAGAAGAGATCAACCTCGTTA'),
    ('1630-18','CAGAAGAGATCAACCTCGTTA'),],
    [('1631-19','CAGAAGAGATCAACCTCGTTA'),
    ('1631-20','CAGAAGAGATCAACCTCGTTA'),],
    [('1631-21','CAGAAGAGATCAACCTCGTTA'),
    ('1631-22','CAGAAGAGATCAACCTCGTTA'),],
    [('1631-23','CAGAAGAGATCAACCTCGTTA'),
    ('1631-24','CAGAAGAGATCAACCTCGTTA'),],
]


def blast_all(oligos): 
    #Generate the oligo temporary file
    fasta = tempfile.NamedTemporaryFile(delete=True)
    #print(fasta)
    print(oligos)
    print(type(oligos))
    for oligo in oligos:
        print(f">{oligo[0]}\n{oligo[1]}\n")
        fasta.write(f">{oligo[0]}\n{oligo[1]}\n".encode())
        #fasta.write(f">{oligo}\n{oligos[oligo]}\n".encode())
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
    #print(data)
    #Split the data
    list_oligo_ids = set(data['qacc'])
    #print(f"qacc set:\n{list_oligo_ids}")
    blast_results = dict()
    for oligo_id in list_oligo_ids: 
        blast_results[str(oligo_id)]=data.loc[data['qacc']==oligo_id]
    return blast_results

def multiblast(oligos): 
    #print(oligos)
    pool = multiprocessing.Pool(4)
    result = pool.map(blast_all, oligos)
    print(result)
    for item in result: 
        print(item)

if __name__ == "__main__":
    # multiblast(test_oligos)
    multiblast(test_oligos2)
