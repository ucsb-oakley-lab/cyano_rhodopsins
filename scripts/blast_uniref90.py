import argparse
import requests
import time
from Bio import SeqIO
import os

def submit_blast(sequence, email, db='uniref90', program='blastp'):
    url = 'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run'
    clean_seq = sequence.replace("\n", "").replace("*", "").strip()
    data = {
        'email': email,
        'sequence': clean_seq,
        'stype': 'protein',
        'database': db,
        'program': program
    }
    response = requests.post(url, data=data)
    response.raise_for_status()
    return response.text.strip()

def check_status(job_id):
    url = f'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}'
    response = requests.get(url)
    response.raise_for_status()
    return response.text.strip()

def get_tabular_results(job_id):
    url = f'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/tab'
    response = requests.get(url)
    if response.status_code == 400:
        print(f"Job {job_id} finished but returned no results in tab format (likely no hits).")
        return None
    response.raise_for_status()
    return response.text

def parse_uniref_ids(tab_text, max_hits):
    ids = []
    for line in tab_text.strip().splitlines():
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        subject_id = fields[1]
        if subject_id.startswith("UniRef90_"):
            ids.append(subject_id)
        if len(ids) >= max_hits:
            break
    return ids

def fetch_fasta(uniref_ids):
    sequences = []
    for uid in uniref_ids:
        url = f"https://rest.uniprot.org/uniref/{uid}.fasta"
        try:
            response = requests.get(url)
            response.raise_for_status()
            sequences.append(response.text)
        except requests.exceptions.HTTPError:
            print(f"Warning: could not fetch {uid}")
        time.sleep(0.3)  # rate limit to avoid server errors
    return sequences

def write_fasta(sequences, out_file):
    with open(out_file, 'w') as f:
        for seq in sequences:
            f.write(seq)

def main(fasta_file, email, max_hits):
    for record in SeqIO.parse(fasta_file, "fasta"):
        print(f"Submitting BLAST for: {record.id}")
        try:
            job_id = submit_blast(str(record.seq), email)
            print(f"Job ID: {job_id}")

            while True:
                status = check_status(job_id)
                print(f"[{time.strftime('%H:%M:%S')}] Status: {status}")
                if status == "FINISHED":
                    break
                elif status == "ERROR":
                    print(f"Job {job_id} failed.")
                    return
                time.sleep(5)

            tab_result = get_tabular_results(job_id)
            if not tab_result:
                print(f"No hits found for {record.id}. Skipping.")
                continue

            uniref_ids = parse_uniref_ids(tab_result, max_hits)
            print(f"Retrieved {len(uniref_ids)} UniRef IDs")

            fasta_seqs = fetch_fasta(uniref_ids)
            output_file = f"{record.id}_top{max_hits}_uniref90.fasta"
            write_fasta(fasta_seqs, output_file)
            print(f"Saved FASTA to: {output_file}")

        except Exception as e:
            print(f"Error with {record.id}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BLAST protein sequences against UniRef90 and download top hits.")
    parser.add_argument("fasta_file", help="Input FASTA file with protein sequences.")
    parser.add_argument("email", help="Your email (required by EBI).")
    parser.add_argument("--max_hits", type=int, default=200, help="Number of top hits to retain (default: 200)")
    args = parser.parse_args()

    main(args.fasta_file, args.email, args.max_hits)

