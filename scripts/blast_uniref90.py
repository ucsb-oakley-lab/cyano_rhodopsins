import argparse
import requests
import time
from Bio import SeqIO
import os
import sys
import logging
import datetime

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

import datetime

def main(fasta_file, email, max_hits, out_file=None, force=False, log_file=None):
    # Configure logging
    if log_file:
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir)  # make logs/ if it doesn’t exist
        logging.basicConfig(
            filename=log_file,
            filemode="a",  # append mode
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        logging.info("==== New BLAST run started ====")
        logging.info(f"Query FASTA: {fasta_file}")
        logging.info(f"Email: {email}")
        logging.info(f"Max hits: {max_hits}")
        if out_file:
            logging.info(f"Output file (user-specified): {out_file}")
        if force:
            logging.info("Force overwrite: True")

    timestamp = datetime.datetime.now().strftime("%Y%m%d")

    for record in SeqIO.parse(fasta_file, "fasta"):
        if log_file:
            logging.info(f"Submitting BLAST for: {record.id}")
        print(f"Submitting BLAST for: {record.id}")
        try:
            job_id = submit_blast(str(record.seq), email)
            print(f"Job ID: {job_id}")
            if log_file:
                logging.info(f"Job ID: {job_id}")

            while True:
                status = check_status(job_id)
                print(f"[{time.strftime('%H:%M:%S')}] Status: {status}")
                if log_file:
                    logging.info(f"Status: {status}")
                if status == "FINISHED":
                    break
                elif status == "ERROR":
                    print(f"Job {job_id} failed.")
                    if log_file:
                        logging.error(f"Job {job_id} failed.")
                    return
                time.sleep(5)

            tab_result = get_tabular_results(job_id)
            if not tab_result:
                print(f"No hits found for {record.id}. Skipping.")
                if log_file:
                    logging.warning(f"No hits found for {record.id}")
                continue

            uniref_ids = parse_uniref_ids(tab_result, max_hits)
            print(f"Retrieved {len(uniref_ids)} UniRef IDs")
            if log_file:
                logging.info(f"Retrieved {len(uniref_ids)} UniRef IDs")

            fasta_seqs = fetch_fasta(uniref_ids)

            # Output filename logic
            if out_file and len(list(SeqIO.parse(fasta_file, "fasta"))) == 1:
                # If only one record in input and --out is given → use it directly
                output_file = out_file
            else:
                # Default: one file per input record
                safe_id = record.id.replace("|", "_").replace(":", "_")
                output_file = f"{safe_id}_top{max_hits}_uniref90_{timestamp}.fasta"

            # Safety check
            if os.path.exists(output_file) and not force:
                msg = f"ERROR: Output file {output_file} already exists. Use --force to overwrite."
                print(msg)
                if log_file:
                    logging.error(msg)
                sys.exit(1)

            write_fasta(fasta_seqs, output_file)
            print(f"Saved FASTA to: {output_file}")
            if log_file:
                logging.info(f"Saved FASTA to: {output_file}")

        except Exception as e:
            print(f"Error with {record.id}: {e}")
            if log_file:
                logging.error(f"Error with {record.id}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="BLAST protein sequences against UniRef90 and download top hits."
    )
    parser.add_argument("fasta_file", help="Input FASTA file with protein sequences.")
    parser.add_argument("email", help="Your email (required by EBI).")
    parser.add_argument("--max_hits", type=int, default=200, help="Number of top hits to retain (default: 200)")
    parser.add_argument("--out", help="Optional output FASTA filename")
    parser.add_argument("--force", action="store_true", help="Allow overwriting existing output file")
    parser.add_argument("--log", help="Optional log file to record run details (e.g., logs/blast_runs.log)")
    args = parser.parse_args()

    main(args.fasta_file, args.email, args.max_hits, args.out, args.force, args.log)


