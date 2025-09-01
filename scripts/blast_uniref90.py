import argparse
import requests
from Bio import SeqIO
import os
import time
import sys
import io
import logging
from Bio import Entrez
import xml.etree.ElementTree as ET

def fetch_sequence_by_accession(accession, email):
    Entrez.email = email
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
    seq = handle.read()
    handle.close()
    if not seq.strip() or not seq.startswith(">"):
        raise ValueError(f"Could not fetch sequence for accession {accession}")
    return seq

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

def get_xml_results(job_id):
    url = f'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/xml'
    response = requests.get(url)
    response.raise_for_status()
    return response.text

def parse_uniref_ids_from_xml(xml_text, max_hits):
    import xml.etree.ElementTree as ET
    ids = []
    root = ET.fromstring(xml_text)
    # The <hit> elements are inside <SequenceSimilaritySearchResult>
    for hit in root.findall(".//{http://www.ebi.ac.uk/schema}hit"):
        uid = hit.attrib.get("id")
        if uid and uid.startswith("UniRef90_"):
            ids.append(uid)
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

def main(accession, email, max_hits, out_file=None, force=False, log_file=None, job_id=None):
    output_file = out_file if out_file else f"{accession}_top{max_hits}_uniref90.fasta"
    # Check for existing output file before starting any BLAST/search
    if os.path.exists(output_file) and not force:
        msg = f"ERROR: Output file {output_file} already exists. Use --force to overwrite."
        print(msg)
        if log_file:
            logging.error(msg)
        sys.exit(1)

    if job_id:
        print(f"Using existing BLAST job ID: {job_id}")
        if log_file:
            logging.info(f"Using existing BLAST job ID: {job_id}")
        # Wait for job to finish
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
        xml_result = get_xml_results(job_id)
        uniref_ids = parse_uniref_ids_from_xml(xml_result, max_hits)
        print(f"Retrieved {len(uniref_ids)} UniRef IDs")
        if log_file:
            logging.info(f"Retrieved {len(uniref_ids)} UniRef IDs")
        fasta_seqs = fetch_fasta(uniref_ids)
        write_fasta(fasta_seqs, output_file)
        print(f"Saved FASTA to: {output_file}")
        if log_file:
            logging.info(f"Saved FASTA to: {output_file}")
        return

    # Download sequence from GenBank
    print(f"Fetching sequence for accession: {accession}")
    if log_file:
        logging.info(f"Fetching sequence for accession: {accession}")
    try:
        fasta_seq = fetch_sequence_by_accession(accession, email)
    except Exception as e:
        print(f"Error fetching sequence: {e}")
        if log_file:
            logging.error(f"Error fetching sequence: {e}")
        sys.exit(1)

    record = next(SeqIO.parse(io.StringIO(fasta_seq), "fasta"))
    sequence = str(record.seq)

    print(f"Submitting BLAST for: {accession}")
    if log_file:
        logging.info(f"Submitting BLAST for: {accession}")
    try:
        job_id = submit_blast(sequence, email)
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

        xml_result = get_xml_results(job_id)
        uniref_ids = parse_uniref_ids_from_xml(xml_result, max_hits)
        print(f"Retrieved {len(uniref_ids)} UniRef IDs")
        if log_file:
            logging.info(f"Retrieved {len(uniref_ids)} UniRef IDs")

        fasta_seqs = fetch_fasta(uniref_ids)
        write_fasta(fasta_seqs, output_file)
        print(f"Saved FASTA to: {output_file}")
        if log_file:
            logging.info(f"Saved FASTA to: {output_file}")

    except Exception as e:
        print(f"Error with {accession}: {e}")
        if log_file:
            logging.error(f"Error with {accession}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="BLAST a protein accession against UniRef90 and download top hits, or use an existing EBI BLAST job ID."
    )
    parser.add_argument("accession", nargs="?", help="GenBank protein accession number.")
    parser.add_argument("email", nargs="?", help="Your email (required by EBI).")
    parser.add_argument("--max_hits", type=int, default=200, help="Number of top hits to retain (default: 200)")
    parser.add_argument("--out", help="Optional output FASTA filename")
    parser.add_argument("--force", action="store_true", help="Allow overwriting existing output file")
    parser.add_argument("--log", help="Optional log file to record run details (e.g., logs/blast_runs.log)")
    parser.add_argument("--jobid", help="Existing EBI BLAST job ID to fetch results (skip BLAST submission)")
    args = parser.parse_args()

    main(args.accession, args.email, args.max_hits, args.out, args.force, args.log, args.jobid)