import requests
import time
from Bio import SeqIO
import os

# EBI BLASTP submission
def submit_blast(sequence, email, db='uniref100', program='blastp'):
    url = 'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run'
    data = {
        'email': email,
        'sequence': sequence,
        'database': db,
        'program': program,
    }
    response = requests.post(url, data=data)
    response.raise_for_status()
    return response.text.strip()

# Check job status
def check_status(job_id):
    url = f'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}'
    response = requests.get(url)
    response.raise_for_status()
    return response.text.strip()

# Retrieve result
def get_results(job_id, result_type='out'):
    url = f'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/{result_type}'
    response = requests.get(url)
    response.raise_for_status()
    return response.text

# Run BLAST for all sequences in a FASTA file
def blast_fasta_file(fasta_path, email, out_dir='blast_results', 
result_type='out'):
    os.makedirs(out_dir, exist_ok=True)
    for record in SeqIO.parse(fasta_path, "fasta"):
        print(f"Submitting BLAST for: {record.id}")
        try:
            job_id = submit_blast(str(record.seq), email)
            print(f"Job ID: {job_id}")
            
            # Poll until done
            while True:
                status = check_status(job_id)
                if status == 'FINISHED':
                    break
                elif status == 'ERROR':
                    print(f"Job {job_id} failed.")
                    break
                time.sleep(5)
            
            # Get and save result
            result = get_results(job_id, result_type=result_type)
            out_path = os.path.join(out_dir, f"{record.id}.{result_type}")
            with open(out_path, 'w') as f:
                f.write(result)
            print(f"Saved result to {out_path}")
        
        except Exception as e:
            print(f"Error with {record.id}: {e}")

# Example usage
if __name__ == "__main__":
    fasta_file = "your_sequences.fasta"  # <- Change to your FASTA path
    user_email = "your_email@domain.com"  # <- Must be valid
    blast_fasta_file(fasta_file, user_email)

