from Bio import SeqIO

def modify_fasta_header(record, prefix):
    # Extract the necessary parts from the description
    parts = record.description.split(" ")
    uniref_id = parts[0]
    taxon = None
    
    # Find the Tax= part to extract the taxon name
    for part in parts:
        if part.startswith("Tax="):
            taxon = part.replace("Tax=", "").replace(" ", "_")
            break
    
    if taxon:
        # Create the new header
        new_id = f"{prefix}__{taxon}__{uniref_id.replace('UniRef50_', 'UR50_')}"
        record.id = new_id
        record.description = ""  # Clear the description to prevent duplication
    return record

def process_fasta(input_file, output_file, prefix):
    # Parse the FASTA file and modify the headers
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            modified_record = modify_fasta_header(record, prefix)
            SeqIO.write(modified_record, out_handle, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modify FASTA headers.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file.")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file.")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix to add to the new FASTA headers.")
    
    args = parser.parse_args()
    process_fasta(args.input, args.output, args.prefix)