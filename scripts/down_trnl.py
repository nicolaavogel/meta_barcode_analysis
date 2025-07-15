import os
from time import sleep
from Bio import Entrez, SeqIO

# ---- Configuration ----
EMAIL = "n.alexandra.vogel@gmail.com"  # REQUIRED: Set your email for Entrez
API_KEY = None  # Optional: add your NCBI API key here
SEARCH_TERM = "chloroplast[filter] AND Viridiplantae[Organism] AND complete genome AND srcdb_refseq[PROP]"
MAX_RECORDS = 5000
OUTPUT_DIR = "viridiplantae_chloroplast"
TRNL_NAME_MATCH = ["trnL", "trnL-UAA", "trnL(CAA)"]

# ---- Entrez Setup ----
Entrez.email = EMAIL
if API_KEY:
    Entrez.api_key = API_KEY

# ---- Create folders ----
fasta_dir = os.path.join(OUTPUT_DIR, "fasta")
gb_dir = os.path.join(OUTPUT_DIR, "gb")
trnl_dir = os.path.join(OUTPUT_DIR, "trnl")

os.makedirs(fasta_dir, exist_ok=True)
os.makedirs(gb_dir, exist_ok=True)
os.makedirs(trnl_dir, exist_ok=True)

START_AT = 5000  # Adjust this to however many you've already processed
CHUNK_SIZE = 5000  # Number of records to fetch at a time

# ---- Step 1: Search NCBI ----
print(f"🔍 Searching NCBI for chloroplast genomes starting at {START_AT}...")
search_handle = Entrez.esearch(
    db="nuccore",
    term=SEARCH_TERM,
    retstart=START_AT,
    retmax=CHUNK_SIZE
)
search_results = Entrez.read(search_handle)
search_handle.close()
accessions = search_results["IdList"]
print(f"✅ Found {len(accessions)} additional genomes.")
# ---- Step 2: Download and parse each record ----
for i, acc in enumerate(accessions, 1):
    print(f"[{i}/{len(accessions)}] Processing accession: {acc}")
    try:
        # Fetch GenBank file
        gb_file = os.path.join(gb_dir, f"{acc}.gb")
        if not os.path.exists(gb_file):
            with Entrez.efetch(db="nuccore", id=acc, rettype="gbwithparts", retmode="text") as handle:
                gb_data = handle.read()
            with open(gb_file, "w") as out:
                out.write(gb_data)
            sleep(0.5)

        # Fetch FASTA file
        fasta_file = os.path.join(fasta_dir, f"{acc}.fasta")
        if not os.path.exists(fasta_file):
            with Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text") as handle:
                fasta_data = handle.read()
            with open(fasta_file, "w") as out:
                out.write(fasta_data)
            sleep(0.5)

        # Parse FASTA header for use in trnL output
        original_header = None
        if os.path.exists(fasta_file):
            with open(fasta_file) as f:
                first_line = f.readline().strip()
                if first_line.startswith(">"):
                    original_header = first_line[1:].replace(" ", "_")

        # Parse GenBank and extract trnL
        trnl_found = False
        for record in SeqIO.parse(gb_file, "genbank"):
            for feature in record.features:
                if feature.type in ["gene", "tRNA"]:
                    gene_names = feature.qualifiers.get("gene", []) + feature.qualifiers.get("product", [])
                    if any(trnl.lower() in name.lower() for name in gene_names for trnl in TRNL_NAME_MATCH):
                        trnl_seq = feature.extract(record.seq)
                        if 300 <= len(trnl_seq) <= 600:
                            header = original_header if original_header else f"{acc}_trnL"
                            out_file = os.path.join(trnl_dir, f"{acc}_trnL.fasta")
                            with open(out_file, "w") as out:
                                out.write(f">{header}\n{trnl_seq}\n")
                            trnl_found = True
                            break
            if trnl_found:
                break

        if not trnl_found:
            print(f"⚠️  trnL gene not found or too short/long in {acc}")

    except Exception as e:
        print(f"❌ Error processing {acc}: {e}")
        continue

print("🎉 Done downloading and extracting trnL genes.")
