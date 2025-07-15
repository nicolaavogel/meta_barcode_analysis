import os
import time
import argparse
from Bio import Entrez, SeqIO

# ---- Setup ----
Entrez.email = "n.alexandra.vogel@gmail.com"

# Default values
TRNL_NAME_MATCH = ["trnL", "trnL-UAA", "trnL(CAA)"]
DEFAULT_SEARCH_TERM = "chloroplast[filter] AND Viridiplantae[Organism] AND complete genome AND srcdb_refseq[PROP]"

# ---- Parse arguments ----
parser = argparse.ArgumentParser(description="Download all matching chloroplast genomes and extract trnL genes.")
parser.add_argument("--search", type=str, default=DEFAULT_SEARCH_TERM, help="Entrez search string")
parser.add_argument("--output_dir", type=str, default="viridiplantae_chloroplast", help="Output folder")
parser.add_argument("--max_records", type=int, default=100000, help="Maximum number of records to fetch")
args = parser.parse_args()

# ---- Output dirs ----
fasta_dir = os.path.join(args.output_dir, "fasta")
gb_dir = os.path.join(args.output_dir, "gb")
trnl_dir = os.path.join(args.output_dir, "trnl")

for path in [fasta_dir, gb_dir, trnl_dir]:
    os.makedirs(path, exist_ok=True)

# ---- Step 1: Search NCBI ----
print(f"🔍 Searching NCBI for: {args.search}")
search_handle = Entrez.esearch(
    db="nuccore",
    term=args.search,
    usehistory="y",
    retmax=args.max_records
)
search_results = Entrez.read(search_handle)
search_handle.close()

id_list = search_results["IdList"]
print(f"✅ Found {len(id_list)} genomes.")

# ---- Step 2: Process each accession ----
for i, acc in enumerate(id_list, 1):
    print(f"[{i}/{len(id_list)}] Accession: {acc}")
    try:
        # Download GenBank
        gb_file = os.path.join(gb_dir, f"{acc}.gb")
        if not os.path.exists(gb_file):
            with Entrez.efetch(db="nuccore", id=acc, rettype="gbwithparts", retmode="text") as handle:
                with open(gb_file, "w") as out:
                    out.write(handle.read())
            time.sleep(0.5)

        # Download FASTA
        fasta_file = os.path.join(fasta_dir, f"{acc}.fasta")
        if not os.path.exists(fasta_file):
            with Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text") as handle:
                with open(fasta_file, "w") as out:
                    out.write(handle.read())
            time.sleep(0.5)

        # Parse FASTA header
        header = acc
        with open(fasta_file) as f:
            line = f.readline()
            if line.startswith(">"):
                header = line[1:].strip().replace(" ", "_")

        # Parse GenBank for trnL
        trnl_found = False
        for record in SeqIO.parse(gb_file, "genbank"):
            for feature in record.features:
                if feature.type in ["gene", "tRNA"]:
                    names = feature.qualifiers.get("gene", []) + feature.qualifiers.get("product", [])
                    if any(trnl.lower() in name.lower() for name in names for trnl in TRNL_NAME_MATCH):
                        trnl_seq = feature.extract(record.seq)
                        if 300 <= len(trnl_seq) <= 600:
                            out_file = os.path.join(trnl_dir, f"{acc}_trnL.fasta")
                            with open(out_file, "w") as out:
                                out.write(f">{header}_trnL\n{trnl_seq}\n")
                            trnl_found = True
                            break
            if trnl_found:
                break
        if not trnl_found:
            print(f"⚠️ trnL not found or invalid length in {acc}")

    except Exception as e:
        print(f"❌ Error processing {acc}: {e}")
        continue

print("🎉 Done! All genomes downloaded and trnL genes extracted.")
