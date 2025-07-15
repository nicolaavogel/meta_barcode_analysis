import os
from time import sleep
from Bio import Entrez, SeqIO

# ---- Configuration ----
EMAIL = "n.alexandra.vogel@gmail.com"
API_KEY = None
FASTA_INPUT = "/projects/caeg/people/bfj994/barcodeMiner/DBs/viridiplantae_chloroplast/fasta/viri_full_chlo.fa"  # Multi-FASTA input
OUTPUT_DIR = "trnl_extracted"
TRNL_NAME_MATCH = ["trnL"]  # Generalized filter

# ---- Entrez Setup ----
Entrez.email = EMAIL
if API_KEY:
    Entrez.api_key = API_KEY

# ---- Create output folder ----
os.makedirs(OUTPUT_DIR, exist_ok=True)
not_found_file = os.path.join(OUTPUT_DIR, "trnl_not_found.txt")
not_found_list = []

# ---- Process each genome ----
for i, record in enumerate(SeqIO.parse(FASTA_INPUT, "fasta"), 1):
    acc = record.id.split('.')[0]  # Strip version suffix
    print(f"[{i}] Fetching GenBank for: {acc}")
    
    try:
        with Entrez.efetch(db="nuccore", id=acc, rettype="gbwithparts", retmode="text") as handle:
            gb_record = SeqIO.read(handle, "genbank")
        sleep(0.5)

        trnl_found = False
        for feature in gb_record.features:
            if feature.type in ["gene", "tRNA"]:
                gene_names = feature.qualifiers.get("gene", []) + feature.qualifiers.get("product", [])
                if any("trnl" in name.lower() for name in gene_names):
                    trnl_seq = feature.extract(gb_record.seq)
                    header = f"{record.id}_trnL"
                    out_file = os.path.join(OUTPUT_DIR, f"{acc}_trnL.fasta")
                    with open(out_file, "w") as out:
                        out.write(f">{header}\n{trnl_seq}\n")
                    trnl_found = True
                    print(f"✅ trnL found and written for {acc}")
                    break

        if not trnl_found:
            print(f"⚠️  trnL not found in {acc}")
            not_found_list.append(acc)

    except Exception as e:
        print(f"❌ Error fetching or processing {acc}: {e}")
        not_found_list.append(f"{acc} - ERROR: {e}")
        continue

# ---- Write missing trnL list ----
if not_found_list:
    with open(not_found_file, "w") as nf:
        nf.write("\n".join(not_found_list))
    print(f"📄 Wrote list of genomes without trnL to: {not_found_file}")

print("🎉 Done extracting trnL regions.")
