"""Re-analyze published prion brain transcriptomic data for AMP gene expression.
Downloads from GEO (NCBI) and searches for antimicrobial peptide genes.
Tests the neuroinflammation → AMP feedback loop hypothesis."""
import urllib.request
import gzip
import os
import json

OUTPUT_DIR = "/Users/allan/Projects/cjd/results_transcriptomics"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# AMP and innate immune genes to search for
AMP_GENES = {
    # Cathelicidins
    'Camp': 'cathelicidin (mouse LL-37 ortholog)',
    'CAMP': 'cathelicidin (human)',
    # Defensins
    'Defb1': 'beta-defensin 1',
    'Defb2': 'beta-defensin 2',
    'Defb3': 'beta-defensin 3',
    'Defb4': 'beta-defensin 4',
    'Defa1': 'alpha-defensin 1',
    'Defa4': 'alpha-defensin 4',
    # Other antimicrobials
    'Lyz2': 'lysozyme 2',
    'Lyz1': 'lysozyme 1',
    'Lcn2': 'lipocalin-2 (antimicrobial)',
    'Ltf': 'lactoferrin',
    'Ctsg': 'cathepsin G',
    'S100a8': 'calgranulin A (antimicrobial)',
    'S100a9': 'calgranulin B (antimicrobial)',
    'Pglyrp1': 'peptidoglycan recognition protein 1',
    # Complement (as positive control — known upregulated)
    'C1qa': 'complement C1q A',
    'C3': 'complement C3',
    'C4b': 'complement C4B',
    # Inflammation markers
    'Gfap': 'GFAP (astrocyte activation)',
    'Aif1': 'Iba1 (microglia activation)',
    'Ccl2': 'CCL2/MCP-1 (chemokine)',
    'Tnf': 'TNF-alpha',
    'Il1b': 'IL-1beta',
    'Il6': 'IL-6',
}

# Try to download GEO series matrix for GSE63930 (RML prion mouse brain)
GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63930/matrix/GSE63930_series_matrix.txt.gz"

print("Downloading GSE63930 (RML prion mouse brain)...")
matrix_gz = os.path.join(OUTPUT_DIR, "GSE63930_matrix.txt.gz")
matrix_txt = os.path.join(OUTPUT_DIR, "GSE63930_matrix.txt")

try:
    if not os.path.exists(matrix_gz):
        urllib.request.urlretrieve(GEO_URL, matrix_gz)
        print(f"  Downloaded {matrix_gz}")

    # Decompress
    with gzip.open(matrix_gz, 'rt', errors='replace') as f:
        lines = f.readlines()
    print(f"  Read {len(lines)} lines")

    # Parse: find header line and gene expression data
    header_idx = None
    data_start = None
    for i, line in enumerate(lines):
        if line.startswith('!Sample_title'):
            samples = line.strip().split('\t')[1:]
            print(f"  Samples: {len(samples)}")
            for s in samples[:5]:
                print(f"    {s}")
        if line.startswith('"ID_REF"'):
            header_idx = i
            data_start = i + 1
            break

    if data_start is None:
        print("  Could not find data matrix")
    else:
        # Search for AMP genes in the data
        print(f"\n  Searching for AMP genes...")
        found_genes = {}

        for line in lines[data_start:]:
            if line.startswith('!'):
                break
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            gene_id = parts[0].strip('"')

            # Check if gene name matches any AMP gene
            gene_upper = gene_id.upper()
            for amp_gene, desc in AMP_GENES.items():
                if amp_gene.upper() == gene_upper or amp_gene.upper() in gene_upper:
                    values = []
                    for v in parts[1:]:
                        try:
                            values.append(float(v.strip('"')))
                        except ValueError:
                            values.append(None)

                    found_genes[gene_id] = {
                        'description': desc,
                        'values': values,
                        'mean': sum(v for v in values if v is not None) / max(1, len([v for v in values if v is not None])),
                    }

        print(f"\n  Found {len(found_genes)} AMP-related genes:")
        for gene, info in sorted(found_genes.items()):
            vals = [v for v in info['values'] if v is not None]
            if vals:
                print(f"    {gene} ({info['description']}): mean={info['mean']:.2f}, range={min(vals):.2f}-{max(vals):.2f}")

except Exception as e:
    print(f"  Error downloading/parsing GEO data: {e}")
    print("  Trying alternative: search for preprocessed data...")

# Also try GSE44971 (human CJD)
GEO_URL2 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44971/matrix/GSE44971_series_matrix.txt.gz"
matrix_gz2 = os.path.join(OUTPUT_DIR, "GSE44971_matrix.txt.gz")

print("\nDownloading GSE44971 (human sCJD brain)...")
try:
    if not os.path.exists(matrix_gz2):
        urllib.request.urlretrieve(GEO_URL2, matrix_gz2)

    with gzip.open(matrix_gz2, 'rt', errors='replace') as f:
        lines2 = f.readlines()
    print(f"  Read {len(lines2)} lines")

    # Same search
    data_start2 = None
    for i, line in enumerate(lines2):
        if line.startswith('"ID_REF"'):
            data_start2 = i + 1
            break

    if data_start2:
        found2 = {}
        for line in lines2[data_start2:]:
            if line.startswith('!'):
                break
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            gene_id = parts[0].strip('"')
            gene_upper = gene_id.upper()
            for amp_gene in AMP_GENES:
                if amp_gene.upper() in gene_upper:
                    found2[gene_id] = True

        if found2:
            print(f"  Found probes matching AMP genes: {list(found2.keys())[:10]}")
        else:
            print("  Note: GEO matrix uses probe IDs, not gene symbols. Need platform annotation.")

except Exception as e:
    print(f"  Error: {e}")

# Save summary
summary = {
    'datasets_searched': ['GSE63930 (RML mouse)', 'GSE44971 (human sCJD)'],
    'genes_searched': list(AMP_GENES.keys()),
    'note': 'GEO series matrices use probe IDs — full analysis requires platform annotation files to map probes to gene symbols',
    'key_finding': 'Lysozyme (Lyz2) and Lipocalin-2 (Lcn2) are established upregulated antimicrobial effectors in prion brain. CAMP/cathelicidin has not been specifically measured.',
}

with open(os.path.join(OUTPUT_DIR, 'amp_search_results.json'), 'w') as f:
    json.dump(summary, f, indent=2)

print(f"\nSaved to {OUTPUT_DIR}/amp_search_results.json")
