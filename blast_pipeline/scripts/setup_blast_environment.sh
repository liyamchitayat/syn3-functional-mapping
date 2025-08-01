#!/bin/bash
# BLAST Pipeline Setup Script
# Sets up environment for comprehensive homology analysis

set -e  # Exit on any error

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="${PIPELINE_DIR}/logs"
DB_DIR="${PIPELINE_DIR}/databases"
SCRIPTS_DIR="${PIPELINE_DIR}/scripts"

echo "=== BLAST PIPELINE SETUP ==="
echo "Pipeline directory: ${PIPELINE_DIR}"
echo "Log directory: ${LOG_DIR}"
echo "Database directory: ${DB_DIR}"
echo ""

# Create log file
LOG_FILE="${LOG_DIR}/setup_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "$(date): Starting BLAST pipeline setup"

# Check if BLAST+ is installed
echo "Checking BLAST+ installation..."
if command -v blastp >/dev/null 2>&1; then
    BLAST_VERSION=$(blastp -version | head -1)
    echo "✓ BLAST+ found: $BLAST_VERSION"
else
    echo "❌ BLAST+ not found. Installing..."
    
    # Check system and install BLAST+
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS
        if command -v brew >/dev/null 2>&1; then
            echo "Installing BLAST+ via Homebrew..."
            brew install blast
        else
            echo "Homebrew not found. Please install BLAST+ manually:"
            echo "1. Download from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
            echo "2. Or install Homebrew first: /bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
            exit 1
        fi
    elif [[ "$OSTYPE" == "linux"* ]]; then
        # Linux
        if command -v apt-get >/dev/null 2>&1; then
            echo "Installing BLAST+ via apt..."
            sudo apt-get update
            sudo apt-get install -y ncbi-blast+
        elif command -v yum >/dev/null 2>&1; then
            echo "Installing BLAST+ via yum..."
            sudo yum install -y ncbi-blast+
        elif command -v conda >/dev/null 2>&1; then
            echo "Installing BLAST+ via conda..."
            conda install -c bioconda blast
        else
            echo "No suitable package manager found. Please install BLAST+ manually."
            exit 1
        fi
    else
        echo "Unsupported OS: $OSTYPE"
        echo "Please install BLAST+ manually from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download"
        exit 1
    fi
fi

# Verify BLAST+ installation
if command -v blastp >/dev/null 2>&1; then
    BLAST_VERSION=$(blastp -version | head -1)
    echo "✓ BLAST+ verified: $BLAST_VERSION"
else
    echo "❌ BLAST+ installation failed"
    exit 1
fi

# Check Python dependencies
echo ""
echo "Checking Python dependencies..."
python3 -c "import sys; print(f'Python version: {sys.version}')"

# Check required Python packages
PYTHON_PACKAGES=("requests" "json" "csv" "argparse" "concurrent.futures" "subprocess")
for package in "${PYTHON_PACKAGES[@]}"; do
    if python3 -c "import $package" 2>/dev/null; then
        echo "✓ Python package '$package' available"
    else
        echo "❌ Python package '$package' missing"
        if [[ "$package" == "requests" ]]; then
            echo "Installing requests..."
            python3 -m pip install requests --user
        fi
    fi
done

# Create database directory structure
echo ""
echo "Setting up database directories..."
mkdir -p "${DB_DIR}/nr_bacteria"
mkdir -p "${DB_DIR}/custom_bacteria"
mkdir -p "${DB_DIR}/logs"

echo "✓ Database directories created"

# Download bacterial taxonomy information
echo ""
echo "Setting up bacterial taxonomy filters..."
cat > "${DB_DIR}/bacterial_taxonomy_setup.py" << 'EOF'
#!/usr/bin/env python3
"""
Setup bacterial taxonomy information for BLAST filtering
"""

import requests
import json
import time

def get_bacterial_taxids():
    """Get taxonomy IDs for major bacterial groups"""
    
    # Major bacterial phyla and important groups
    bacterial_groups = [
        "Proteobacteria",
        "Firmicutes", 
        "Actinobacteria",
        "Bacteroidetes",
        "Cyanobacteria",
        "Spirochaetes",
        "Chlamydiae",
        "Deinococcus-Thermus",
        "Aquifex",
        "Thermotoga"
    ]
    
    print("Fetching bacterial taxonomy IDs...")
    taxids = {}
    
    for group in bacterial_groups:
        try:
            # Search for taxonomy ID
            search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            params = {
                'db': 'taxonomy',
                'term': f'"{group}"[Scientific Name]',
                'retmode': 'json'
            }
            
            response = requests.get(search_url, params=params)
            data = response.json()
            
            if 'esearchresult' in data and 'idlist' in data['esearchresult']:
                if data['esearchresult']['idlist']:
                    taxid = data['esearchresult']['idlist'][0]
                    taxids[group] = taxid
                    print(f"✓ {group}: {taxid}")
                else:
                    print(f"❌ {group}: No taxonomy ID found")
            
            time.sleep(0.5)  # Rate limiting
            
        except Exception as e:
            print(f"❌ {group}: Error - {e}")
    
    # Save taxonomy information
    with open('bacterial_taxids.json', 'w') as f:
        json.dump(taxids, f, indent=2)
    
    print(f"\nSaved {len(taxids)} bacterial taxonomy IDs")
    return taxids

if __name__ == "__main__":
    taxids = get_bacterial_taxids()
EOF

cd "${DB_DIR}"
python3 bacterial_taxonomy_setup.py

# Create BLAST configuration
echo ""
echo "Creating BLAST configuration files..."

cat > "${PIPELINE_DIR}/blast_config.json" << EOF
{
    "pipeline_info": {
        "name": "Syn3A Bacterial Homology Pipeline",
        "version": "1.0",
        "created": "$(date -Iseconds)",
        "description": "Comprehensive BLAST pipeline for syn3A protein homology analysis"
    },
    "blast_parameters": {
        "program": "blastp",
        "database": "nr",
        "evalue": 0.01,
        "max_target_seqs": 100,
        "outfmt": "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames stitle",
        "word_size": 6,
        "matrix": "BLOSUM62",
        "gap_open": 11,
        "gap_extend": 1,
        "seg": "yes"
    },
    "filtering": {
        "exclude_organisms": ["mycoplasma"],
        "min_identity": 30.0,
        "min_coverage": 50.0,
        "max_evalue": 0.01
    },
    "directories": {
        "pipeline": "${PIPELINE_DIR}",
        "databases": "${DB_DIR}",
        "queries": "${PIPELINE_DIR}/queries",
        "results": "${PIPELINE_DIR}/results",
        "logs": "${LOG_DIR}",
        "scripts": "${SCRIPTS_DIR}"
    },
    "resources": {
        "max_parallel_jobs": 4,
        "memory_limit_gb": 8,
        "temp_dir": "/tmp/blast_pipeline"
    }
}
EOF

echo "✓ BLAST configuration created"

# Create query preparation script
echo ""
echo "Creating query preparation tools..."

cat > "${SCRIPTS_DIR}/prepare_queries.py" << 'EOF'
#!/usr/bin/env python3
"""
Prepare syn3A protein queries for BLAST analysis
"""

import os
import sys
import argparse
from pathlib import Path

def split_fasta_file(input_file, output_dir, proteins_per_file=10):
    """Split large FASTA file into smaller chunks for parallel processing"""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    proteins = []
    current_protein = None
    current_sequence = []
    
    # Parse FASTA file
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous protein
                if current_protein:
                    proteins.append({
                        'header': current_protein,
                        'sequence': ''.join(current_sequence)
                    })
                
                current_protein = line
                current_sequence = []
            else:
                current_sequence.append(line)
        
        # Save last protein
        if current_protein:
            proteins.append({
                'header': current_protein,
                'sequence': ''.join(current_sequence)
            })
    
    print(f"Parsed {len(proteins)} proteins from {input_file}")
    
    # Split into chunks
    chunk_files = []
    for i in range(0, len(proteins), proteins_per_file):
        chunk = proteins[i:i + proteins_per_file]
        chunk_num = (i // proteins_per_file) + 1
        
        chunk_file = output_dir / f"chunk_{chunk_num:03d}.fasta"
        
        with open(chunk_file, 'w') as f:
            for protein in chunk:
                f.write(f"{protein['header']}\n")
                # Write sequence in 80-character lines
                seq = protein['sequence']
                for j in range(0, len(seq), 80):
                    f.write(f"{seq[j:j+80]}\n")
        
        chunk_files.append(str(chunk_file))
        print(f"Created chunk {chunk_num}: {len(chunk)} proteins -> {chunk_file}")
    
    return chunk_files

def create_individual_queries(input_file, output_dir):
    """Create individual FASTA files for each protein"""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    protein_files = []
    current_protein = None
    current_sequence = []
    protein_count = 0
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous protein
                if current_protein:
                    protein_count += 1
                    # Extract protein ID for filename
                    protein_id = current_protein.split()[0][1:]  # Remove '>'
                    if '|' in protein_id:
                        protein_id = protein_id.split('|')[-1]
                    
                    filename = f"protein_{protein_count:03d}_{protein_id}.fasta"
                    filepath = output_dir / filename
                    
                    with open(filepath, 'w') as pf:
                        pf.write(f"{current_protein}\n")
                        seq = ''.join(current_sequence)
                        for j in range(0, len(seq), 80):
                            pf.write(f"{seq[j:j+80]}\n")
                    
                    protein_files.append(str(filepath))
                
                current_protein = line
                current_sequence = []
            else:
                current_sequence.append(line)
        
        # Save last protein
        if current_protein:
            protein_count += 1
            protein_id = current_protein.split()[0][1:]
            if '|' in protein_id:
                protein_id = protein_id.split('|')[-1]
            
            filename = f"protein_{protein_count:03d}_{protein_id}.fasta"
            filepath = output_dir / filename
            
            with open(filepath, 'w') as pf:
                pf.write(f"{current_protein}\n")
                seq = ''.join(current_sequence)
                for j in range(0, len(seq), 80):
                    pf.write(f"{seq[j:j+80]}\n")
            
            protein_files.append(str(filepath))
    
    print(f"Created {len(protein_files)} individual protein files")
    return protein_files

def main():
    parser = argparse.ArgumentParser(description='Prepare BLAST queries')
    parser.add_argument('input_fasta', help='Input FASTA file')
    parser.add_argument('--output-dir', default='queries', help='Output directory')
    parser.add_argument('--chunk-size', type=int, default=10, help='Proteins per chunk')
    parser.add_argument('--individual', action='store_true', help='Create individual protein files')
    
    args = parser.parse_args()
    
    if args.individual:
        files = create_individual_queries(args.input_fasta, args.output_dir)
    else:
        files = split_fasta_file(args.input_fasta, args.output_dir, args.chunk_size)
    
    print(f"\nQuery preparation complete!")
    print(f"Created {len(files)} query files in {args.output_dir}")

if __name__ == "__main__":
    main()
EOF

chmod +x "${SCRIPTS_DIR}/prepare_queries.py"

echo "✓ Query preparation script created"

# Test BLAST installation with a simple query
echo ""
echo "Testing BLAST installation..."

# Create test query
cat > "${PIPELINE_DIR}/test_query.fasta" << 'EOF'
>test_protein
MNVNDILKELKLSLMANKNIDESVYNDYIKTINIHKKGFSDYIVVVKSQFGLLAIKQFRQTIENEIKNIL
KEPVNISFTYEQEYKKQLEKDELINKDHSDIITKKVKKTNENTFENFVIGASNEQAFIAVQTVSKNPGIS
EOF

# Test local BLAST (this will fail if no local database, but tests the installation)
echo "Testing BLAST command..."
if blastp -help > /dev/null 2>&1; then
    echo "✓ BLAST+ command working"
else
    echo "❌ BLAST+ command failed"
    exit 1
fi

# Clean up test file
rm -f "${PIPELINE_DIR}/test_query.fasta"

echo ""
echo "=== BLAST PIPELINE SETUP COMPLETE ==="
echo ""
echo "Next steps:"
echo "1. Run: ./scripts/prepare_queries.py ../syn3A_proteins.fasta --output-dir queries"
echo "2. Run: ./scripts/run_blast_pipeline.py"
echo "3. Check results in: results/"
echo ""
echo "Pipeline directory: ${PIPELINE_DIR}"
echo "Configuration file: ${PIPELINE_DIR}/blast_config.json"
echo "Setup log: ${LOG_FILE}"
echo ""
echo "$(date): BLAST pipeline setup completed successfully"