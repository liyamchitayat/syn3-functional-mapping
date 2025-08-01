#!/usr/bin/env python3
"""
Download syn3A proteins from NCBI
"""

import requests
import time

def download_syn3a_proteins():
    """Download syn3A protein sequences from NCBI"""
    
    print("Downloading syn3A proteins from NCBI...")
    
    # NCBI E-utilities
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    # Search for syn3A proteins
    search_url = f"{base_url}esearch.fcgi"
    search_params = {
        'db': 'protein',
        'term': 'Mesoplasma florum JCVI-syn3A[Organism]',
        'retmode': 'json',
        'retmax': 500
    }
    
    response = requests.get(search_url, params=search_params)
    data = response.json()
    
    id_list = data['esearchresult']['idlist']
    print(f"Found {len(id_list)} proteins")
    
    # Fetch sequences
    fetch_url = f"{base_url}efetch.fcgi"
    fetch_params = {
        'db': 'protein',
        'id': ','.join(id_list),
        'rettype': 'fasta',
        'retmode': 'text'
    }
    
    response = requests.get(fetch_url, params=fetch_params)
    
    # Save to file
    with open('syn3A_proteins.fasta', 'w') as f:
        f.write(response.text)
    
    print(f"✅ Downloaded {len(id_list)} proteins to syn3A_proteins.fasta")

if __name__ == "__main__":
    download_syn3a_proteins()