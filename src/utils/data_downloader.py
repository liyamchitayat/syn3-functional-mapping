#!/usr/bin/env python3
"""
Utility to download various datasets for the pipeline
"""

import os
import requests
import time
from pathlib import Path

class DataDownloader:
    def __init__(self, data_dir="data"):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(exist_ok=True)
        
    def download_syn3a_proteins(self):
        """Download syn3A protein sequences from NCBI"""
        
        print("📥 Downloading syn3A proteins from NCBI...")
        
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
        output_file = self.data_dir / 'syn3A_proteins.fasta'
        with open(output_file, 'w') as f:
            f.write(response.text)
        
        print(f"✅ Downloaded {len(id_list)} proteins to {output_file}")
        return output_file

if __name__ == "__main__":
    downloader = DataDownloader()
    downloader.download_syn3a_proteins()