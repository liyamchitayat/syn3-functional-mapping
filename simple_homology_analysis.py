#!/usr/bin/env python3
"""
Simple homology analysis using NCBI BLAST web API
"""

import requests
import time
import re
import json

def parse_fasta(filename):
    """Simple FASTA parser"""
    sequences = []
    with open(filename, 'r') as f:
        current_seq = ""
        current_id = ""
        current_desc = ""
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq and current_id:
                    sequences.append({
                        'id': current_id,
                        'description': current_desc,
                        'sequence': current_seq
                    })
                
                parts = line[1:].split(' ', 1)
                current_id = parts[0]
                current_desc = parts[1] if len(parts) > 1 else ""
                current_seq = ""
            else:
                current_seq += line
        
        if current_seq and current_id:
            sequences.append({
                'id': current_id,
                'description': current_desc,
                'sequence': current_seq
            })
    
    return sequences

def blast_search_simple(sequence, database="nr", program="blastp"):
    """Simple BLAST search using web API"""
    
    # Submit BLAST job
    submit_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    submit_params = {
        'CMD': 'Put',
        'PROGRAM': program,
        'DATABASE': database,
        'QUERY': sequence,
        'FORMAT_TYPE': 'XML',
        'HITLIST_SIZE': 50,
        'EXPECT': 0.01,
        'WORD_SIZE': 3,
        'COMPOSITION_BASED_STATISTICS': 'on'
    }
    
    print("Submitting BLAST search...")
    response = requests.post(submit_url, data=submit_params)
    
    if response.status_code != 200:
        print(f"Error submitting BLAST: {response.status_code}")
        return None
    
    # Extract RID (Request ID)
    rid_match = re.search(r'RID = (\w+)', response.text)
    if not rid_match:
        print("Could not extract RID from response")
        return None
    
    rid = rid_match.group(1)
    print(f"BLAST job submitted with RID: {rid}")
    
    # Poll for results
    check_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    check_params = {
        'CMD': 'Get',
        'FORMAT_TYPE': 'XML',
        'RID': rid
    }
    
    max_attempts = 60  # Maximum 10 minutes waiting
    attempt = 0
    
    while attempt < max_attempts:
        print(f"Checking results... (attempt {attempt + 1})")
        time.sleep(10)  # Wait 10 seconds between checks
        
        response = requests.get(check_url, params=check_params)
        
        if "Status=WAITING" in response.text:
            attempt += 1
            continue
        elif "Status=FAILED" in response.text:
            print("BLAST search failed")
            return None
        elif "Status=UNKNOWN" in response.text:
            print("BLAST RID expired or unknown")
            return None
        elif response.status_code == 200 and len(response.text) > 1000:
            print("BLAST results received!")
            return response.text
        else:
            attempt += 1
    
    print("Timeout waiting for BLAST results")
    return None

def analyze_first_protein():
    """Analyze the first syn3A protein as a test"""
    print("Reading syn3A protein sequences...")
    proteins = parse_fasta('syn3A_proteins.fasta')
    print(f"Found {len(proteins)} proteins")
    
    if not proteins:
        print("No proteins found!")
        return
    
    # Analyze first protein
    first_protein = proteins[0]
    print(f"\nAnalyzing first protein: {first_protein['id']}")
    print(f"Description: {first_protein['description']}")
    print(f"Sequence length: {len(first_protein['sequence'])} amino acids")
    
    # Perform BLAST search
    blast_result = blast_search_simple(first_protein['sequence'])
    
    if blast_result:
        # Save raw XML result
        with open('blast_result_protein1.xml', 'w') as f:
            f.write(blast_result)
        print("BLAST results saved to blast_result_protein1.xml")
        
        # Basic parsing of results
        mycoplasma_count = blast_result.lower().count('mycoplasma')
        print(f"Found {mycoplasma_count} hits containing 'mycoplasma'")
        
        return blast_result
    else:
        print("BLAST search failed")
        return None

if __name__ == "__main__":
    analyze_first_protein()