#!/usr/bin/env python3
"""
Test bacterial homology analysis on first 20 syn3A proteins
"""

import time
import re
import csv
import json
from urllib.parse import urlencode
from urllib.request import urlopen, Request

def parse_first_n_proteins(n=20):
    """Parse first N syn3A proteins for testing"""
    proteins = []
    with open('syn3A_proteins.fasta', 'r') as f:
        current_seq = ""
        current_header = ""
        protein_count = 0
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq and current_header and protein_count < n:
                    protein_id = extract_protein_id(current_header)
                    protein_name = extract_protein_name(current_header)
                    
                    proteins.append({
                        'protein_id': protein_id,
                        'protein_name': protein_name,
                        'sequence': current_seq,
                        'length': len(current_seq)
                    })
                    protein_count += 1
                
                if protein_count >= n:
                    break
                    
                current_header = line[1:]
                current_seq = ""
            else:
                current_seq += line
        
        # Add last protein if within limit
        if current_seq and current_header and protein_count < n:
            protein_id = extract_protein_id(current_header)
            protein_name = extract_protein_name(current_header)
            
            proteins.append({
                'protein_id': protein_id,
                'protein_name': protein_name,
                'sequence': current_seq,
                'length': len(current_seq)
            })
    
    return proteins

def extract_protein_id(header):
    match = re.search(r'protein_id=([^\]]+)', header)
    return match.group(1) if match else "Unknown"

def extract_protein_name(header):
    match = re.search(r'protein=([^\]]+)', header)
    return match.group(1) if match else "hypothetical protein"

def quick_blast_test(sequence):
    """Quick BLAST test to see if we get bacterial hits"""
    
    submit_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    
    params = {
        'CMD': 'Put',
        'PROGRAM': 'blastp',
        'DATABASE': 'nr',
        'QUERY': sequence,
        'FORMAT_TYPE': 'XML',
        'HITLIST_SIZE': 10,
        'EXPECT': 0.01,
        'ENTREZ_QUERY': 'NOT mycoplasma[Organism]'
    }
    
    try:
        data = urlencode(params).encode('utf-8')
        request = Request(submit_url, data=data)
        
        with urlopen(request, timeout=30) as response:
            result = response.read().decode('utf-8')
        
        rid_match = re.search(r'RID = (\w+)', result)
        if not rid_match:
            return None, "No RID found"
        
        rid = rid_match.group(1)
        
        # Wait for results
        check_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID={rid}"
        
        for attempt in range(30):  # 7.5 minute max wait
            time.sleep(15)
            
            with urlopen(check_url, timeout=20) as response:
                xml_result = response.read().decode('utf-8')
            
            if "Status=WAITING" in xml_result:
                continue
            elif "Status=FAILED" in xml_result:
                return None, "BLAST failed"
            elif len(xml_result) > 1000 and "<BlastOutput>" in xml_result:
                # Count hits quickly
                hit_count = xml_result.count('<Hit>')
                return hit_count, "Success"
        
        return None, "Timeout"
        
    except Exception as e:
        return None, f"Error: {str(e)}"

def test_analysis():
    """Test analysis on first few proteins"""
    
    print("=== TESTING BACTERIAL HOMOLOGY ANALYSIS ===")
    print("Testing on first 5 syn3A proteins")
    print()
    
    # Load test proteins
    test_proteins = parse_first_n_proteins(5)
    print(f"Loaded {len(test_proteins)} test proteins")
    
    results = []
    
    for i, protein in enumerate(test_proteins):
        print(f"\n--- Testing protein {i+1}/5 ---")
        print(f"ID: {protein['protein_id']}")
        print(f"Name: {protein['protein_name'][:60]}...")
        print(f"Length: {protein['length']} AA")
        
        start_time = time.time()
        hit_count, status = quick_blast_test(protein['sequence'])
        elapsed = time.time() - start_time
        
        result = {
            'protein_id': protein['protein_id'],
            'protein_name': protein['protein_name'],
            'length': protein['length'],
            'bacterial_hits': hit_count if hit_count else 0,
            'has_bacterial_homolog': hit_count > 0 if hit_count else False,
            'status': status,
            'analysis_time_seconds': round(elapsed, 1)
        }
        
        results.append(result)
        
        print(f"Result: {hit_count} bacterial hits ({status})")
        print(f"Time: {elapsed:.1f} seconds")
        
        # Rate limiting
        if i < len(test_proteins) - 1:
            print("Waiting 30 seconds...")
            time.sleep(30)
    
    # Save test results
    with open('test_bacterial_homology_results.csv', 'w', newline='') as csvfile:
        fieldnames = ['protein_id', 'protein_name', 'length', 'bacterial_hits', 
                     'has_bacterial_homolog', 'status', 'analysis_time_seconds']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    # Summary
    successful = [r for r in results if r['status'] == 'Success']
    with_homologs = [r for r in successful if r['has_bacterial_homolog']]
    
    print(f"\n=== TEST RESULTS SUMMARY ===")
    print(f"Total proteins tested: {len(results)}")
    print(f"Successful analyses: {len(successful)}")
    print(f"With bacterial homologs: {len(with_homologs)}")
    
    if successful:
        success_rate = len(with_homologs) / len(successful) * 100
        print(f"Success rate for homolog detection: {success_rate:.1f}%")
        
        avg_time = sum(r['analysis_time_seconds'] for r in successful) / len(successful)
        print(f"Average analysis time: {avg_time:.1f} seconds")
        
        total_time_estimate = avg_time * 438 / 3600  # Convert to hours
        print(f"Estimated time for all 438 proteins: {total_time_estimate:.1f} hours")
    
    print(f"\nResults saved to: test_bacterial_homology_results.csv")
    
    return results

if __name__ == "__main__":
    test_analysis()