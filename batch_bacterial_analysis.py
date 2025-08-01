#!/usr/bin/env python3
"""
Batch bacterial homology analysis for all syn3A proteins
Efficient approach with batch processing and progress tracking
"""

import time
import re
import csv
import json
import os
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

def parse_all_syn3a_proteins():
    """Parse all syn3A protein sequences"""
    proteins = []
    with open('syn3A_proteins.fasta', 'r') as f:
        current_seq = ""
        current_header = ""
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq and current_header:
                    protein_id = extract_protein_id(current_header)
                    protein_name = extract_protein_name(current_header)
                    locus_tag = extract_locus_tag(current_header)
                    
                    proteins.append({
                        'protein_id': protein_id,
                        'locus_tag': locus_tag,
                        'protein_name': protein_name,
                        'sequence': current_seq,
                        'length': len(current_seq),
                        'index': len(proteins)
                    })
                
                current_header = line[1:]
                current_seq = ""
            else:
                current_seq += line
        
        # Add last protein
        if current_seq and current_header:
            protein_id = extract_protein_id(current_header)
            protein_name = extract_protein_name(current_header)
            locus_tag = extract_locus_tag(current_header)
            
            proteins.append({
                'protein_id': protein_id,
                'locus_tag': locus_tag,
                'protein_name': protein_name,
                'sequence': current_seq,
                'length': len(current_seq),
                'index': len(proteins)
            })
    
    return proteins

def extract_protein_id(header):
    match = re.search(r'protein_id=([^\]]+)', header)
    return match.group(1) if match else "Unknown"

def extract_protein_name(header):
    match = re.search(r'protein=([^\]]+)', header)
    return match.group(1) if match else "hypothetical protein"

def extract_locus_tag(header):
    match = re.search(r'locus_tag=([^\]]+)', header)
    return match.group(1) if match else "Unknown"

def batch_blast_analysis(proteins, batch_size=20, start_index=0):
    """Process proteins in batches with progress tracking"""
    
    print(f"Starting batch analysis from protein {start_index}")
    print(f"Batch size: {batch_size} proteins")
    print(f"Total proteins to process: {len(proteins) - start_index}")
    
    results = []
    
    # Load existing results if available
    progress_file = 'batch_progress.json'
    if os.path.exists(progress_file):
        with open(progress_file, 'r') as f:
            progress_data = json.load(f)
            results = progress_data.get('results', [])
            start_index = progress_data.get('last_completed_index', 0) + 1
        print(f"Resuming from protein {start_index}")
    
    for i in range(start_index, len(proteins), batch_size):
        batch = proteins[i:i + batch_size]
        batch_num = (i // batch_size) + 1
        total_batches = (len(proteins) - start_index + batch_size - 1) // batch_size
        
        print(f"\n=== BATCH {batch_num}/{total_batches} ===")
        print(f"Processing proteins {i} to {min(i + batch_size - 1, len(proteins) - 1)}")
        
        batch_results = process_batch(batch, i)
        results.extend(batch_results)
        
        # Save progress
        progress_data = {
            'last_completed_index': i + len(batch) - 1,
            'total_proteins': len(proteins),
            'results': results,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        with open(progress_file, 'w') as f:
            json.dump(progress_data, f, indent=2)
        
        # Save intermediate results
        save_batch_results(results, f"batch_results_{len(results)}_proteins.csv")
        
        print(f"Batch {batch_num} complete. Total processed: {len(results)}")
        
        # Rate limiting between batches
        if i + batch_size < len(proteins):
            print("Waiting 60 seconds before next batch...")
            time.sleep(60)
    
    return results

def process_batch(batch_proteins, start_index):
    """Process a batch of proteins"""
    batch_results = []
    
    for i, protein in enumerate(batch_proteins):
        global_index = start_index + i
        print(f"  Protein {global_index + 1}: {protein['protein_id']} ({protein['protein_name'][:50]}...)")
        
        try:
            result = analyze_protein_fast(protein, global_index)
            batch_results.append(result)
            
            # Short delay between proteins
            time.sleep(10)
            
        except Exception as e:
            print(f"    Error: {e}")
            batch_results.append({
                'protein_id': protein['protein_id'],
                'locus_tag': protein['locus_tag'],
                'protein_name': protein['protein_name'],
                'length': protein['length'],
                'index': global_index,
                'blast_status': f'error: {str(e)}',
                'bacterial_homologs_found': 0,
                'has_non_mycoplasma_homolog': False,
                'analysis_time': time.strftime('%Y-%m-%d %H:%M:%S')
            })
    
    return batch_results

def analyze_protein_fast(protein, index):
    """Fast analysis of single protein"""
    
    # Submit BLAST with timeout protection
    try:
        rid = submit_blast_fast(protein['sequence'])
        if not rid:
            return create_failed_result(protein, index, 'failed_submit')
        
        print(f"    BLAST submitted (RID: {rid})")
        
        # Wait for results with shorter timeout
        xml_result = check_blast_results_fast(rid, max_attempts=20)  # 5 minute max wait
        
        if not xml_result:
            return create_failed_result(protein, index, 'timeout')
        
        # Quick parse for bacterial hits
        bacterial_count = count_bacterial_hits(xml_result)
        has_homolog = bacterial_count > 0
        
        print(f"    Found {bacterial_count} bacterial hits")
        
        return {
            'protein_id': protein['protein_id'],
            'locus_tag': protein['locus_tag'],
            'protein_name': protein['protein_name'],
            'length': protein['length'],
            'index': index,
            'blast_status': 'success',
            'bacterial_homologs_found': bacterial_count,
            'has_non_mycoplasma_homolog': has_homolog,
            'analysis_time': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
    except Exception as e:
        return create_failed_result(protein, index, f'error: {str(e)}')

def submit_blast_fast(sequence):
    """Submit BLAST with minimal parameters for speed"""
    
    submit_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    
    params = {
        'CMD': 'Put',
        'PROGRAM': 'blastp',
        'DATABASE': 'nr',
        'QUERY': sequence,
        'FORMAT_TYPE': 'XML',
        'HITLIST_SIZE': 20,  # Fewer hits for speed
        'EXPECT': 0.001,     # Stricter cutoff
        'WORD_SIZE': 6,      # Larger word size for speed
        'ENTREZ_QUERY': 'NOT mycoplasma[Organism]'  # Exclude mycoplasma
    }
    
    try:
        data = urlencode(params).encode('utf-8')
        request = Request(submit_url, data=data)
        
        with urlopen(request, timeout=20) as response:
            result = response.read().decode('utf-8')
        
        rid_match = re.search(r'RID = (\w+)', result)
        return rid_match.group(1) if rid_match else None
        
    except Exception as e:
        print(f"    Submit error: {e}")
        return None

def check_blast_results_fast(rid, max_attempts=20):
    """Check results with shorter timeouts"""
    
    check_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID={rid}"
    
    for attempt in range(max_attempts):
        try:
            time.sleep(15)  # 15 second intervals
            
            with urlopen(check_url, timeout=20) as response:
                result = response.read().decode('utf-8')
            
            if "Status=WAITING" in result:
                if attempt % 5 == 0:  # Print every 5th attempt
                    print(f"    Waiting... ({attempt + 1}/{max_attempts})")
                continue
            elif "Status=FAILED" in result:
                return None
            elif "Status=UNKNOWN" in result:
                return None
            elif len(result) > 1000 and "<BlastOutput>" in result:
                return result
                
        except Exception as e:
            if attempt > 10:  # Only print errors after many attempts
                print(f"    Check error: {e}")
            continue
    
    return None

def count_bacterial_hits(xml_content):
    """Quickly count significant bacterial hits"""
    
    # Count hits with reasonable E-values
    hit_count = 0
    
    # Simple pattern to find hits
    hit_pattern = r'<Hit>(.*?)</Hit>'
    hits = re.findall(hit_pattern, xml_content, re.DOTALL)
    
    for hit_content in hits:
        # Check if it's not mycoplasma (double-check)
        hit_def = re.search(r'<Hit_def>(.*?)</Hit_def>', hit_content)
        if hit_def and 'mycoplasma' in hit_def.group(1).lower():
            continue
        
        # Check E-value
        evalue_match = re.search(r'<Hsp_evalue>([\d.e-]+)</Hsp_evalue>', hit_content)
        if evalue_match:
            evalue = float(evalue_match.group(1))
            if evalue < 0.001:  # Significant hit
                hit_count += 1
    
    return hit_count

def create_failed_result(protein, index, status):
    """Create result for failed analysis"""
    return {
        'protein_id': protein['protein_id'],
        'locus_tag': protein['locus_tag'],
        'protein_name': protein['protein_name'],
        'length': protein['length'],
        'index': index,
        'blast_status': status,
        'bacterial_homologs_found': 0,
        'has_non_mycoplasma_homolog': False,
        'analysis_time': time.strftime('%Y-%m-%d %H:%M:%S')
    }

def save_batch_results(results, filename):
    """Save current batch results"""
    
    with open(filename, 'w', newline='') as csvfile:
        fieldnames = [
            'index', 'protein_id', 'locus_tag', 'protein_name', 'length',
            'blast_status', 'bacterial_homologs_found', 'has_non_mycoplasma_homolog',
            'analysis_time'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

def main():
    """Main analysis function"""
    
    print("=== COMPREHENSIVE BACTERIAL HOMOLOGY ANALYSIS ===")
    print("Analyzing ALL 438 syn3A proteins for bacterial homologs")
    print("Excluding mycoplasma organisms")
    print()
    
    # Load proteins
    proteins = parse_all_syn3a_proteins()
    print(f"Loaded {len(proteins)} syn3A proteins")
    
    # Run batch analysis
    try:
        results = batch_blast_analysis(proteins, batch_size=10)  # Smaller batches
        
        # Generate final summary
        total = len(results)
        with_homologs = sum(1 for r in results if r['has_non_mycoplasma_homolog'])
        successful = sum(1 for r in results if r['blast_status'] == 'success')
        
        summary = {
            'total_proteins': total,
            'successful_analyses': successful,
            'with_bacterial_homologs': with_homologs,
            'percentage_with_homologs': (with_homologs / total * 100) if total > 0 else 0,
            'unique_to_mycoplasma': total - with_homologs,
            'analysis_date': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        with open('final_bacterial_homology_summary.json', 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Save final results
        save_batch_results(results, 'final_all_syn3a_bacterial_homology.csv')
        
        print(f"\n=== FINAL RESULTS ===")
        print(f"Total proteins analyzed: {total}")
        print(f"Successful BLAST searches: {successful}")
        print(f"Proteins with bacterial homologs: {with_homologs} ({summary['percentage_with_homologs']:.1f}%)")
        print(f"Proteins unique to mycoplasma: {total - with_homologs}")
        
    except KeyboardInterrupt:
        print("\nAnalysis interrupted. Progress saved in batch_progress.json")
    except Exception as e:
        print(f"Analysis error: {e}")

if __name__ == "__main__":
    main()