#!/usr/bin/env python3
"""
Comprehensive bacterial homology analysis for ALL 438 syn3A proteins
Goal: Determine which proteins have at least one non-mycoplasmal bacterial homolog
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
    """Parse all syn3A protein sequences from FASTA file"""
    print("Parsing ALL syn3A proteins...")
    
    proteins = []
    with open('syn3A_proteins.fasta', 'r') as f:
        current_seq = ""
        current_header = ""
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous protein
                if current_seq and current_header:
                    protein_id = extract_protein_id(current_header)
                    protein_name = extract_protein_name(current_header)
                    locus_tag = extract_locus_tag(current_header)
                    
                    proteins.append({
                        'protein_id': protein_id,
                        'locus_tag': locus_tag,
                        'protein_name': protein_name,
                        'header': current_header,
                        'sequence': current_seq,
                        'length': len(current_seq)
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
                'header': current_header,
                'sequence': current_seq,
                'length': len(current_seq)
            })
    
    print(f"Loaded {len(proteins)} syn3A proteins")
    return proteins

def extract_protein_id(header):
    """Extract protein ID from FASTA header"""
    match = re.search(r'protein_id=([^\]]+)', header)
    return match.group(1) if match else "Unknown"

def extract_protein_name(header):
    """Extract protein name from FASTA header"""
    match = re.search(r'protein=([^\]]+)', header)
    return match.group(1) if match else "hypothetical protein"

def extract_locus_tag(header):
    """Extract locus tag from FASTA header"""
    match = re.search(r'locus_tag=([^\]]+)', header)
    return match.group(1) if match else "Unknown"

def submit_blast_search(sequence, database="nr", program="blastp"):
    """Submit BLAST search to NCBI"""
    
    submit_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    
    # BLAST parameters
    params = {
        'CMD': 'Put',
        'PROGRAM': program,
        'DATABASE': database,
        'QUERY': sequence,
        'FORMAT_TYPE': 'XML',
        'HITLIST_SIZE': 100,  # More hits for comprehensive analysis  
        'EXPECT': 0.01,       # Stricter E-value
        'WORD_SIZE': 3,
        'COMPOSITION_BASED_STATISTICS': 'on',
        'FILTER': 'L',        # Low complexity filter
        'ENTREZ_QUERY': 'NOT mycoplasma[Organism]'  # Exclude mycoplasma
    }
    
    try:
        data = urlencode(params).encode('utf-8')
        request = Request(submit_url, data=data)
        
        with urlopen(request, timeout=30) as response:
            result = response.read().decode('utf-8')
        
        # Extract RID (Request ID)
        rid_match = re.search(r'RID = (\w+)', result)
        if rid_match:
            return rid_match.group(1)
        else:
            return None
            
    except (HTTPError, URLError, Exception) as e:
        print(f"Error submitting BLAST: {e}")
        return None

def check_blast_results(rid, max_attempts=60):
    """Check BLAST results and retrieve when ready"""
    
    check_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    params = {
        'CMD': 'Get',
        'FORMAT_TYPE': 'XML',
        'RID': rid
    }
    
    url = f"{check_url}?{urlencode(params)}"
    
    for attempt in range(max_attempts):
        try:
            time.sleep(15)  # Wait 15 seconds between checks
            
            with urlopen(url, timeout=30) as response:
                result = response.read().decode('utf-8')
            
            if "Status=WAITING" in result:
                print(f"  Still waiting... (attempt {attempt + 1}/{max_attempts})")
                continue
            elif "Status=FAILED" in result:
                print("  BLAST search failed")
                return None
            elif "Status=UNKNOWN" in result:
                print("  BLAST RID expired")
                return None
            elif len(result) > 1000 and "<BlastOutput>" in result:
                print("  BLAST results ready!")
                return result
            else:
                continue
                
        except Exception as e:
            print(f"  Error checking results: {e}")
            continue
    
    print("  Timeout waiting for results")
    return None

def parse_blast_xml_simple(xml_content):
    """Simple XML parsing to extract bacterial hits"""
    
    bacterial_hits = []
    
    # Extract all hit entries
    hit_pattern = r'<Hit>(.*?)</Hit>'
    hits = re.findall(hit_pattern, xml_content, re.DOTALL)
    
    for hit_content in hits:
        # Extract hit information
        hit_def = re.search(r'<Hit_def>(.*?)</Hit_def>', hit_content)
        hit_len = re.search(r'<Hit_len>(\d+)</Hit_len>', hit_content)
        
        # Extract best HSP (High-scoring Segment Pair)
        hsp_pattern = r'<Hsp>(.*?)</Hsp>'
        hsps = re.findall(hsp_pattern, hit_content, re.DOTALL)
        
        if hit_def and hsps:
            hit_description = hit_def.group(1)
            hit_length = int(hit_len.group(1)) if hit_len else 0
            
            # Get best HSP
            best_hsp = hsps[0]
            
            # Extract HSP details
            evalue = re.search(r'<Hsp_evalue>([\d.e-]+)</Hsp_evalue>', best_hsp)
            identity = re.search(r'<Hsp_identity>(\d+)</Hsp_identity>', best_hsp)
            align_len = re.search(r'<Hsp_align-len>(\d+)</Hsp_align-len>', best_hsp)
            query_len = re.search(r'<Hsp_query-len>(\d+)</Hsp_query-len>', best_hsp)
            
            if evalue and identity and align_len:
                # Skip mycoplasma hits (double-check)
                if 'mycoplasma' in hit_description.lower():
                    continue
                
                bacterial_hits.append({
                    'description': hit_description,
                    'length': hit_length,
                    'evalue': float(evalue.group(1)),
                    'identity': int(identity.group(1)),
                    'align_length': int(align_len.group(1)),
                    'query_length': int(query_len.group(1)) if query_len else 0,
                    'identity_percent': (int(identity.group(1)) / int(align_len.group(1))) * 100 if align_len else 0
                })
    
    return bacterial_hits

def analyze_single_protein(protein, protein_index, total_proteins):
    """Analyze a single protein for bacterial homologs"""
    
    print(f"\n=== Protein {protein_index + 1}/{total_proteins} ===")
    print(f"ID: {protein['protein_id']}")
    print(f"Name: {protein['protein_name']}")
    print(f"Length: {protein['length']} amino acids")
    
    # Submit BLAST search
    print("Submitting BLAST search...")
    rid = submit_blast_search(protein['sequence'])
    
    if not rid:
        print("Failed to submit BLAST search")
        return {
            'protein_id': protein['protein_id'],
            'protein_name': protein['protein_name'],
            'length': protein['length'],
            'blast_status': 'failed_submit',
            'bacterial_homologs_found': 0,
            'has_non_mycoplasma_homolog': False,
            'best_bacterial_hits': []
        }
    
    print(f"BLAST submitted with RID: {rid}")
    
    # Wait for results
    xml_result = check_blast_results(rid)
    
    if not xml_result:
        print("Failed to retrieve BLAST results")
        return {
            'protein_id': protein['protein_id'],
            'protein_name': protein['protein_name'],
            'length': protein['length'],
            'blast_status': 'failed_retrieve',
            'bacterial_homologs_found': 0,
            'has_non_mycoplasma_homolog': False,
            'best_bacterial_hits': []
        }
    
    # Parse results
    bacterial_hits = parse_blast_xml_simple(xml_result)
    
    # Filter significant hits (E-value < 0.01, identity > 30%)
    significant_hits = [
        hit for hit in bacterial_hits 
        if hit['evalue'] < 0.01 and hit['identity_percent'] > 30.0
    ]
    
    # Sort by E-value (best first)
    significant_hits.sort(key=lambda x: x['evalue'])
    
    print(f"Found {len(significant_hits)} significant bacterial hits")
    
    # Save raw XML for this protein
    xml_filename = f"blast_results/protein_{protein['protein_id']}_blast.xml"
    os.makedirs('blast_results', exist_ok=True)
    with open(xml_filename, 'w') as f:
        f.write(xml_result)
    
    return {
        'protein_id': protein['protein_id'],
        'locus_tag': protein['locus_tag'],
        'protein_name': protein['protein_name'],
        'length': protein['length'],
        'blast_status': 'success',
        'bacterial_homologs_found': len(significant_hits),
        'has_non_mycoplasma_homolog': len(significant_hits) > 0,
        'best_bacterial_hits': significant_hits[:10],  # Top 10 hits
        'blast_xml_file': xml_filename
    }

def comprehensive_bacterial_analysis():
    """Main function to analyze all syn3A proteins"""
    
    print("=== COMPREHENSIVE BACTERIAL HOMOLOGY ANALYSIS ===")
    print("Analyzing ALL 438 syn3A proteins against ALL bacterial proteins")
    print("Excluding mycoplasma organisms")
    print()
    
    # Load all proteins
    all_proteins = parse_all_syn3a_proteins()
    
    # Initialize results
    all_results = []
    
    # Process each protein
    start_time = time.time()
    
    for i, protein in enumerate(all_proteins):
        try:
            result = analyze_single_protein(protein, i, len(all_proteins))
            all_results.append(result)
            
            # Save progress periodically
            if (i + 1) % 10 == 0:
                save_intermediate_results(all_results, i + 1)
                
                elapsed = time.time() - start_time
                avg_time = elapsed / (i + 1)
                remaining = avg_time * (len(all_proteins) - i - 1)
                
                print(f"\n--- Progress Update ---")
                print(f"Completed: {i + 1}/{len(all_proteins)} proteins")
                print(f"Time elapsed: {elapsed/60:.1f} minutes")
                print(f"Estimated remaining: {remaining/60:.1f} minutes")
                print(f"With homologs: {sum(1 for r in all_results if r['has_non_mycoplasma_homolog'])}")
                print("--- --- --- --- --- ---\n")
            
            # Rate limiting - wait between searches
            time.sleep(5)
            
        except KeyboardInterrupt:
            print("\nAnalysis interrupted by user")
            print(f"Completed {len(all_results)} proteins so far")
            break
        except Exception as e:
            print(f"Error analyzing protein {protein['protein_id']}: {e}")
            # Add failed result
            all_results.append({
                'protein_id': protein['protein_id'],
                'protein_name': protein['protein_name'],
                'length': protein['length'],
                'blast_status': f'error: {str(e)}',
                'bacterial_homologs_found': 0,
                'has_non_mycoplasma_homolog': False,
                'best_bacterial_hits': []
            })
    
    return all_results

def save_intermediate_results(results, num_completed):
    """Save intermediate results"""
    filename = f'intermediate_results_{num_completed}_proteins.csv'
    
    with open(filename, 'w', newline='') as csvfile:
        fieldnames = [
            'protein_id', 'locus_tag', 'protein_name', 'length', 
            'blast_status', 'bacterial_homologs_found', 'has_non_mycoplasma_homolog',
            'best_hit_description', 'best_hit_evalue', 'best_hit_identity_percent'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for result in results:
            best_hit = result['best_bacterial_hits'][0] if result['best_bacterial_hits'] else {}
            
            writer.writerow({
                'protein_id': result['protein_id'],
                'locus_tag': result.get('locus_tag', ''),
                'protein_name': result['protein_name'],
                'length': result['length'],
                'blast_status': result['blast_status'],
                'bacterial_homologs_found': result['bacterial_homologs_found'],
                'has_non_mycoplasma_homolog': result['has_non_mycoplasma_homolog'],
                'best_hit_description': best_hit.get('description', ''),
                'best_hit_evalue': best_hit.get('evalue', ''),
                'best_hit_identity_percent': best_hit.get('identity_percent', '')
            })

def save_final_results(results):
    """Save comprehensive final results"""
    
    print("\n=== SAVING FINAL RESULTS ===")
    
    # Summary statistics
    total_proteins = len(results)
    with_homologs = sum(1 for r in results if r['has_non_mycoplasma_homolog'])
    successful_searches = sum(1 for r in results if r['blast_status'] == 'success')
    
    summary = {
        'analysis_date': time.strftime('%Y-%m-%d %H:%M:%S'),
        'total_syn3a_proteins': total_proteins,
        'successful_blast_searches': successful_searches,
        'proteins_with_bacterial_homologs': with_homologs,
        'percentage_with_homologs': (with_homologs / total_proteins) * 100,
        'proteins_without_bacterial_homologs': total_proteins - with_homologs,
        'percentage_without_homologs': ((total_proteins - with_homologs) / total_proteins) * 100
    }
    
    # Save summary
    with open('comprehensive_bacterial_homology_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Save detailed results
    with open('all_syn3a_bacterial_homology.csv', 'w', newline='') as csvfile:
        fieldnames = [
            'protein_id', 'locus_tag', 'protein_name', 'length', 
            'blast_status', 'bacterial_homologs_found', 'has_non_mycoplasma_homolog',
            'best_hit_description', 'best_hit_evalue', 'best_hit_identity_percent',
            'second_hit_description', 'second_hit_evalue', 'second_hit_identity_percent'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for result in results:
            best_hit = result['best_bacterial_hits'][0] if len(result['best_bacterial_hits']) > 0 else {}
            second_hit = result['best_bacterial_hits'][1] if len(result['best_bacterial_hits']) > 1 else {}
            
            writer.writerow({
                'protein_id': result['protein_id'],
                'locus_tag': result.get('locus_tag', ''),
                'protein_name': result['protein_name'],
                'length': result['length'],
                'blast_status': result['blast_status'],
                'bacterial_homologs_found': result['bacterial_homologs_found'],
                'has_non_mycoplasma_homolog': result['has_non_mycoplasma_homolog'],
                'best_hit_description': best_hit.get('description', ''),
                'best_hit_evalue': best_hit.get('evalue', ''),
                'best_hit_identity_percent': best_hit.get('identity_percent', ''),
                'second_hit_description': second_hit.get('description', ''),
                'second_hit_evalue': second_hit.get('evalue', ''),
                'second_hit_identity_percent': second_hit.get('identity_percent', '')
            })
    
    # Save proteins WITHOUT bacterial homologs (interesting for minimal genome analysis)
    without_homologs = [r for r in results if not r['has_non_mycoplasma_homolog']]
    
    with open('syn3a_proteins_no_bacterial_homologs.csv', 'w', newline='') as csvfile:
        fieldnames = ['protein_id', 'locus_tag', 'protein_name', 'length', 'blast_status']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for result in without_homologs:
            writer.writerow({
                'protein_id': result['protein_id'],
                'locus_tag': result.get('locus_tag', ''),
                'protein_name': result['protein_name'],
                'length': result['length'],
                'blast_status': result['blast_status']
            })
    
    print(f"Analysis complete!")
    print(f"Total proteins analyzed: {total_proteins}")
    print(f"Successful BLAST searches: {successful_searches}")
    print(f"Proteins with bacterial homologs: {with_homologs} ({summary['percentage_with_homologs']:.1f}%)")
    print(f"Proteins without bacterial homologs: {total_proteins - with_homologs} ({summary['percentage_without_homologs']:.1f}%)")
    
    print(f"\nFiles created:")
    print(f"- all_syn3a_bacterial_homology.csv (detailed results)")
    print(f"- syn3a_proteins_no_bacterial_homologs.csv (proteins unique to mycoplasma)")
    print(f"- comprehensive_bacterial_homology_summary.json (summary statistics)")
    print(f"- blast_results/ directory (raw BLAST XML files)")

if __name__ == "__main__":
    try:
        results = comprehensive_bacterial_analysis()
        save_final_results(results)
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user")
    except Exception as e:
        print(f"Analysis failed: {e}")