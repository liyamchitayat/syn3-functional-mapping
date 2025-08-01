#!/usr/bin/env python3
"""
Homology analysis script for syn3A proteins against mycoplasma and bacterial databases
"""

import requests
import time
import re
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd

def read_fasta_sequences(filename):
    """Read sequences from FASTA file"""
    sequences = []
    with open(filename, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append({
                'id': record.id,
                'description': record.description,
                'sequence': str(record.seq)
            })
    return sequences

def blast_search_web(sequence, database="nr", program="blastp", hitlist_size=50):
    """Perform BLAST search via NCBI web service"""
    try:
        print(f"Submitting BLAST search for sequence...")
        result_handle = NCBIWWW.qblast(program, database, sequence, hitlist_size=hitlist_size)
        blast_record = NCBIXML.read(result_handle)
        result_handle.close()
        return blast_record
    except Exception as e:
        print(f"Error in BLAST search: {e}")
        return None

def filter_hits_by_organism(blast_record, include_terms=None, exclude_terms=None):
    """Filter BLAST hits by organism names"""
    filtered_hits = []
    
    for alignment in blast_record.alignments:
        title = alignment.title.lower()
        
        # Check if we should include this hit
        include = True
        if include_terms:
            include = any(term.lower() in title for term in include_terms)
        
        # Check if we should exclude this hit
        if exclude_terms and include:
            include = not any(term.lower() in title for term in exclude_terms)
        
        if include and alignment.hsps:
            best_hsp = alignment.hsps[0]  # Best HSP
            filtered_hits.append({
                'title': alignment.title,
                'length': alignment.length,
                'e_value': best_hsp.expect,
                'identity': best_hsp.identities,
                'positives': best_hsp.positives,
                'query_length': blast_record.query_length,
                'identity_percent': (best_hsp.identities / blast_record.query_length) * 100,
                'coverage_percent': (len(best_hsp.query) / blast_record.query_length) * 100
            })
    
    return filtered_hits

def analyze_syn3a_homology():
    """Main analysis function"""
    print("Reading syn3A protein sequences...")
    syn3a_proteins = read_fasta_sequences('syn3A_proteins.fasta')
    print(f"Found {len(syn3a_proteins)} proteins in syn3A")
    
    results = {
        'mycoplasma_hits': [],
        'bacterial_hits': []
    }
    
    # Analyze first 5 proteins as a test
    sample_proteins = syn3a_proteins[:5]
    
    for i, protein in enumerate(sample_proteins):
        print(f"\nAnalyzing protein {i+1}/{len(sample_proteins)}: {protein['id']}")
        
        # Search against mycoplasma proteins (excluding capri)
        print("Searching against mycoplasma proteins...")
        blast_result = blast_search_web(protein['sequence'])
        
        if blast_result:
            # Filter for mycoplasma hits (excluding capri subspecies)
            mycoplasma_hits = filter_hits_by_organism(
                blast_result, 
                include_terms=['mycoplasma'], 
                exclude_terms=['capri']
            )
            
            # Filter for bacterial hits (excluding mycoplasma)
            bacterial_hits = filter_hits_by_organism(
                blast_result,
                exclude_terms=['mycoplasma']
            )
            
            results['mycoplasma_hits'].extend([
                {**hit, 'query_protein': protein['id']} 
                for hit in mycoplasma_hits[:10]  # Top 10 hits
            ])
            
            results['bacterial_hits'].extend([
                {**hit, 'query_protein': protein['id']} 
                for hit in bacterial_hits[:10]  # Top 10 hits
            ])
        
        # Rate limiting
        time.sleep(10)  # Wait 10 seconds between searches
    
    return results

def save_results(results):
    """Save results to CSV files"""
    # Save mycoplasma hits
    if results['mycoplasma_hits']:
        df_myco = pd.DataFrame(results['mycoplasma_hits'])
        df_myco.to_csv('syn3a_mycoplasma_homology.csv', index=False)
        print(f"Saved {len(results['mycoplasma_hits'])} mycoplasma hits to syn3a_mycoplasma_homology.csv")
    
    # Save bacterial hits
    if results['bacterial_hits']:
        df_bact = pd.DataFrame(results['bacterial_hits'])
        df_bact.to_csv('syn3a_bacterial_homology.csv', index=False)
        print(f"Saved {len(results['bacterial_hits'])} bacterial hits to syn3a_bacterial_homology.csv")

if __name__ == "__main__":
    try:
        results = analyze_syn3a_homology()
        save_results(results)
        print("\nHomology analysis completed!")
    except Exception as e:
        print(f"Error in analysis: {e}")