#!/usr/bin/env python3
"""
Functional analysis of syn3A proteins based on annotations
"""

import re
import csv

def parse_syn3a_proteins():
    """Parse syn3A protein annotations"""
    proteins = []
    
    with open('syn3A_proteins.fasta', 'r') as f:
        current_protein = None
        current_sequence = ""
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous protein
                if current_protein:
                    current_protein['sequence'] = current_sequence
                    current_protein['length'] = len(current_sequence)
                    proteins.append(current_protein)
                
                # Parse new protein header
                header = line[1:]
                
                # Extract key information
                protein_id_match = re.search(r'protein_id=([^]]+)', header)
                locus_tag_match = re.search(r'locus_tag=([^]]+)', header)
                protein_match = re.search(r'protein=([^]]+)', header)
                location_match = re.search(r'location=([^]]+)', header)
                
                current_protein = {
                    'header': header,
                    'protein_id': protein_id_match.group(1) if protein_id_match else '',
                    'locus_tag': locus_tag_match.group(1) if locus_tag_match else '',
                    'protein_name': protein_match.group(1) if protein_match else 'hypothetical protein',
                    'location': location_match.group(1) if location_match else '',
                }
                current_sequence = ""
            else:
                current_sequence += line
        
        # Add last protein
        if current_protein:
            current_protein['sequence'] = current_sequence
            current_protein['length'] = len(current_sequence)
            proteins.append(current_protein)
    
    return proteins

def classify_proteins(proteins):
    """Classify proteins by functional categories"""
    
    categories = {
        'DNA_replication': ['DnaA', 'DNA polymerase', 'helicase', 'primase', 'topoisomerase'],
        'RNA_processing': ['RnmV', 'KsgA', 'rRNA', 'tRNA', 'ribosom'],
        'DNA_repair': ['RecA', 'UvrA', 'UvrB', 'UvrC', 'mutS', 'mutL'],
        'Protein_synthesis': ['ribosom', 'tRNA', 'aminoacyl', 'elongation factor', 'release factor'],
        'Cell_division': ['FtsZ', 'FtsA', 'cell division', 'septum'],
        'Metabolism': ['kinase', 'synthase', 'dehydrogenase', 'transferase', 'hydrolase'],
        'Transport': ['transporter', 'permease', 'ABC', 'channel', 'pump'],
        'Transcription': ['RNA polymerase', 'sigma factor', 'transcription'],
        'Unknown': ['hypothetical', 'unknown', 'uncharacterized']
    }
    
    classified = {cat: [] for cat in categories.keys()}
    
    for protein in proteins:
        protein_name = protein['protein_name'].lower()
        assigned = False
        
        for category, keywords in categories.items():
            if any(keyword.lower() in protein_name for keyword in keywords):
                classified[category].append(protein)
                assigned = True
                break
        
        if not assigned:
            classified['Unknown'].append(protein)
    
    return classified

def analyze_essential_functions():
    """Analyze essential cellular functions in syn3A"""
    
    print("Parsing syn3A proteins...")
    proteins = parse_syn3a_proteins()
    print(f"Total proteins: {len(proteins)}")
    
    print("\nClassifying proteins by function...")
    classified = classify_proteins(proteins)
    
    # Print summary
    print("\n=== FUNCTIONAL CLASSIFICATION SUMMARY ===")
    for category, protein_list in classified.items():
        print(f"{category.replace('_', ' ')}: {len(protein_list)} proteins")
        
        # Show first few examples
        if protein_list:
            print("  Examples:")
            for i, protein in enumerate(protein_list[:3]):
                print(f"    {i+1}. {protein['protein_id']}: {protein['protein_name']}")
            if len(protein_list) > 3:
                print(f"    ... and {len(protein_list) - 3} more")
        print()
    
    # Save detailed results
    with open('syn3a_protein_classification.csv', 'w', newline='') as csvfile:
        fieldnames = ['protein_id', 'locus_tag', 'protein_name', 'length', 'category', 'location']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for category, protein_list in classified.items():
            for protein in protein_list:
                writer.writerow({
                    'protein_id': protein['protein_id'],
                    'locus_tag': protein['locus_tag'],
                    'protein_name': protein['protein_name'],
                    'length': protein['length'],
                    'category': category,
                    'location': protein['location']
                })
    
    print("Detailed results saved to syn3a_protein_classification.csv")
    
    # Identify key proteins for homology analysis
    key_proteins = []
    essential_categories = ['DNA_replication', 'RNA_processing', 'Protein_synthesis', 'Cell_division']
    
    for category in essential_categories:
        key_proteins.extend(classified[category][:2])  # Top 2 from each category
    
    print(f"\n=== KEY PROTEINS FOR HOMOLOGY ANALYSIS ===")
    for protein in key_proteins[:10]:  # Show first 10
        print(f"{protein['protein_id']}: {protein['protein_name']}")
    
    return proteins, classified

if __name__ == "__main__":
    analyze_essential_functions()