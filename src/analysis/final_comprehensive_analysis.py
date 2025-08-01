#!/usr/bin/env python3
"""
Final comprehensive bacterial homology analysis for all 438 syn3A proteins
Using functional annotation patterns and known protein families
"""

import csv
import json
import re
import time

def parse_all_syn3a_proteins():
    """Parse all syn3A proteins with detailed annotations"""
    proteins = []
    
    with open('syn3A_proteins.fasta', 'r') as f:
        current_seq = ""
        current_header = ""
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq and current_header:
                    protein_info = parse_protein_header(current_header)
                    protein_info['sequence'] = current_seq
                    protein_info['length'] = len(current_seq)
                    proteins.append(protein_info)
                
                current_header = line[1:]
                current_seq = ""
            else:
                current_seq += line
        
        # Add last protein
        if current_seq and current_header:
            protein_info = parse_protein_header(current_header)
            protein_info['sequence'] = current_seq
            protein_info['length'] = len(current_seq)
            proteins.append(protein_info)
    
    return proteins

def parse_protein_header(header):
    """Parse protein header to extract all information"""
    
    protein_id_match = re.search(r'protein_id=([^\]]+)', header)
    locus_tag_match = re.search(r'locus_tag=([^\]]+)', header)
    protein_match = re.search(r'protein=([^\]]+)', header)
    location_match = re.search(r'location=([^\]]+)', header)
    
    return {
        'protein_id': protein_id_match.group(1) if protein_id_match else 'Unknown',
        'locus_tag': locus_tag_match.group(1) if locus_tag_match else 'Unknown',
        'protein_name': protein_match.group(1) if protein_match else 'hypothetical protein',
        'location': location_match.group(1) if location_match else 'Unknown',
        'full_header': header
    }

def classify_bacterial_homology_likelihood(protein):
    """Classify likelihood of bacterial homology based on protein function"""
    
    protein_name = protein['protein_name'].lower()
    
    # Universal/Essential proteins - very likely to have bacterial homologs
    universal_keywords = [
        'dnaa', 'dna polymerase', 'gyrase', 'rna polymerase', 'ribosom', 'trna', 'rrna',
        'elongation factor', 'release factor', 'aminoacyl', 'peptidyl', 'initiation factor',
        'chaperone', 'hsp', 'groel', 'ftsa', 'ftsz', 'reca', 'uvr', 'muts', 'mutl',
        'pyruvate', 'glycolysis', 'tca cycle', 'atp synthase', 'nadh', 'thioredoxin'
    ]
    
    # Metabolic enzymes - likely to have bacterial homologs
    metabolic_keywords = [
        'kinase', 'synthase', 'synthetase', 'dehydrogenase', 'reductase', 'oxidase',
        'transferase', 'hydrolase', 'lyase', 'isomerase', 'ligase', 'phosphatase',
        'aldolase', 'enolase', 'mutase', 'epimerase', 'carboxylase'
    ]
    
    # Transport proteins - moderately likely
    transport_keywords = [
        'transporter', 'permease', 'abc', 'channel', 'pump', 'porter', 'symporter',
        'antiporter', 'exchanger'
    ]
    
    # DNA/RNA processing - very likely
    nucleic_acid_keywords = [
        'helicase', 'primase', 'ligase', 'topoisomerase', 'exonuclease', 'endonuclease',
        'methyltransferase', 'demethylase', 'pseudouridine', 'modification'
    ]
    
    # Transcription factors - moderately likely (can be specific)
    transcription_keywords = [
        'transcription', 'regulator', 'repressor', 'activator', 'sigma', 'rho'
    ]
    
    # Cell wall/membrane - variable (mycoplasma lack cell walls)
    cell_structure_keywords = [
        'cell wall', 'peptidoglycan', 'murein', 'lipopolysaccharide'
    ]
    
    # Mycoplasma-specific or unique
    mycoplasma_specific = [
        'lipoprotein', 'vlp', 'variable', 'surface', 'adhesin', 'cytadherence'
    ]
    
    # Check categories
    homology_score = 0
    homology_category = "Unknown"
    reasoning = []
    
    if any(keyword in protein_name for keyword in universal_keywords):
        homology_score = 95
        homology_category = "Universal/Essential"
        reasoning.append("Universal cellular function")
    
    elif any(keyword in protein_name for keyword in nucleic_acid_keywords):
        homology_score = 90
        homology_category = "DNA/RNA Processing"
        reasoning.append("Essential DNA/RNA processing")
    
    elif any(keyword in protein_name for keyword in metabolic_keywords):
        homology_score = 85
        homology_category = "Metabolic Enzyme"
        reasoning.append("Metabolic enzyme")
    
    elif any(keyword in protein_name for keyword in transport_keywords):
        homology_score = 70
        homology_category = "Transport"
        reasoning.append("Transport function")
    
    elif any(keyword in protein_name for keyword in transcription_keywords):
        homology_score = 65
        homology_category = "Transcription"
        reasoning.append("Transcriptional control")
    
    elif any(keyword in protein_name for keyword in mycoplasma_specific):
        homology_score = 20
        homology_category = "Mycoplasma-specific"
        reasoning.append("Mycoplasma-specific function")
    
    elif any(keyword in protein_name for keyword in cell_structure_keywords):
        homology_score = 30
        homology_category = "Cell Structure"
        reasoning.append("Cell structure (mycoplasma lack cell walls)")
    
    elif 'hypothetical' in protein_name or 'unknown' in protein_name:
        homology_score = 40
        homology_category = "Hypothetical"
        reasoning.append("Unknown function")
    
    else:
        homology_score = 50
        homology_category = "Other"
        reasoning.append("Unclassified function")
    
    # Additional scoring based on protein characteristics
    if protein['length'] > 500:
        homology_score += 5  # Larger proteins more likely conserved
        reasoning.append("Large protein")
    elif protein['length'] < 100:
        homology_score -= 10  # Very small proteins may be unique
        reasoning.append("Small protein")
    
    # Cap score at 100
    homology_score = min(homology_score, 100)
    
    return {
        'homology_likelihood_score': homology_score,
        'homology_category': homology_category,
        'predicted_has_bacterial_homolog': homology_score >= 60,
        'confidence_level': 'High' if homology_score >= 80 else 'Medium' if homology_score >= 60 else 'Low',
        'reasoning': '; '.join(reasoning)
    }

def comprehensive_analysis():
    """Perform comprehensive analysis of all syn3A proteins"""
    
    print("=== COMPREHENSIVE BACTERIAL HOMOLOGY ANALYSIS ===")
    print("Analyzing ALL 438 syn3A proteins for bacterial homology likelihood")
    print("Using functional annotation patterns and protein family analysis")
    print()
    
    # Load all proteins
    proteins = parse_all_syn3a_proteins()
    print(f"Loaded {len(proteins)} syn3A proteins")
    
    # Load existing functional classification
    classifications = {}
    try:
        with open('syn3a_protein_classification.csv', 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                classifications[row['protein_id']] = row['category']
    except FileNotFoundError:
        print("Using default classifications")
    
    # Analyze each protein
    results = []
    
    for i, protein in enumerate(proteins):
        print(f"Analyzing protein {i+1}/{len(proteins)}: {protein['protein_id']}")
        
        # Get functional category
        functional_category = classifications.get(protein['protein_id'], 'Unknown')
        
        # Predict bacterial homology
        homology_analysis = classify_bacterial_homology_likelihood(protein)
        
        result = {
            'index': i,
            'protein_id': protein['protein_id'],
            'locus_tag': protein['locus_tag'],
            'protein_name': protein['protein_name'],
            'length': protein['length'],
            'functional_category': functional_category,
            **homology_analysis,
            'analysis_method': 'functional_annotation_pattern',
            'analysis_date': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        results.append(result)
    
    return results

def save_comprehensive_results(results):
    """Save comprehensive analysis results"""
    
    print(f"\n=== SAVING COMPREHENSIVE RESULTS ===")
    
    # Calculate summary statistics
    total_proteins = len(results)
    predicted_with_homologs = sum(1 for r in results if r['predicted_has_bacterial_homolog'])
    high_confidence = sum(1 for r in results if r['confidence_level'] == 'High')
    medium_confidence = sum(1 for r in results if r['confidence_level'] == 'Medium')
    low_confidence = sum(1 for r in results if r['confidence_level'] == 'Low')
    
    # Category breakdown
    category_stats = {}
    for result in results:
        category = result['homology_category']
        if category not in category_stats:
            category_stats[category] = {'count': 0, 'with_homologs': 0}
        category_stats[category]['count'] += 1
        if result['predicted_has_bacterial_homolog']:
            category_stats[category]['with_homologs'] += 1
    
    # Summary statistics
    summary = {
        'analysis_date': time.strftime('%Y-%m-%d %H:%M:%S'),
        'analysis_method': 'functional_annotation_based_prediction',
        'total_syn3a_proteins': total_proteins,
        'predicted_with_bacterial_homologs': predicted_with_homologs,
        'percentage_with_homologs': (predicted_with_homologs / total_proteins) * 100,
        'predicted_mycoplasma_specific': total_proteins - predicted_with_homologs,
        'percentage_mycoplasma_specific': ((total_proteins - predicted_with_homologs) / total_proteins) * 100,
        'confidence_levels': {
            'high_confidence_predictions': high_confidence,
            'medium_confidence_predictions': medium_confidence,
            'low_confidence_predictions': low_confidence
        },
        'homology_categories': category_stats
    }
    
    # Save summary
    with open('comprehensive_bacterial_homology_analysis_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Save detailed results
    with open('all_syn3a_bacterial_homology_predictions.csv', 'w', newline='') as csvfile:
        fieldnames = [
            'index', 'protein_id', 'locus_tag', 'protein_name', 'length',
            'functional_category', 'homology_likelihood_score', 'homology_category',
            'predicted_has_bacterial_homolog', 'confidence_level', 'reasoning',
            'analysis_method', 'analysis_date'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    # Save proteins predicted to be mycoplasma-specific
    mycoplasma_specific = [r for r in results if not r['predicted_has_bacterial_homolog']]
    
    with open('predicted_mycoplasma_specific_proteins.csv', 'w', newline='') as csvfile:
        fieldnames = ['protein_id', 'locus_tag', 'protein_name', 'length', 
                     'homology_category', 'homology_likelihood_score', 'reasoning']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for result in mycoplasma_specific:
            writer.writerow({
                'protein_id': result['protein_id'],
                'locus_tag': result['locus_tag'],
                'protein_name': result['protein_name'],
                'length': result['length'],
                'homology_category': result['homology_category'],
                'homology_likelihood_score': result['homology_likelihood_score'],
                'reasoning': result['reasoning']
            })
    
    # Print summary
    print(f"Analysis completed for {total_proteins} proteins")
    print(f"Predicted to have bacterial homologs: {predicted_with_homologs} ({summary['percentage_with_homologs']:.1f}%)")
    print(f"Predicted mycoplasma-specific: {len(mycoplasma_specific)} ({summary['percentage_mycoplasma_specific']:.1f}%)")
    print(f"\nConfidence levels:")
    print(f"  High confidence: {high_confidence} proteins")
    print(f"  Medium confidence: {medium_confidence} proteins") 
    print(f"  Low confidence: {low_confidence} proteins")
    
    print(f"\nTop homology categories:")
    for category, stats in sorted(category_stats.items(), key=lambda x: x[1]['count'], reverse=True)[:5]:
        percentage = (stats['with_homologs'] / stats['count']) * 100
        print(f"  {category}: {stats['count']} proteins ({stats['with_homologs']} with predicted homologs, {percentage:.1f}%)")
    
    print(f"\nFiles created:")
    print(f"  - all_syn3a_bacterial_homology_predictions.csv (detailed results)")
    print(f"  - predicted_mycoplasma_specific_proteins.csv ({len(mycoplasma_specific)} proteins)")
    print(f"  - comprehensive_bacterial_homology_analysis_summary.json (summary)")

def main():
    """Main analysis function"""
    
    try:
        results = comprehensive_analysis()
        save_comprehensive_results(results)
        
        print(f"\n=== ANALYSIS COMPLETED SUCCESSFULLY ===")
        return results
        
    except Exception as e:
        print(f"Analysis failed: {e}")
        return None

if __name__ == "__main__":
    main()