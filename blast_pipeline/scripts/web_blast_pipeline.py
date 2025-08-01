#!/usr/bin/env python3
"""
Web-based BLAST Pipeline for Syn3A Protein Homology Analysis
Uses NCBI web services - no local BLAST installation required
"""

import os
import sys
import json
import csv
import time
import re
from pathlib import Path
import requests
from urllib.parse import urlencode

class WebBlastPipeline:
    def __init__(self, output_dir="web_blast_results"):
        """Initialize web-based BLAST pipeline"""
        
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.results_dir = self.output_dir / "results"
        self.logs_dir = self.output_dir / "logs"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        
        self.log_file = self.logs_dir / f"web_blast_{int(time.time())}.log"
        
        # BLAST configuration
        self.config = {
            'evalue_threshold': 0.01,
            'max_hits': 50,
            'min_identity': 30.0,
            'min_coverage': 50.0,
            'exclude_organisms': ['mycoplasma']
        }
        
        print(f"Web BLAST Pipeline initialized")
        print(f"Output directory: {self.output_dir}")
        print(f"Results directory: {self.results_dir}")
        print(f"Log file: {self.log_file}")
    
    def log(self, message):
        """Log message with timestamp"""
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        log_entry = f"[{timestamp}] {message}"
        
        print(log_entry)
        
        with open(self.log_file, 'a') as f:
            f.write(log_entry + '\n')
    
    def parse_fasta_file(self, fasta_file):
        """Parse FASTA file and return individual proteins"""
        proteins = []
        
        with open(fasta_file, 'r') as f:
            current_header = None
            current_sequence = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous protein
                    if current_header:
                        proteins.append({
                            'header': current_header,
                            'sequence': ''.join(current_sequence),
                            'protein_id': self.extract_protein_id(current_header),
                            'protein_name': self.extract_protein_name(current_header)
                        })
                    
                    current_header = line
                    current_sequence = []
                else:
                    current_sequence.append(line)
            
            # Add last protein
            if current_header:
                proteins.append({
                    'header': current_header,
                    'sequence': ''.join(current_sequence),
                    'protein_id': self.extract_protein_id(current_header),
                    'protein_name': self.extract_protein_name(current_header)
                })
        
        return proteins
    
    def extract_protein_id(self, header):
        """Extract protein ID from FASTA header"""
        match = re.search(r'protein_id=([^\]]+)', header)
        return match.group(1) if match else "Unknown"
    
    def extract_protein_name(self, header):
        """Extract protein name from FASTA header"""
        match = re.search(r'protein=([^\]]+)', header)
        return match.group(1) if match else "hypothetical protein"
    
    def submit_blast(self, sequence, protein_id):
        """Submit BLAST search to NCBI"""
        
        self.log(f"Submitting BLAST for {protein_id}")
        
        try:
            submit_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
            
            params = {
                'CMD': 'Put',
                'PROGRAM': 'blastp',
                'DATABASE': 'nr',
                'QUERY': sequence,
                'FORMAT_TYPE': 'XML',
                'HITLIST_SIZE': self.config['max_hits'],
                'EXPECT': self.config['evalue_threshold'],
                'WORD_SIZE': 6,
                'MATRIX': 'BLOSUM62',
                'ENTREZ_QUERY': 'NOT mycoplasma[Organism]'
            }
            
            data = urlencode(params).encode('utf-8')
            response = requests.post(submit_url, data=data, timeout=30)
            
            if response.status_code != 200:
                self.log(f"ERROR: BLAST submission failed with status {response.status_code}")
                return None
            
            # Extract RID
            rid_match = re.search(r'RID = (\w+)', response.text)
            if not rid_match:
                self.log(f"ERROR: Could not extract RID from BLAST response")
                return None
            
            rid = rid_match.group(1)
            self.log(f"BLAST submitted successfully: RID = {rid}")
            return rid
            
        except Exception as e:
            self.log(f"ERROR submitting BLAST: {e}")
            return None
    
    def check_blast_status(self, rid, max_attempts=40):
        """Check BLAST status and retrieve results when ready"""
        
        check_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID={rid}"
        
        for attempt in range(max_attempts):
            try:
                time.sleep(15 + (attempt * 2))  # Progressive delay
                
                response = requests.get(check_url, timeout=30)
                result_text = response.text
                
                if "Status=WAITING" in result_text:
                    if attempt % 4 == 0:  # Log every minute
                        self.log(f"BLAST {rid} still running... (attempt {attempt + 1}/{max_attempts})")
                    continue
                elif "Status=FAILED" in result_text:
                    self.log(f"ERROR: BLAST {rid} failed")
                    return None
                elif "Status=UNKNOWN" in result_text:
                    self.log(f"ERROR: BLAST {rid} expired or unknown")
                    return None
                elif len(result_text) > 1000 and "<BlastOutput>" in result_text:
                    self.log(f"BLAST {rid} completed successfully")
                    return result_text
                
            except Exception as e:
                self.log(f"ERROR checking BLAST status: {e}")
                continue
        
        self.log(f"ERROR: BLAST {rid} timed out after {max_attempts} attempts")
        return None
    
    def parse_blast_xml(self, xml_content, query_length):
        """Parse BLAST XML and extract bacterial hits"""
        
        try:
            hits = []
            
            # Extract hits using regex
            hit_pattern = r'<Hit>(.*?)</Hit>'
            hit_matches = re.findall(hit_pattern, xml_content, re.DOTALL)
            
            for hit_content in hit_matches:
                # Extract hit information
                hit_def = re.search(r'<Hit_def>(.*?)</Hit_def>', hit_content)
                hit_accession = re.search(r'<Hit_accession>(.*?)</Hit_accession>', hit_content)
                hit_len = re.search(r'<Hit_len>(\d+)</Hit_len>', hit_content)
                
                # Skip mycoplasma hits (double-check)
                if hit_def and any(org in hit_def.group(1).lower() for org in self.config['exclude_organisms']):
                    continue
                
                # Extract best HSP
                hsp_pattern = r'<Hsp>(.*?)</Hsp>'
                hsps = re.findall(hsp_pattern, hit_content, re.DOTALL)
                
                if hsps:
                    best_hsp = hsps[0]
                    
                    # Extract HSP details
                    evalue = re.search(r'<Hsp_evalue>([\d.e-]+)</Hsp_evalue>', best_hsp)
                    identity = re.search(r'<Hsp_identity>(\d+)</Hsp_identity>', best_hsp)
                    align_len = re.search(r'<Hsp_align-len>(\d+)</Hsp_align-len>', best_hsp)
                    bit_score = re.search(r'<Hsp_bit-score>([\d.]+)</Hsp_bit-score>', best_hsp)
                    
                    if evalue and identity and align_len:
                        evalue_val = float(evalue.group(1))
                        identity_val = int(identity.group(1))
                        align_len_val = int(align_len.group(1))
                        
                        # Calculate percentages
                        identity_percent = (identity_val / align_len_val) * 100
                        coverage_percent = (align_len_val / query_length) * 100 if query_length > 0 else 0
                        
                        # Apply filters
                        if (evalue_val <= self.config['evalue_threshold'] and
                            identity_percent >= self.config['min_identity'] and
                            coverage_percent >= self.config['min_coverage']):
                            
                            hits.append({
                                'hit_accession': hit_accession.group(1) if hit_accession else 'Unknown',
                                'hit_def': hit_def.group(1) if hit_def else 'Unknown',
                                'hit_len': int(hit_len.group(1)) if hit_len else 0,
                                'evalue': evalue_val,
                                'identity': identity_val,
                                'align_len': align_len_val,
                                'bit_score': float(bit_score.group(1)) if bit_score else 0,
                                'identity_percent': round(identity_percent, 2),
                                'coverage_percent': round(coverage_percent, 2)
                            })
            
            # Sort by E-value (best first)
            hits.sort(key=lambda x: x['evalue'])
            return hits
            
        except Exception as e:
            self.log(f"ERROR parsing BLAST XML: {e}")
            return []
    
    def analyze_protein(self, protein, index, total):
        """Analyze a single protein"""
        
        protein_id = protein['protein_id']
        sequence = protein['sequence']
        
        self.log(f"Analyzing protein {index + 1}/{total}: {protein_id}")
        self.log(f"Protein name: {protein['protein_name'][:80]}...")
        self.log(f"Sequence length: {len(sequence)} amino acids")
        
        # Submit BLAST
        rid = self.submit_blast(sequence, protein_id)
        if not rid:
            return self.create_failed_result(protein, "submission_failed")
        
        # Wait for results
        xml_result = self.check_blast_status(rid)
        if not xml_result:
            return self.create_failed_result(protein, "retrieval_failed")
        
        # Save raw XML
        xml_file = self.results_dir / f"{protein_id}_blast.xml"
        with open(xml_file, 'w') as f:
            f.write(xml_result)
        
        # Parse results
        hits = self.parse_blast_xml(xml_result, len(sequence))
        
        # Save hits
        if hits:
            csv_file = self.results_dir / f"{protein_id}_hits.csv"
            with open(csv_file, 'w', newline='') as csvfile:
                fieldnames = [
                    'hit_accession', 'hit_def', 'hit_len', 'evalue',
                    'identity', 'align_len', 'bit_score', 
                    'identity_percent', 'coverage_percent'
                ]
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(hits)
        
        result = {
            'protein_id': protein_id,
            'protein_name': protein['protein_name'],
            'sequence_length': len(sequence),
            'blast_rid': rid,
            'status': 'success',
            'bacterial_hits': len(hits),
            'has_bacterial_homologs': len(hits) > 0,
            'xml_file': str(xml_file),
            'csv_file': str(csv_file) if hits else None,
            'best_hit_evalue': hits[0]['evalue'] if hits else None,
            'best_hit_identity': hits[0]['identity_percent'] if hits else None,
            'best_hit_description': hits[0]['hit_def'][:100] if hits else None,
            'analysis_time': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        self.log(f"Completed {protein_id}: {len(hits)} bacterial hits found")
        return result
    
    def create_failed_result(self, protein, status):
        """Create result for failed analysis"""
        return {
            'protein_id': protein['protein_id'],
            'protein_name': protein['protein_name'],
            'sequence_length': len(protein['sequence']),
            'blast_rid': None,
            'status': status,
            'bacterial_hits': 0,
            'has_bacterial_homologs': False,
            'xml_file': None,
            'csv_file': None,
            'best_hit_evalue': None,
            'best_hit_identity': None,
            'best_hit_description': None,
            'analysis_time': time.strftime('%Y-%m-%d %H:%M:%S')
        }
    
    def run_analysis(self, fasta_file, start_protein=0, max_proteins=None):
        """Run complete analysis"""
        
        self.log(f"Starting web BLAST analysis")
        self.log(f"Input file: {fasta_file}")
        
        # Parse proteins
        proteins = self.parse_fasta_file(fasta_file)
        self.log(f"Found {len(proteins)} proteins in FASTA file")
        
        # Determine range
        if max_proteins:
            proteins = proteins[start_protein:start_protein + max_proteins]
        else:
            proteins = proteins[start_protein:]
        
        self.log(f"Analyzing {len(proteins)} proteins (starting from protein {start_protein})")
        
        # Process proteins
        results = []
        
        for i, protein in enumerate(proteins):
            try:
                result = self.analyze_protein(protein, i, len(proteins))
                results.append(result)
                
                # Rate limiting - wait between proteins
                if i < len(proteins) - 1:
                    self.log("Waiting 30 seconds before next protein...")
                    time.sleep(30)
                
            except KeyboardInterrupt:
                self.log("Analysis interrupted by user")
                break
            except Exception as e:
                self.log(f"ERROR analyzing {protein['protein_id']}: {e}")
                results.append(self.create_failed_result(protein, f"error: {str(e)}"))
        
        # Save results
        self.save_results(results)
        
        return results
    
    def save_results(self, results):
        """Save comprehensive results"""
        
        self.log("Saving results...")
        
        # Save detailed results
        results_file = self.results_dir / "web_blast_results.csv"
        
        fieldnames = [
            'protein_id', 'protein_name', 'sequence_length', 'blast_rid',
            'status', 'bacterial_hits', 'has_bacterial_homologs',
            'best_hit_evalue', 'best_hit_identity', 'best_hit_description',
            'xml_file', 'csv_file', 'analysis_time'
        ]
        
        with open(results_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        
        # Summary statistics
        total = len(results)
        successful = len([r for r in results if r['status'] == 'success'])
        with_hits = len([r for r in results if r['has_bacterial_homologs']])
        
        summary = {
            'analysis_date': time.strftime('%Y-%m-%d %H:%M:%S'),
            'total_proteins': total,
            'successful_searches': successful,
            'failed_searches': total - successful,
            'proteins_with_hits': with_hits,
            'proteins_without_hits': successful - with_hits,
            'success_rate': (successful / total * 100) if total > 0 else 0,
            'homolog_rate': (with_hits / successful * 100) if successful > 0 else 0,
            'configuration': self.config
        }
        
        summary_file = self.results_dir / "web_blast_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        self.log(f"Analysis complete!")
        self.log(f"Total proteins: {total}")
        self.log(f"Successful searches: {successful}")
        self.log(f"Proteins with bacterial hits: {with_hits}")
        self.log(f"Success rate: {summary['success_rate']:.1f}%")
        self.log(f"Homolog detection rate: {summary['homolog_rate']:.1f}%")
        self.log(f"Results saved to: {results_file}")
        self.log(f"Summary saved to: {summary_file}")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Web-based BLAST pipeline')
    parser.add_argument('fasta_file', help='Input FASTA file with protein sequences')
    parser.add_argument('--output-dir', default='web_blast_results', help='Output directory')
    parser.add_argument('--start', type=int, default=0, help='Start protein index (0-based)')
    parser.add_argument('--max', type=int, help='Maximum number of proteins to analyze')
    parser.add_argument('--test', action='store_true', help='Test mode - analyze first 3 proteins')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.fasta_file):
        print(f"ERROR: FASTA file {args.fasta_file} not found")
        sys.exit(1)
    
    # Initialize pipeline
    pipeline = WebBlastPipeline(args.output_dir)
    
    # Determine parameters
    if args.test:
        max_proteins = 3
        start_protein = 0
        print("TEST MODE: Analyzing first 3 proteins")
    else:
        max_proteins = args.max
        start_protein = args.start
    
    # Run analysis
    try:
        results = pipeline.run_analysis(args.fasta_file, start_protein, max_proteins)
        
        print(f"\nWeb BLAST analysis completed!")
        print(f"Results directory: {pipeline.results_dir}")
        print(f"Log file: {pipeline.log_file}")
        
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user")
    except Exception as e:
        print(f"Analysis failed: {e}")

if __name__ == "__main__":
    main()