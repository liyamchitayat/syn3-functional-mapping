#!/usr/bin/env python3
"""
Comprehensive BLAST Pipeline for Syn3A Protein Homology Analysis
Performs true sequence homology searches against bacterial proteins
"""

import os
import sys
import json
import csv
import time
import argparse
import subprocess
import tempfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests
from urllib.parse import urlencode

class BlastPipeline:
    def __init__(self, config_file):
        """Initialize BLAST pipeline with configuration"""
        
        with open(config_file, 'r') as f:
            self.config = json.load(f)
        
        self.pipeline_dir = Path(self.config['directories']['pipeline'])
        self.results_dir = Path(self.config['directories']['results'])
        self.logs_dir = Path(self.config['directories']['logs'])
        self.queries_dir = Path(self.config['directories']['queries'])
        
        # Create directories
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        self.log_file = self.logs_dir / f"blast_pipeline_{int(time.time())}.log"
        
        print(f"BLAST Pipeline initialized")
        print(f"Results directory: {self.results_dir}")
        print(f"Log file: {self.log_file}")
    
    def log(self, message):
        """Log message to file and console"""
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        log_entry = f"[{timestamp}] {message}"
        
        print(log_entry)
        
        with open(self.log_file, 'a') as f:
            f.write(log_entry + '\n')
    
    def run_remote_blast(self, query_file, output_file):
        """Run BLAST search via NCBI web service"""
        
        self.log(f"Starting remote BLAST for {query_file}")
        
        try:
            # Read query sequence
            with open(query_file, 'r') as f:
                fasta_content = f.read()
            
            # Extract sequence (remove header lines)
            sequence_lines = [line for line in fasta_content.split('\n') if not line.startswith('>') and line.strip()]
            query_sequence = ''.join(sequence_lines)
            
            if not query_sequence:
                self.log(f"ERROR: No sequence found in {query_file}")
                return False
            
            # Submit BLAST job
            submit_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
            
            params = {
                'CMD': 'Put',
                'PROGRAM': self.config['blast_parameters']['program'],
                'DATABASE': self.config['blast_parameters']['database'],
                'QUERY': query_sequence,
                'FORMAT_TYPE': 'XML',
                'HITLIST_SIZE': self.config['blast_parameters']['max_target_seqs'],
                'EXPECT': self.config['blast_parameters']['evalue'],
                'WORD_SIZE': self.config['blast_parameters']['word_size'],
                'MATRIX': self.config['blast_parameters']['matrix'],
                'ENTREZ_QUERY': 'NOT mycoplasma[Organism]'  # Exclude mycoplasma
            }
            
            # Submit request
            data = urlencode(params).encode('utf-8')
            response = requests.post(submit_url, data=data, timeout=30)
            
            if response.status_code != 200:
                self.log(f"ERROR: BLAST submission failed with status {response.status_code}")
                return False
            
            # Extract RID (Request ID)
            import re
            rid_match = re.search(r'RID = (\w+)', response.text)
            if not rid_match:
                self.log(f"ERROR: Could not extract RID from BLAST response")
                return False
            
            rid = rid_match.group(1)
            self.log(f"BLAST job submitted with RID: {rid}")
            
            # Wait for results
            check_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID={rid}"
            
            max_attempts = 60  # 15 minutes maximum wait
            for attempt in range(max_attempts):
                time.sleep(15)  # Wait 15 seconds between checks
                
                try:
                    result_response = requests.get(check_url, timeout=30)
                    result_text = result_response.text
                    
                    if "Status=WAITING" in result_text:
                        if attempt % 4 == 0:  # Log every minute
                            self.log(f"BLAST {rid} still running... (attempt {attempt + 1}/{max_attempts})")
                        continue
                    elif "Status=FAILED" in result_text:
                        self.log(f"ERROR: BLAST {rid} failed")
                        return False
                    elif "Status=UNKNOWN" in result_text:
                        self.log(f"ERROR: BLAST {rid} RID expired or unknown")
                        return False
                    elif len(result_text) > 1000 and "<BlastOutput>" in result_text:
                        # Results ready
                        self.log(f"BLAST {rid} completed successfully")
                        
                        # Save raw XML results
                        with open(output_file, 'w') as f:
                            f.write(result_text)
                        
                        return True
                
                except Exception as e:
                    self.log(f"ERROR checking BLAST results: {e}")
                    continue
            
            self.log(f"ERROR: BLAST {rid} timed out after {max_attempts} attempts")
            return False
            
        except Exception as e:
            self.log(f"ERROR in remote BLAST: {e}")
            return False
    
    def parse_blast_xml(self, xml_file):
        """Parse BLAST XML results and extract bacterial hits"""
        
        try:
            with open(xml_file, 'r') as f:
                xml_content = f.read()
            
            hits = []
            
            # Simple XML parsing using regex (for robustness without external dependencies)
            import re
            
            # Extract query information
            query_def = re.search(r'<BlastOutput_query-def>(.*?)</BlastOutput_query-def>', xml_content)
            query_len = re.search(r'<BlastOutput_query-len>(\d+)</BlastOutput_query-len>', xml_content)
            
            query_info = {
                'query_def': query_def.group(1) if query_def else 'Unknown',
                'query_len': int(query_len.group(1)) if query_len else 0
            }
            
            # Extract hits
            hit_pattern = r'<Hit>(.*?)</Hit>'
            hit_matches = re.findall(hit_pattern, xml_content, re.DOTALL)
            
            for hit_content in hit_matches:
                # Extract hit information
                hit_def = re.search(r'<Hit_def>(.*?)</Hit_def>', hit_content)
                hit_len = re.search(r'<Hit_len>(\d+)</Hit_len>', hit_content)
                hit_accession = re.search(r'<Hit_accession>(.*?)</Hit_accession>', hit_content)
                
                # Skip mycoplasma hits (double-check)
                if hit_def and 'mycoplasma' in hit_def.group(1).lower():
                    continue
                
                # Extract best HSP (High-scoring Segment Pair)
                hsp_pattern = r'<Hsp>(.*?)</Hsp>'
                hsps = re.findall(hsp_pattern, hit_content, re.DOTALL)
                
                if hsps:
                    best_hsp = hsps[0]  # Take the best (first) HSP
                    
                    # Extract HSP details
                    evalue = re.search(r'<Hsp_evalue>([\d.e-]+)</Hsp_evalue>', best_hsp)
                    identity = re.search(r'<Hsp_identity>(\d+)</Hsp_identity>', best_hsp)
                    positive = re.search(r'<Hsp_positive>(\d+)</Hsp_positive>', best_hsp)
                    align_len = re.search(r'<Hsp_align-len>(\d+)</Hsp_align-len>', best_hsp)
                    bit_score = re.search(r'<Hsp_bit-score>([\d.]+)</Hsp_bit-score>', best_hsp)
                    
                    if evalue and identity and align_len:
                        evalue_val = float(evalue.group(1))
                        identity_val = int(identity.group(1))
                        align_len_val = int(align_len.group(1))
                        
                        # Calculate identity percentage
                        identity_percent = (identity_val / align_len_val) * 100
                        
                        # Calculate coverage percentage
                        coverage_percent = (align_len_val / query_info['query_len']) * 100 if query_info['query_len'] > 0 else 0
                        
                        # Apply filtering criteria
                        if (evalue_val <= self.config['filtering']['max_evalue'] and
                            identity_percent >= self.config['filtering']['min_identity'] and
                            coverage_percent >= self.config['filtering']['min_coverage']):
                            
                            hits.append({
                                'hit_accession': hit_accession.group(1) if hit_accession else 'Unknown',
                                'hit_def': hit_def.group(1) if hit_def else 'Unknown',
                                'hit_len': int(hit_len.group(1)) if hit_len else 0,
                                'evalue': evalue_val,
                                'identity': identity_val,
                                'positive': int(positive.group(1)) if positive else 0,
                                'align_len': align_len_val,
                                'bit_score': float(bit_score.group(1)) if bit_score else 0,
                                'identity_percent': round(identity_percent, 2),
                                'coverage_percent': round(coverage_percent, 2)
                            })
            
            # Sort hits by E-value (best first)
            hits.sort(key=lambda x: x['evalue'])
            
            return query_info, hits
            
        except Exception as e:
            self.log(f"ERROR parsing BLAST XML {xml_file}: {e}")
            return None, []
    
    def process_single_query(self, query_file):
        """Process a single protein query"""
        
        query_name = query_file.stem
        self.log(f"Processing query: {query_name}")
        
        # Setup output files
        xml_output = self.results_dir / f"{query_name}_blast.xml"
        csv_output = self.results_dir / f"{query_name}_hits.csv"
        
        # Run BLAST
        success = self.run_remote_blast(query_file, xml_output)
        
        if not success:
            self.log(f"BLAST failed for {query_name}")
            return {
                'query_name': query_name,
                'query_file': str(query_file),
                'status': 'failed',
                'hit_count': 0,
                'has_bacterial_homologs': False
            }
        
        # Parse results
        query_info, hits = self.parse_blast_xml(xml_output)
        
        if query_info is None:
            self.log(f"Failed to parse results for {query_name}")
            return {
                'query_name': query_name,
                'query_file': str(query_file),
                'status': 'parse_failed',
                'hit_count': 0,
                'has_bacterial_homologs': False
            }
        
        # Save hits to CSV
        if hits:
            with open(csv_output, 'w', newline='') as csvfile:
                fieldnames = [
                    'hit_accession', 'hit_def', 'hit_len', 'evalue', 'identity',
                    'positive', 'align_len', 'bit_score', 'identity_percent', 'coverage_percent'
                ]
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(hits)
        
        result = {
            'query_name': query_name,
            'query_file': str(query_file),
            'query_def': query_info['query_def'],
            'query_len': query_info['query_len'],
            'status': 'success',
            'hit_count': len(hits),
            'has_bacterial_homologs': len(hits) > 0,
            'xml_output': str(xml_output),
            'csv_output': str(csv_output) if hits else None,
            'best_hit_evalue': hits[0]['evalue'] if hits else None,
            'best_hit_identity': hits[0]['identity_percent'] if hits else None,
            'best_hit_description': hits[0]['hit_def'] if hits else None
        }
        
        self.log(f"Completed {query_name}: {len(hits)} bacterial hits found")
        return result
    
    def run_pipeline(self, max_workers=2):
        """Run the complete BLAST pipeline"""
        
        self.log("Starting BLAST pipeline")
        
        # Find query files
        query_files = list(self.queries_dir.glob("*.fasta"))
        
        if not query_files:
            self.log("ERROR: No query files found in queries directory")
            return []
        
        self.log(f"Found {len(query_files)} query files")
        
        # Process queries
        results = []
        
        if max_workers == 1:
            # Sequential processing
            for query_file in query_files:
                result = self.process_single_query(query_file)
                results.append(result)
                
                # Rate limiting between queries
                time.sleep(10)
        else:
            # Parallel processing (limited by NCBI rate limits)
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                # Submit all jobs
                future_to_query = {
                    executor.submit(self.process_single_query, query_file): query_file 
                    for query_file in query_files
                }
                
                # Collect results as they complete
                for future in as_completed(future_to_query):
                    query_file = future_to_query[future]
                    try:
                        result = future.result()
                        results.append(result)
                    except Exception as e:
                        self.log(f"ERROR processing {query_file}: {e}")
                        results.append({
                            'query_name': query_file.stem,
                            'query_file': str(query_file),
                            'status': 'error',
                            'hit_count': 0,
                            'has_bacterial_homologs': False,
                            'error': str(e)
                        })
        
        # Save comprehensive results
        self.save_pipeline_results(results)
        
        return results
    
    def save_pipeline_results(self, results):
        """Save comprehensive pipeline results"""
        
        self.log("Saving pipeline results")
        
        # Save detailed results
        results_file = self.results_dir / "blast_pipeline_results.csv"
        
        fieldnames = [
            'query_name', 'query_file', 'query_def', 'query_len', 'status',
            'hit_count', 'has_bacterial_homologs', 'best_hit_evalue',
            'best_hit_identity', 'best_hit_description', 'xml_output', 'csv_output'
        ]
        
        with open(results_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        
        # Generate summary statistics
        total_queries = len(results)
        successful = len([r for r in results if r['status'] == 'success'])
        with_hits = len([r for r in results if r['has_bacterial_homologs']])
        
        summary = {
            'pipeline_info': self.config['pipeline_info'],
            'analysis_date': time.strftime('%Y-%m-%d %H:%M:%S'),
            'total_queries': total_queries,
            'successful_searches': successful,
            'failed_searches': total_queries - successful,
            'queries_with_bacterial_hits': with_hits,
            'queries_without_bacterial_hits': successful - with_hits,
            'percentage_with_hits': (with_hits / successful * 100) if successful > 0 else 0,
            'blast_parameters': self.config['blast_parameters'],
            'filtering_criteria': self.config['filtering']
        }
        
        # Save summary
        summary_file = self.results_dir / "blast_pipeline_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        self.log(f"Pipeline completed!")
        self.log(f"Total queries: {total_queries}")
        self.log(f"Successful searches: {successful}")
        self.log(f"Queries with bacterial hits: {with_hits}")
        self.log(f"Success rate: {(with_hits / successful * 100):.1f}%" if successful > 0 else "0%")
        self.log(f"Results saved to: {results_file}")
        self.log(f"Summary saved to: {summary_file}")

def main():
    parser = argparse.ArgumentParser(description='Run BLAST pipeline for bacterial homology analysis')
    parser.add_argument('--config', default='blast_config.json', help='Configuration file')
    parser.add_argument('--workers', type=int, default=2, help='Number of parallel workers (max 4 recommended)')
    parser.add_argument('--test', action='store_true', help='Run on first 3 queries only (for testing)')
    
    args = parser.parse_args()
    
    # Check configuration file
    if not os.path.exists(args.config):
        print(f"ERROR: Configuration file {args.config} not found")
        print("Run setup_blast_environment.sh first")
        sys.exit(1)
    
    # Initialize pipeline
    pipeline = BlastPipeline(args.config)
    
    # Limit workers for NCBI rate limiting
    max_workers = min(args.workers, 4)
    
    if args.test:
        # Test mode - process only first few queries
        query_files = list(pipeline.queries_dir.glob("*.fasta"))[:3]
        if query_files:
            pipeline.log(f"TEST MODE: Processing {len(query_files)} queries")
            results = []
            for query_file in query_files:
                result = pipeline.process_single_query(query_file)
                results.append(result)
                time.sleep(30)  # Extra delay for testing
            
            pipeline.save_pipeline_results(results)
        else:
            pipeline.log("ERROR: No query files found for testing")
    else:
        # Full pipeline
        results = pipeline.run_pipeline(max_workers=max_workers)
    
    print(f"\nBLAST pipeline completed!")
    print(f"Check results in: {pipeline.results_dir}")
    print(f"Check logs in: {pipeline.log_file}")

if __name__ == "__main__":
    main()