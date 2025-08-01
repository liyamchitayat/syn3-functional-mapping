#!/usr/bin/env python3
"""
Parallel BLAST Pipeline - Optimized for Speed
Runs multiple BLAST searches concurrently with intelligent batching
"""

import os
import sys
import json
import csv
import time
import re
import asyncio
import aiohttp
import concurrent.futures
from pathlib import Path
from urllib.parse import urlencode
import threading
import queue
import signal

class ParallelBlastPipeline:
    def __init__(self, output_dir="parallel_blast_results", max_concurrent=8):
        """Initialize parallel BLAST pipeline"""
        
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.results_dir = self.output_dir / "results"
        self.logs_dir = self.output_dir / "logs"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        
        # Performance parameters
        self.max_concurrent = min(max_concurrent, 8)  # NCBI rate limit compliance
        self.submission_delay = 2  # Seconds between submissions
        self.check_interval = 10   # Seconds between status checks
        self.max_attempts = 60     # Maximum wait attempts (10 minutes)
        
        # Progress tracking
        self.completed_count = 0
        self.total_count = 0
        self.start_time = time.time()
        self.results_lock = threading.Lock()
        
        # Results storage
        self.all_results = []
        
        self.log_file = self.logs_dir / f"parallel_blast_{int(time.time())}.log"
        
        print(f"🚀 Parallel BLAST Pipeline initialized")
        print(f"   Max concurrent jobs: {self.max_concurrent}")
        print(f"   Results directory: {self.results_dir}")
        print(f"   Log file: {self.log_file}")
    
    def log(self, message):
        """Thread-safe logging"""
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        log_entry = f"[{timestamp}] {message}"
        
        print(log_entry)
        
        with open(self.log_file, 'a') as f:
            f.write(log_entry + '\n')
    
    def parse_fasta_file(self, fasta_file):
        """Parse FASTA file into protein list"""
        proteins = []
        
        with open(fasta_file, 'r') as f:
            current_header = None
            current_sequence = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_header:
                        proteins.append({
                            'header': current_header,
                            'sequence': ''.join(current_sequence),
                            'protein_id': self.extract_protein_id(current_header),
                            'protein_name': self.extract_protein_name(current_header),
                            'index': len(proteins)
                        })
                    
                    current_header = line
                    current_sequence = []
                else:
                    current_sequence.append(line)
            
            if current_header:
                proteins.append({
                    'header': current_header,
                    'sequence': ''.join(current_sequence),
                    'protein_id': self.extract_protein_id(current_header),
                    'protein_name': self.extract_protein_name(current_header),
                    'index': len(proteins)
                })
        
        return proteins
    
    def extract_protein_id(self, header):
        """Extract protein ID from header"""
        match = re.search(r'protein_id=([^\]]+)', header)
        return match.group(1) if match else f"protein_{hash(header) % 10000}"
    
    def extract_protein_name(self, header):
        """Extract protein name from header"""
        match = re.search(r'protein=([^\]]+)', header)
        return match.group(1) if match else "hypothetical protein"
    
    def submit_blast_optimized(self, sequence, protein_id):
        """Submit BLAST with optimized parameters for speed"""
        
        try:
            submit_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
            
            # Optimized parameters for speed
            params = {
                'CMD': 'Put',
                'PROGRAM': 'blastp',
                'DATABASE': 'nr',
                'QUERY': sequence,
                'FORMAT_TYPE': 'XML',
                'HITLIST_SIZE': 20,        # Fewer hits for speed
                'EXPECT': 0.001,           # Stricter E-value for speed
                'WORD_SIZE': 6,            # Larger word size for speed
                'MATRIX': 'BLOSUM62',
                'COMPOSITION_BASED_STATISTICS': '0',  # Disable for speed
                'FILTER': 'F',             # No low complexity filter for speed
                'ENTREZ_QUERY': 'NOT mycoplasma[Organism]'
            }
            
            data = urlencode(params).encode('utf-8')
            
            import requests
            response = requests.post(submit_url, data=data, timeout=20)
            
            if response.status_code != 200:
                return None, f"HTTP error {response.status_code}"
            
            # Extract RID
            rid_match = re.search(r'RID = (\w+)', response.text)
            if not rid_match:
                return None, "No RID found"
            
            return rid_match.group(1), "success"
            
        except Exception as e:
            return None, str(e)
    
    def check_blast_status_fast(self, rid):
        """Fast BLAST status checking"""
        
        try:
            check_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID={rid}"
            
            import requests
            response = requests.get(check_url, timeout=15)
            
            if response.status_code != 200:
                return None, "HTTP error"
            
            result_text = response.text
            
            if "Status=WAITING" in result_text:
                return None, "waiting"
            elif "Status=FAILED" in result_text:
                return None, "failed"
            elif "Status=UNKNOWN" in result_text:
                return None, "expired"
            elif len(result_text) > 1000 and "<BlastOutput>" in result_text:
                return result_text, "ready"
            else:
                return None, "unknown"
                
        except Exception as e:
            return None, str(e)
    
    def parse_blast_xml_fast(self, xml_content, query_length):
        """Fast XML parsing focused on key metrics"""
        
        try:
            hits = []
            
            # Extract hits using optimized regex
            hit_pattern = r'<Hit>(.*?)</Hit>'
            hit_matches = re.findall(hit_pattern, xml_content, re.DOTALL)
            
            for hit_content in hit_matches[:10]:  # Only process top 10 hits for speed
                # Skip mycoplasma (double-check)
                hit_def = re.search(r'<Hit_def>(.*?)</Hit_def>', hit_content)
                if hit_def and 'mycoplasma' in hit_def.group(1).lower():
                    continue
                
                # Extract key HSP data
                evalue = re.search(r'<Hsp_evalue>([\d.e-]+)</Hsp_evalue>', hit_content)
                identity = re.search(r'<Hsp_identity>(\d+)</Hsp_identity>', hit_content)
                align_len = re.search(r'<Hsp_align-len>(\d+)</Hsp_align-len>', hit_content)
                
                if evalue and identity and align_len:
                    evalue_val = float(evalue.group(1))
                    identity_val = int(identity.group(1))
                    align_len_val = int(align_len.group(1))
                    
                    identity_percent = (identity_val / align_len_val) * 100
                    coverage_percent = (align_len_val / query_length) * 100 if query_length > 0 else 0
                    
                    # Quick filtering
                    if evalue_val <= 0.01 and identity_percent >= 30.0 and coverage_percent >= 30.0:
                        hit_accession = re.search(r'<Hit_accession>(.*?)</Hit_accession>', hit_content)
                        
                        hits.append({
                            'hit_accession': hit_accession.group(1) if hit_accession else 'Unknown',
                            'hit_def': hit_def.group(1) if hit_def else 'Unknown',
                            'evalue': evalue_val,
                            'identity_percent': round(identity_percent, 2),
                            'coverage_percent': round(coverage_percent, 2)
                        })
            
            return hits
            
        except Exception as e:
            self.log(f"ERROR parsing XML: {e}")
            return []
    
    def process_single_protein_worker(self, protein):
        """Worker function for processing a single protein"""
        
        protein_id = protein['protein_id']
        sequence = protein['sequence']
        
        try:
            # Submit BLAST
            rid, status = self.submit_blast_optimized(sequence, protein_id)
            
            if not rid:
                return self.create_failed_result(protein, f"submit_failed: {status}")
            
            # Wait for results with timeout
            xml_result = None
            for attempt in range(self.max_attempts):
                time.sleep(self.check_interval)
                
                xml_result, check_status = self.check_blast_status_fast(rid)
                
                if check_status == "ready":
                    break
                elif check_status in ["failed", "expired"]:
                    return self.create_failed_result(protein, f"blast_{check_status}")
                elif check_status == "waiting":
                    continue
                else:
                    # Unknown status, continue trying
                    continue
            
            if not xml_result:
                return self.create_failed_result(protein, "timeout")
            
            # Save raw XML (optional, for debugging)
            xml_file = self.results_dir / f"{protein_id}_blast.xml"
            with open(xml_file, 'w') as f:
                f.write(xml_result)
            
            # Parse results
            hits = self.parse_blast_xml_fast(xml_result, len(sequence))
            
            # Create result
            result = {
                'protein_id': protein_id,
                'protein_name': protein['protein_name'],
                'sequence_length': len(sequence),
                'blast_rid': rid,
                'status': 'success',
                'bacterial_hits': len(hits),
                'has_bacterial_homologs': len(hits) > 0,
                'best_hit_evalue': hits[0]['evalue'] if hits else None,
                'best_hit_identity': hits[0]['identity_percent'] if hits else None,
                'best_hit_coverage': hits[0]['coverage_percent'] if hits else None,
                'best_hit_description': hits[0]['hit_def'][:100] if hits else None,
                'processing_time': time.strftime('%Y-%m-%d %H:%M:%S')
            }
            
            # Update progress
            with self.results_lock:
                self.completed_count += 1
                elapsed = time.time() - self.start_time
                rate = self.completed_count / elapsed * 60  # proteins per minute
                remaining = (self.total_count - self.completed_count) / (self.completed_count / elapsed) / 3600  # hours
                
                self.log(f"✅ Completed {protein_id} ({self.completed_count}/{self.total_count}) - "
                        f"Rate: {rate:.1f}/min, ETA: {remaining:.1f}h, Hits: {len(hits)}")
            
            return result
            
        except Exception as e:
            self.log(f"❌ Error processing {protein_id}: {e}")
            return self.create_failed_result(protein, f"error: {str(e)}")
    
    def create_failed_result(self, protein, status):
        """Create failed result"""
        return {
            'protein_id': protein['protein_id'],
            'protein_name': protein['protein_name'],
            'sequence_length': len(protein['sequence']),
            'blast_rid': None,
            'status': status,
            'bacterial_hits': 0,
            'has_bacterial_homologs': False,
            'best_hit_evalue': None,
            'best_hit_identity': None,
            'best_hit_coverage': None,
            'best_hit_description': None,
            'processing_time': time.strftime('%Y-%m-%d %H:%M:%S')
        }
    
    def run_parallel_analysis(self, fasta_file, start_protein=0, max_proteins=None):
        """Run parallel BLAST analysis"""
        
        self.log(f"🚀 Starting parallel BLAST analysis")
        self.log(f"Input file: {fasta_file}")
        self.log(f"Max concurrent jobs: {self.max_concurrent}")
        
        # Parse proteins
        proteins = self.parse_fasta_file(fasta_file)
        self.log(f"Found {len(proteins)} proteins")
        
        # Select range
        if max_proteins:
            proteins = proteins[start_protein:start_protein + max_proteins]
        else:
            proteins = proteins[start_protein:]
        
        self.total_count = len(proteins)
        self.log(f"Processing {self.total_count} proteins (starting from {start_protein})")
        
        # Process proteins in parallel
        self.all_results = []
        
        # Use ThreadPoolExecutor for I/O bound tasks
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_concurrent) as executor:
            # Submit all jobs with staggered start times
            futures = []
            
            for i, protein in enumerate(proteins):
                # Stagger submissions to avoid overwhelming NCBI
                delay = i * self.submission_delay
                
                def delayed_submit(p, d):
                    time.sleep(d)
                    return self.process_single_protein_worker(p)
                
                future = executor.submit(delayed_submit, protein, delay)
                futures.append(future)
            
            # Collect results as they complete
            try:
                for future in concurrent.futures.as_completed(futures, timeout=3600*5):  # 5 hour timeout
                    try:
                        result = future.result()
                        self.all_results.append(result)
                    except Exception as e:
                        self.log(f"❌ Future error: {e}")
                        
            except KeyboardInterrupt:
                self.log("🛑 Analysis interrupted by user")
                executor.shutdown(wait=False)
            except concurrent.futures.TimeoutError:
                self.log("⏰ Analysis timed out")
                executor.shutdown(wait=False)
        
        # Save results
        self.save_parallel_results()
        
        return self.all_results
    
    def save_parallel_results(self):
        """Save comprehensive parallel results"""
        
        self.log("💾 Saving parallel analysis results...")
        
        # Save detailed results
        results_file = self.results_dir / "parallel_blast_results.csv"
        
        fieldnames = [
            'protein_id', 'protein_name', 'sequence_length', 'blast_rid',
            'status', 'bacterial_hits', 'has_bacterial_homologs',
            'best_hit_evalue', 'best_hit_identity', 'best_hit_coverage',
            'best_hit_description', 'processing_time'
        ]
        
        with open(results_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(self.all_results)
        
        # Generate summary
        total = len(self.all_results)
        successful = len([r for r in self.all_results if r['status'] == 'success'])
        with_hits = len([r for r in self.all_results if r['has_bacterial_homologs']])
        
        elapsed = time.time() - self.start_time
        rate = total / elapsed * 3600  # proteins per hour
        
        summary = {
            'analysis_date': time.strftime('%Y-%m-%d %H:%M:%S'),
            'analysis_type': 'parallel_blast',
            'max_concurrent_jobs': self.max_concurrent,
            'total_proteins': total,
            'successful_searches': successful,
            'failed_searches': total - successful,
            'proteins_with_hits': with_hits,
            'proteins_without_hits': successful - with_hits,
            'success_rate': (successful / total * 100) if total > 0 else 0,
            'homolog_detection_rate': (with_hits / successful * 100) if successful > 0 else 0,
            'total_time_hours': elapsed / 3600,
            'processing_rate_per_hour': rate,
            'average_time_per_protein_minutes': (elapsed / total * 60) if total > 0 else 0
        }
        
        summary_file = self.results_dir / "parallel_blast_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Log final statistics
        self.log(f"🎉 Parallel analysis completed!")
        self.log(f"   Total proteins: {total}")
        self.log(f"   Successful: {successful} ({summary['success_rate']:.1f}%)")
        self.log(f"   With bacterial hits: {with_hits} ({summary['homolog_detection_rate']:.1f}%)")
        self.log(f"   Processing rate: {rate:.1f} proteins/hour")
        self.log(f"   Total time: {elapsed/3600:.2f} hours")
        self.log(f"   Results: {results_file}")
        self.log(f"   Summary: {summary_file}")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Parallel BLAST pipeline for maximum speed')
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('--output-dir', default='parallel_blast_results', help='Output directory')
    parser.add_argument('--start', type=int, default=0, help='Start protein index')
    parser.add_argument('--max', type=int, help='Maximum proteins to process')
    parser.add_argument('--concurrent', type=int, default=6, help='Max concurrent jobs (1-8)')
    parser.add_argument('--test', action='store_true', help='Test mode - 5 proteins only')
    
    args = parser.parse_args()
    
    # Validate concurrent jobs
    args.concurrent = max(1, min(args.concurrent, 8))
    
    if not os.path.exists(args.fasta_file):
        print(f"❌ ERROR: FASTA file {args.fasta_file} not found")
        sys.exit(1)
    
    # Initialize pipeline
    pipeline = ParallelBlastPipeline(args.output_dir, args.concurrent)
    
    # Set parameters
    if args.test:
        max_proteins = 5
        start_protein = 0
        print("🧪 TEST MODE: Processing 5 proteins with parallel pipeline")
    else:
        max_proteins = args.max
        start_protein = args.start
    
    # Run analysis
    try:
        results = pipeline.run_parallel_analysis(args.fasta_file, start_protein, max_proteins)
        
        print(f"\n🎉 Parallel BLAST analysis completed!")
        print(f"📁 Results directory: {pipeline.results_dir}")
        print(f"📄 Log file: {pipeline.log_file}")
        
    except KeyboardInterrupt:
        print("\n🛑 Analysis interrupted by user")
    except Exception as e:
        print(f"❌ Analysis failed: {e}")

if __name__ == "__main__":
    main()