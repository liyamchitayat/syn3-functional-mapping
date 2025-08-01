#!/usr/bin/env python3
"""
Fast BLAST Pipeline - Optimized for Speed (Standard Libraries Only)
Uses threading and optimized parameters for maximum performance
"""

import os
import sys
import json
import csv
import time
import re
import threading
import concurrent.futures
from pathlib import Path
from urllib.parse import urlencode
import requests
import queue

class FastBlastPipeline:
    def __init__(self, output_dir="fast_blast_results", max_workers=6):
        """Initialize fast BLAST pipeline"""
        
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.results_dir = self.output_dir / "results"
        self.logs_dir = self.output_dir / "logs"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        
        # Optimized parameters
        self.max_workers = min(max_workers, 8)  # NCBI compliance
        self.submission_delay = 1.5  # Faster submission rate
        self.check_interval = 8     # Faster checking
        self.max_wait_time = 600    # 10 minutes max per protein
        
        # Performance tracking
        self.progress_lock = threading.Lock()
        self.completed = 0
        self.total = 0
        self.start_time = time.time()
        self.results = []
        
        self.log_file = self.logs_dir / f"fast_blast_{int(time.time())}.log"
        
        print(f"🚀 Fast BLAST Pipeline initialized")
        print(f"   Workers: {self.max_workers}")
        print(f"   Results: {self.results_dir}")
        print(f"   Logs: {self.log_file}")
    
    def log(self, message):
        """Thread-safe logging"""
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        log_entry = f"[{timestamp}] {message}"
        
        print(log_entry)
        with open(self.log_file, 'a') as f:
            f.write(log_entry + '\n')
    
    def parse_proteins(self, fasta_file):
        """Parse FASTA file efficiently"""
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
                            'protein_name': self.extract_protein_name(current_header)
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
                    'protein_name': self.extract_protein_name(current_header)
                })
        
        return proteins
    
    def extract_protein_id(self, header):
        """Extract protein ID"""
        match = re.search(r'protein_id=([^\]]+)', header)
        return match.group(1) if match else f"protein_{abs(hash(header)) % 10000}"
    
    def extract_protein_name(self, header):
        """Extract protein name"""
        match = re.search(r'protein=([^\]]+)', header)
        return match.group(1) if match else "hypothetical protein"
    
    def submit_blast_fast(self, sequence, protein_id):
        """Submit BLAST with speed-optimized parameters"""
        
        try:
            submit_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
            
            # Speed-optimized parameters
            params = {
                'CMD': 'Put',
                'PROGRAM': 'blastp',
                'DATABASE': 'nr',
                'QUERY': sequence,
                'FORMAT_TYPE': 'XML',
                'HITLIST_SIZE': 15,           # Fewer hits for speed
                'EXPECT': 0.001,              # Stricter for speed
                'WORD_SIZE': 6,               # Larger word size
                'MATRIX': 'BLOSUM62',
                'COMPOSITION_BASED_STATISTICS': '0',  # Disable
                'FILTER': 'F',                # No filtering
                'ENTREZ_QUERY': 'NOT mycoplasma[Organism]'
            }
            
            data = urlencode(params).encode('utf-8')
            response = requests.post(submit_url, data=data, timeout=15)
            
            if response.status_code != 200:
                return None, f"HTTP_{response.status_code}"
            
            rid_match = re.search(r'RID = (\w+)', response.text)
            return rid_match.group(1) if rid_match else None, "success"
            
        except Exception as e:
            return None, str(e)
    
    def check_blast_fast(self, rid):
        """Fast BLAST status checking"""
        
        try:
            url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID={rid}"
            response = requests.get(url, timeout=10)
            
            if response.status_code != 200:
                return None, "http_error"
            
            text = response.text
            
            if "Status=WAITING" in text:
                return None, "waiting"
            elif "Status=FAILED" in text:
                return None, "failed"
            elif "Status=UNKNOWN" in text:
                return None, "expired"
            elif len(text) > 1000 and "<BlastOutput>" in text:
                return text, "ready"
            else:
                return None, "unknown"
                
        except Exception as e:
            return None, str(e)
    
    def parse_xml_fast(self, xml_content, query_length):
        """Fast XML parsing for key results"""
        
        try:
            hits = []
            
            # Fast regex-based parsing
            hit_pattern = r'<Hit>(.*?)</Hit>'
            matches = re.findall(hit_pattern, xml_content, re.DOTALL)
            
            for hit_content in matches[:8]:  # Top 8 hits only
                # Skip mycoplasma
                hit_def = re.search(r'<Hit_def>(.*?)</Hit_def>', hit_content)
                if hit_def and 'mycoplasma' in hit_def.group(1).lower():
                    continue
                
                # Extract key data
                evalue = re.search(r'<Hsp_evalue>([\d.e-]+)</Hsp_evalue>', hit_content)
                identity = re.search(r'<Hsp_identity>(\d+)</Hsp_identity>', hit_content)
                align_len = re.search(r'<Hsp_align-len>(\d+)</Hsp_align-len>', hit_content)
                
                if evalue and identity and align_len:
                    e_val = float(evalue.group(1))
                    ident = int(identity.group(1))
                    alen = int(align_len.group(1))
                    
                    identity_pct = (ident / alen) * 100
                    coverage_pct = (alen / query_length) * 100 if query_length > 0 else 0
                    
                    # Quick filtering
                    if e_val <= 0.01 and identity_pct >= 25.0 and coverage_pct >= 25.0:
                        hits.append({
                            'evalue': e_val,
                            'identity_percent': round(identity_pct, 1),
                            'coverage_percent': round(coverage_pct, 1),
                            'description': hit_def.group(1)[:80] if hit_def else 'Unknown'
                        })
            
            return hits
            
        except Exception as e:
            return []
    
    def process_protein_worker(self, protein):
        """Worker function for processing single protein"""
        
        protein_id = protein['protein_id']
        sequence = protein['sequence']
        
        try:
            # Submit BLAST
            rid, status = self.submit_blast_fast(sequence, protein_id)
            
            if not rid:
                return self.create_result(protein, f"submit_failed_{status}", 0, [])
            
            # Wait for results
            start_wait = time.time()
            xml_result = None
            
            while (time.time() - start_wait) < self.max_wait_time:
                time.sleep(self.check_interval)
                
                xml_result, check_status = self.check_blast_fast(rid)
                
                if check_status == "ready":
                    break
                elif check_status in ["failed", "expired"]:
                    return self.create_result(protein, f"blast_{check_status}", 0, [])
                # Continue waiting for "waiting" or other statuses
            
            if not xml_result:
                return self.create_result(protein, "timeout", 0, [])
            
            # Parse results
            hits = self.parse_xml_fast(xml_result, len(sequence))
            
            # Update progress
            with self.progress_lock:
                self.completed += 1
                elapsed = time.time() - self.start_time
                rate = self.completed / elapsed * 3600  # per hour
                remaining = (self.total - self.completed) / (self.completed / elapsed) / 3600 if self.completed > 0 else 0
                
                self.log(f"✅ {protein_id} ({self.completed}/{self.total}) - "
                        f"{rate:.0f}/hr, ETA: {remaining:.1f}h, Hits: {len(hits)}")
            
            return self.create_result(protein, "success", len(hits), hits)
            
        except Exception as e:
            return self.create_result(protein, f"error_{str(e)[:20]}", 0, [])
    
    def create_result(self, protein, status, hit_count, hits):
        """Create standardized result"""
        return {
            'protein_id': protein['protein_id'],
            'protein_name': protein['protein_name'],
            'sequence_length': len(protein['sequence']),
            'status': status,
            'bacterial_hits': hit_count,
            'has_bacterial_homologs': hit_count > 0,
            'best_hit_evalue': hits[0]['evalue'] if hits else None,
            'best_hit_identity': hits[0]['identity_percent'] if hits else None,
            'best_hit_coverage': hits[0]['coverage_percent'] if hits else None,
            'best_hit_description': hits[0]['description'] if hits else None,
            'processing_time': time.strftime('%Y-%m-%d %H:%M:%S')
        }
    
    def run_fast_analysis(self, fasta_file, start_index=0, max_proteins=None):
        """Run fast parallel analysis"""
        
        self.log(f"🚀 Starting fast BLAST analysis")
        self.log(f"Input: {fasta_file}")
        self.log(f"Workers: {self.max_workers}")
        
        # Parse proteins
        proteins = self.parse_proteins(fasta_file)
        
        # Select range
        if max_proteins:
            proteins = proteins[start_index:start_index + max_proteins]
        else:
            proteins = proteins[start_index:]
        
        self.total = len(proteins)
        self.log(f"Processing {self.total} proteins (from index {start_index})")
        
        # Process with ThreadPoolExecutor
        results = []
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit jobs with staggered timing
            futures = []
            
            for i, protein in enumerate(proteins):
                # Stagger submissions
                delay = i * self.submission_delay
                
                def delayed_process(p, d):
                    time.sleep(d)
                    return self.process_protein_worker(p)
                
                future = executor.submit(delayed_process, protein, delay)
                futures.append(future)
            
            # Collect results
            try:
                for future in concurrent.futures.as_completed(futures, timeout=7200):  # 2 hour total timeout
                    try:
                        result = future.result()
                        results.append(result)
                    except Exception as e:
                        self.log(f"❌ Future error: {e}")
                        
            except KeyboardInterrupt:
                self.log("🛑 Interrupted by user")
                executor.shutdown(wait=False)
        
        # Save results
        self.save_results(results)
        return results
    
    def save_results(self, results):
        """Save comprehensive results"""
        
        self.log("💾 Saving results...")
        
        # CSV results
        csv_file = self.results_dir / "fast_blast_results.csv"
        fieldnames = [
            'protein_id', 'protein_name', 'sequence_length', 'status',
            'bacterial_hits', 'has_bacterial_homologs', 
            'best_hit_evalue', 'best_hit_identity', 'best_hit_coverage',
            'best_hit_description', 'processing_time'
        ]
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        
        # Summary statistics
        total = len(results)
        successful = len([r for r in results if r['status'] == 'success'])
        with_hits = len([r for r in results if r['has_bacterial_homologs']])
        
        elapsed = time.time() - self.start_time
        rate = total / elapsed * 3600 if elapsed > 0 else 0
        
        summary = {
            'analysis_date': time.strftime('%Y-%m-%d %H:%M:%S'),
            'pipeline_type': 'fast_blast',
            'max_workers': self.max_workers,
            'total_proteins': total,
            'successful_searches': successful,
            'proteins_with_hits': with_hits,
            'success_rate_percent': round(successful / total * 100, 1) if total > 0 else 0,
            'homolog_rate_percent': round(with_hits / successful * 100, 1) if successful > 0 else 0,
            'total_time_hours': round(elapsed / 3600, 2),
            'processing_rate_per_hour': round(rate, 1),
            'average_time_per_protein_minutes': round(elapsed / total * 60, 1) if total > 0 else 0
        }
        
        json_file = self.results_dir / "fast_blast_summary.json"
        with open(json_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Final stats
        self.log(f"🎉 Analysis complete!")
        self.log(f"   Total: {total}, Success: {successful} ({summary['success_rate_percent']}%)")
        self.log(f"   With hits: {with_hits} ({summary['homolog_rate_percent']}%)")
        self.log(f"   Rate: {summary['processing_rate_per_hour']} proteins/hour")
        self.log(f"   Time: {summary['total_time_hours']} hours")
        self.log(f"   Results: {csv_file}")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Fast BLAST pipeline with threading')
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('--output-dir', default='fast_blast_results', help='Output directory')
    parser.add_argument('--start', type=int, default=0, help='Start protein index')
    parser.add_argument('--max', type=int, help='Max proteins to process')
    parser.add_argument('--workers', type=int, default=6, help='Number of worker threads (1-8)')
    parser.add_argument('--test', action='store_true', help='Test mode - 5 proteins')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.fasta_file):
        print(f"❌ FASTA file not found: {args.fasta_file}")
        sys.exit(1)
    
    args.workers = max(1, min(args.workers, 8))
    
    # Initialize pipeline
    pipeline = FastBlastPipeline(args.output_dir, args.workers)
    
    # Set parameters
    if args.test:
        max_proteins = 5
        start_index = 0
        print("🧪 TEST MODE: 5 proteins with fast pipeline")
    else:
        max_proteins = args.max
        start_index = args.start
    
    # Run analysis
    try:
        results = pipeline.run_fast_analysis(args.fasta_file, start_index, max_proteins)
        
        print(f"\n🎉 Fast BLAST analysis completed!")
        print(f"Results: {pipeline.results_dir}")
        print(f"Logs: {pipeline.log_file}")
        
    except KeyboardInterrupt:
        print("\n🛑 Analysis interrupted")
    except Exception as e:
        print(f"❌ Analysis failed: {e}")

if __name__ == "__main__":
    main()