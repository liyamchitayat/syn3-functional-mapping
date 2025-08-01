#!/usr/bin/env python3
"""
Main entry point for the Syn3A BLAST pipeline
"""

import sys
import argparse
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from pipeline.fast_blast_pipeline import FastBlastPipeline

def main():
    parser = argparse.ArgumentParser(
        description='Syn3A BLAST Pipeline - Find bacterial homologs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Quick test (5 proteins)
  python run_pipeline.py data/syn3A_proteins.fasta --test
  
  # Full analysis with 8 workers
  python run_pipeline.py data/syn3A_proteins.fasta --workers 8
  
  # Resume from protein 100
  python run_pipeline.py data/syn3A_proteins.fasta --start 100
        """
    )
    
    parser.add_argument('fasta_file', help='Input FASTA file with syn3A proteins')
    parser.add_argument('--output-dir', default='results', help='Output directory')
    parser.add_argument('--workers', type=int, default=6, help='Number of parallel workers (1-50)')
    parser.add_argument('--start', type=int, default=0, help='Start from protein index')
    parser.add_argument('--max', type=int, help='Maximum proteins to process')
    parser.add_argument('--test', action='store_true', help='Test mode - process 5 proteins')
    
    args = parser.parse_args()
    
    # Validate input
    if not Path(args.fasta_file).exists():
        print(f"❌ Error: FASTA file not found: {args.fasta_file}")
        print("💡 Tip: Download syn3A proteins first:")
        print("   python src/utils/download_syn3a.py")
        sys.exit(1)
    
    # Initialize pipeline
    pipeline = FastBlastPipeline(
        output_dir=args.output_dir,
        max_workers=min(args.workers, 50)
    )
    
    # Run analysis
    if args.test:
        print("🧪 Running in TEST MODE (5 proteins)")
        results = pipeline.run_fast_analysis(args.fasta_file, 0, 5)
    else:
        results = pipeline.run_fast_analysis(
            args.fasta_file, 
            args.start, 
            args.max
        )
    
    print(f"\n✅ Pipeline completed!")
    print(f"📊 Results saved to: {args.output_dir}/")

if __name__ == "__main__":
    main()