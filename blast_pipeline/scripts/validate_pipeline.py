#!/usr/bin/env python3
"""
Comprehensive Pipeline Validation Script
Tests all components before running the full analysis
"""

import os
import sys
import json
import time
import requests
from pathlib import Path
import re

def test_environment():
    """Test basic environment setup"""
    print("🔧 Testing Environment Setup...")
    
    tests = []
    
    # Python version
    if sys.version_info >= (3, 6):
        tests.append("✅ Python version OK (3.6+)")
    else:
        tests.append("❌ Python version too old")
    
    # Required modules
    try:
        import requests
        tests.append("✅ requests module available")
    except ImportError:
        tests.append("❌ requests module missing")
    
    try:
        import json, csv, pathlib, re, urllib.parse
        tests.append("✅ Standard library modules OK")
    except ImportError as e:
        tests.append(f"❌ Missing standard library: {e}")
    
    # File system permissions
    try:
        test_dir = Path("test_permissions")
        test_dir.mkdir(exist_ok=True)
        test_file = test_dir / "test.txt"
        test_file.write_text("test")
        test_file.unlink()
        test_dir.rmdir()
        tests.append("✅ File system permissions OK")
    except Exception as e:
        tests.append(f"❌ File system error: {e}")
    
    for test in tests:
        print(f"   {test}")
    
    return all("✅" in test for test in tests)

def test_input_files():
    """Test input file availability and format"""
    print("\n📄 Testing Input Files...")
    
    tests = []
    
    # Check if syn3A_proteins.fasta exists
    fasta_file = Path("../syn3A_proteins.fasta")
    if fasta_file.exists():
        tests.append("✅ syn3A_proteins.fasta found")
        
        # Check file size
        size = fasta_file.stat().st_size
        if size > 100000:  # Should be ~235KB
            tests.append(f"✅ File size OK ({size:,} bytes)")
        else:
            tests.append(f"❌ File too small ({size} bytes)")
        
        # Check FASTA format
        try:
            with open(fasta_file, 'r') as f:
                content = f.read(1000)  # Read first 1KB
                
            if content.startswith('>'):
                tests.append("✅ FASTA format OK")
            else:
                tests.append("❌ Invalid FASTA format")
            
            # Count proteins
            with open(fasta_file, 'r') as f:
                protein_count = sum(1 for line in f if line.startswith('>'))
            
            if protein_count == 438:
                tests.append(f"✅ Protein count correct ({protein_count})")
            else:
                tests.append(f"❌ Wrong protein count ({protein_count}, expected 438)")
            
        except Exception as e:
            tests.append(f"❌ Error reading FASTA: {e}")
    else:
        tests.append("❌ syn3A_proteins.fasta not found")
    
    for test in tests:
        print(f"   {test}")
    
    return all("✅" in test for test in tests)

def test_network_connectivity():
    """Test network connectivity to NCBI"""
    print("\n🌐 Testing Network Connectivity...")
    
    tests = []
    
    # Test basic internet
    try:
        response = requests.get("https://www.google.com", timeout=5)
        if response.status_code == 200:
            tests.append("✅ Internet connectivity OK")
        else:
            tests.append(f"❌ Internet issue (status {response.status_code})")
    except Exception as e:
        tests.append(f"❌ Internet connectivity failed: {e}")
    
    # Test NCBI BLAST service
    try:
        response = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi", timeout=10)
        if response.status_code == 200:
            tests.append("✅ NCBI BLAST service accessible")
        else:
            tests.append(f"❌ NCBI BLAST issue (status {response.status_code})")
    except Exception as e:
        tests.append(f"❌ NCBI BLAST not accessible: {e}")
    
    # Test NCBI E-utilities
    try:
        response = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi", timeout=10)
        if response.status_code == 200:
            tests.append("✅ NCBI E-utilities accessible")
        else:
            tests.append(f"❌ NCBI E-utilities issue")
    except Exception as e:
        tests.append(f"❌ NCBI E-utilities not accessible: {e}")
    
    for test in tests:
        print(f"   {test}")
    
    return all("✅" in test for test in tests)

def test_pipeline_components():
    """Test pipeline script components"""
    print("\n🔧 Testing Pipeline Components...")
    
    tests = []
    
    # Check script exists
    script_path = Path("scripts/web_blast_pipeline.py")
    if script_path.exists():
        tests.append("✅ web_blast_pipeline.py found")
        
        # Check if executable
        if os.access(script_path, os.X_OK):
            tests.append("✅ Script is executable")
        else:
            tests.append("❌ Script not executable")
        
        # Test help command
        try:
            import subprocess
            result = subprocess.run([sys.executable, str(script_path), "--help"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0 and "usage:" in result.stdout:
                tests.append("✅ Script help command works")
            else:
                tests.append(f"❌ Script help failed: {result.stderr}")
        except Exception as e:
            tests.append(f"❌ Script test failed: {e}")
    else:
        tests.append("❌ web_blast_pipeline.py not found")
    
    for test in tests:
        print(f"   {test}")
    
    return all("✅" in test for test in tests)

def test_blast_submission():
    """Test BLAST submission (without waiting for results)"""
    print("\n🧬 Testing BLAST Submission...")
    
    tests = []
    
    try:
        # Simple test sequence (short)
        test_sequence = "MKKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDCLTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGCDITIILS"
        
        submit_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
        
        params = {
            'CMD': 'Put',
            'PROGRAM': 'blastp',
            'DATABASE': 'nr',
            'QUERY': test_sequence,
            'FORMAT_TYPE': 'XML',
            'HITLIST_SIZE': 10,
            'EXPECT': 0.01,
            'ENTREZ_QUERY': 'NOT mycoplasma[Organism]'
        }
        
        from urllib.parse import urlencode
        data = urlencode(params).encode('utf-8')
        
        print("   📤 Submitting test BLAST...")
        response = requests.post(submit_url, data=data, timeout=30)
        
        if response.status_code == 200:
            tests.append("✅ BLAST submission successful")
            
            # Check for RID
            rid_match = re.search(r'RID = (\w+)', response.text)
            if rid_match:
                rid = rid_match.group(1)
                tests.append(f"✅ BLAST RID obtained: {rid}")
                
                # Test status check (don't wait for completion)
                check_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID={rid}"
                status_response = requests.get(check_url, timeout=10)
                
                if status_response.status_code == 200:
                    if "Status=WAITING" in status_response.text:
                        tests.append("✅ BLAST status check works (job queued)")
                    elif "Status=READY" in status_response.text:
                        tests.append("✅ BLAST status check works (job ready)")
                    else:
                        tests.append("✅ BLAST status check works")
                else:
                    tests.append("❌ BLAST status check failed")
            else:
                tests.append("❌ No BLAST RID in response")
        else:
            tests.append(f"❌ BLAST submission failed (status {response.status_code})")
            
    except Exception as e:
        tests.append(f"❌ BLAST submission error: {e}")
    
    for test in tests:
        print(f"   {test}")
    
    return all("✅" in test for test in tests)

def test_pipeline_parsing():
    """Test pipeline parsing functions"""
    print("\n📊 Testing Pipeline Parsing...")
    
    tests = []
    
    try:
        # Test FASTA parsing
        test_fasta = """>test_protein
MKKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDCLTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGCDITIILS
>protein2
TESTSEQUENCE"""
        
        with open("test_parse.fasta", "w") as f:
            f.write(test_fasta)
        
        # Import and test parsing function
        sys.path.append("scripts")
        from web_blast_pipeline import WebBlastPipeline
        
        pipeline = WebBlastPipeline("test_output")
        proteins = pipeline.parse_fasta_file("test_parse.fasta")
        
        if len(proteins) == 2:
            tests.append("✅ FASTA parsing works")
        else:
            tests.append(f"❌ FASTA parsing error (got {len(proteins)} proteins)")
        
        if proteins[0]['sequence'] == "MKKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDCLTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGCDITIILS":
            tests.append("✅ Sequence extraction correct")
        else:
            tests.append("❌ Sequence extraction error")
        
        # Clean up
        os.unlink("test_parse.fasta")
        
    except Exception as e:
        tests.append(f"❌ Parsing test error: {e}")
    
    for test in tests:
        print(f"   {test}")
    
    return all("✅" in test for test in tests)

def main():
    """Run all validation tests"""
    print("=" * 60)
    print("🧪 BLAST PIPELINE VALIDATION TESTS")
    print("=" * 60)
    
    results = []
    
    # Run all tests
    results.append(("Environment", test_environment()))
    results.append(("Input Files", test_input_files()))
    results.append(("Network Connectivity", test_network_connectivity()))
    results.append(("Pipeline Components", test_pipeline_components()))
    results.append(("BLAST Submission", test_blast_submission()))
    results.append(("Pipeline Parsing", test_pipeline_parsing()))
    
    # Summary
    print("\n" + "=" * 60)
    print("📋 VALIDATION SUMMARY")
    print("=" * 60)
    
    all_passed = True
    for test_name, passed in results:
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{test_name:20} {status}")
        if not passed:
            all_passed = False
    
    print("\n" + "=" * 60)
    if all_passed:
        print("🎉 ALL TESTS PASSED! Pipeline is ready to run.")
        print("\n📋 RECOMMENDED COMMANDS:")
        print("   # Test with 3 proteins:")
        print("   python3 scripts/web_blast_pipeline.py ../syn3A_proteins.fasta --test")
        print("\n   # Full analysis (background):")
        print("   nohup python3 scripts/web_blast_pipeline.py ../syn3A_proteins.fasta --output-dir full_results > blast_full.log 2>&1 &")
        print("\n   # Monitor progress:")
        print("   tail -f blast_full.log")
    else:
        print("❌ SOME TESTS FAILED! Fix issues before running full pipeline.")
    
    print("=" * 60)
    
    return all_passed

if __name__ == "__main__":
    main()