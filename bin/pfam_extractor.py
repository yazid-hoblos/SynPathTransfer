#!/usr/bin/env python3

import os
import csv
import sys
import subprocess
from pathlib import Path
from typing import Set, List, Dict, Optional

class PfamHMMExtractor:
    def __init__(self, pfam_db_path: str = "Pfam/Pfam-A.hmm"):
        """
        Extract Pfam HMM profiles using hmmfetch, analyze with pyhmmer
        
        Args:
            pfam_db_path: Path to Pfam-A.hmm database file
        """
        self.pfam_db_path = Path(pfam_db_path)
        self.output_dir = Path("pfam_profiles")
        
        # Ensure output directory exists
        try:
            self.output_dir.mkdir(exist_ok=True, parents=True)
            print(f"✓ Output directory ready: {self.output_dir.absolute()}")
        except Exception as e:
            print(f"Error creating output directory: {e}")
            sys.exit(1)
        
        # Check if Pfam database exists
        if not self.pfam_db_path.exists():
            print(f"Error: Pfam database not found at {self.pfam_db_path}")
            print("Please download Pfam-A.hmm from: http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/")
            sys.exit(1)
        
        # Check if HMMER tools are available
        if not self.check_hmmer_tools():
            sys.exit(1)
    
    def check_hmmer_tools(self) -> bool:
        """
        Check if HMMER tools are available
        """
        try:
            subprocess.run(['hmmfetch', '-h'], capture_output=True, check=True)
            print("✓ hmmfetch found")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("✗ hmmfetch not found. Please install HMMER tools:")
            print("  Ubuntu/Debian: sudo apt-get install hmmer")
            print("  MacOS: brew install hmmer")
            print("  Or download from: http://hmmer.org/download.html")
            return False
    
    def extract_pfam_ids_from_csv(self, csv_file: str) -> Set[str]:
        """
        Extract unique Pfam IDs from your CSV output
        
        Args:
            csv_file: Path to your CSV file with format: EC_number,UniProt_ID,Protein_Name,Pfam_ID,Pfam_Description
            
        Returns:
            Set of unique Pfam IDs
        """
        pfam_ids = set()
        
        try:
            with open(csv_file, 'r') as f:
                reader = csv.reader(f)
                
                # Skip header if present
                first_row = next(reader)
                if not first_row[3].startswith('PF'):  # Assuming Pfam IDs start with PF
                    # First row is header, continue with next rows
                    pass
                else:
                    # First row contains data, process it
                    pfam_ids.add(first_row[3].strip())
                
                # Process remaining rows
                for row in reader:
                    if len(row) >= 4 and row[3].strip():  # Make sure Pfam_ID field exists and is not empty
                        pfam_ids.add(row[3].strip())
            
            print(f"Found {len(pfam_ids)} unique Pfam IDs: {sorted(pfam_ids)}")
            return pfam_ids
            
        except Exception as e:
            print(f"Error reading CSV file {csv_file}: {e}")
            return set()
    
    def create_pfam_list_file(self, pfam_ids: Set[str], filename: str = "pfam_ids.txt") -> str:
        """
        Create a text file with Pfam IDs for hmmfetch
        
        Args:
            pfam_ids: Set of Pfam IDs
            filename: Output filename
            
        Returns:
            Path to created file
        """
        list_file = self.output_dir / filename
        
        with open(list_file, 'w') as f:
            for pfam_id in sorted(pfam_ids):
                f.write(f"{pfam_id}\n")
        
        print(f"Created Pfam ID list: {list_file}")
        return str(list_file)
    
    def find_versioned_pfam_ids(self, pfam_ids: Set[str]) -> Dict[str, str]:
        # find each base Pfam ID in Pfam-A.hmm to get its specific version
        
        list_file = self.create_pfam_list_file(pfam_ids, "base_pfam_ids.txt")
        # grep -f list_file Pfam-A.hmm | awk '/^ACC/ {print $2}'
        
        subprocess_args = [
            'grep', '-f', list_file, str(self.pfam_db_path)
        ]
        versioned_ids = set()
        try:
            grep_process = subprocess.Popen(subprocess_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            awk_process = subprocess.Popen(['awk', '/^ACC/ {print $2}'], stdin=grep_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            grep_process.stdout.close()  # Allow grep to receive a SIGPIPE if awk exits
            
            out, err = awk_process.communicate()
            if awk_process.returncode == 0:
                for line in out.strip().split('\n'):
                    if line:
                        versioned_ids.add(line.strip())
            else:
                print(f"Error in awk processing: {err}")
        except Exception as e:
            print(f"Error finding versioned Pfam IDs: {e}")
        
        print(f"Resolved {len(versioned_ids)} versioned Pfam IDs")
        print(f"Versioned IDs: {sorted(versioned_ids)}")
        
        return versioned_ids
        
    def extract_hmm_profiles_hmmfetch(self, pfam_ids: Set[str], output_file: str = "pathway_profiles.hmm") -> bool:
        """
        Extract HMM profiles using hmmfetch (reliable approach)
        
        Args:
            pfam_ids: Set of Pfam IDs to extract
            output_file: Output HMM file name
            
        Returns:
            True if successful, False otherwise
        """
        output_path = self.output_dir / output_file
        
        # Find versioned IDs first
        versioned_ids = self.find_versioned_pfam_ids(pfam_ids)
        
        # Create temporary file with versioned Pfam IDs
        pfam_list_file = self.create_pfam_list_file(versioned_ids, "versioned_pfam_ids.txt")
                
        try:
            print(f"Extracting {len(versioned_ids)} HMM profiles using hmmfetch...")
            print(f"Command: hmmfetch -f {self.pfam_db_path} {pfam_list_file}")
            
            # Run hmmfetch
            with open(output_path, 'w') as outf:
                result = subprocess.run([
                    'hmmfetch', 
                    '-f',  # Fetch multiple profiles from file
                    str(self.pfam_db_path), 
                    pfam_list_file
                ], stdout=outf, stderr=subprocess.PIPE, text=True)
            
            if result.returncode == 0:
                print(f"✓ Successfully extracted HMM profiles to: {output_path}")
                
                # Check output file size
                file_size = output_path.stat().st_size
                if file_size > 0:
                    print(f"  Output file size: {file_size:,} bytes")
                    
                    # Count how many profiles were extracted
                    profile_count = self.count_hmm_profiles(output_path)
                    print(f"  Profiles extracted: {profile_count}")
                    
                    return True
                else:
                    print("✗ Output file is empty - no profiles found")
                    return False
            else:
                print(f"✗ hmmfetch failed with error: {result.stderr}")
                return False
                
        except Exception as e:
            print(f"Error running hmmfetch: {e}")
            return False
    
    def extract_individual_profiles_hmmfetch(self, pfam_ids: Set[str]) -> Dict[str, bool]:
        """
        Extract each Pfam profile as a separate file using hmmfetch with version resolution
        
        Args:
            pfam_ids: Set of base Pfam IDs to extract
            
        Returns:
            Dictionary mapping base Pfam ID to success status
        """
        results = {}
        individual_dir = self.output_dir / "individual_profiles"
        individual_dir.mkdir(exist_ok=True)
        
        # Find versioned IDs first
        pfam_mapping = self.find_versioned_pfam_ids(pfam_ids)
        
        for base_id, versioned_id in pfam_mapping.items():
            output_file = individual_dir / f"{base_id}.hmm"  # Use base ID for filename
            
            try:
                with open(output_file, 'w') as outf:
                    result = subprocess.run([
                        'hmmfetch', 
                        str(self.pfam_db_path), 
                        versioned_id  # Use versioned ID for hmmfetch
                    ], stdout=outf, stderr=subprocess.PIPE, text=True)
                
                if result.returncode == 0 and output_file.stat().st_size > 0:
                    print(f"✓ {base_id} ({versioned_id}) -> {output_file}")
                    results[base_id] = True
                else:
                    print(f"✗ {base_id} ({versioned_id}) - not found or empty")
                    results[base_id] = False
                    
            except Exception as e:
                print(f"✗ {base_id} - error: {e}")
                results[base_id] = False
        
        # Handle IDs that weren't found in version mapping
        missing_ids = pfam_ids - set(pfam_mapping.keys())
        for missing_id in missing_ids:
            results[missing_id] = False
            print(f"✗ {missing_id} - no version found in database")
        
        success_count = sum(results.values())
        print(f"\nExtracted {success_count}/{len(pfam_ids)} profiles successfully")
        
        return results
    
    def count_hmm_profiles(self, hmm_file: Path) -> int:
        """
        Count number of HMM profiles in a file
        """
        try:
            with open(hmm_file, 'r') as f:
                count = 0
                for line in f:
                    if line.startswith('HMMER3/f'):
                        count += 1
            return count
        except Exception:
            return 0
    
    def search_target_sequences_pyhmmer(self, hmm_file: str, target_fasta: str, 
                                       e_value: float = 1e-5) -> Optional[List]:
        """
        Search target sequences using pyhmmer (for analysis only)
        
        Args:
            hmm_file: Path to HMM file created by hmmfetch
            target_fasta: Path to target genome FASTA file
            e_value: E-value threshold for hits
            
        Returns:
            List of hits or None if error
        """
        try:
            import pyhmmer
            print(f"Searching target sequences with pyhmmer...")
            print(f"HMM file: {hmm_file}")
            print(f"Target file: {target_fasta}")
            print(f"E-value threshold: {e_value}")
            
            # Load HMMs from hmmfetch output
            with pyhmmer.plan7.HMMFile(hmm_file) as hmm_file_obj:
                hmms = list(hmm_file_obj)
            print(f"✓ Loaded {len(hmms)} HMM profiles")
            
            # Load target sequences
            with pyhmmer.easel.SequenceFile(target_fasta, digital=True) as seq_file:
                target_sequences = list(seq_file)
            print(f"✓ Loaded {len(target_sequences)} target sequences")
            
            # Run hmmsearch
            print("Running hmmsearch with pyhmmer...")
            all_hits = []
            for hits in pyhmmer.hmmsearch(hmms, target_sequences, E=e_value):
                all_hits.extend(hits)
            
            print(f"✓ Search complete! Found {len(all_hits)} hits above E-value threshold")
            return all_hits
            
        except ImportError:
            print("pyhmmer not available for search. Install with: pip install pyhmmer")
            return None
        except Exception as e:
            print(f"Error during pyhmmer search: {e}")
            return None
    
    def analyze_pathway_coverage(self, pfam_ids: Set[str], target_fasta: str, 
                                e_value: float = 1e-5) -> Dict[str, List]:
        """
        Analyze which pathway components are present in target genome
        
        Args:
            pfam_ids: Set of Pfam IDs representing pathway components
            target_fasta: Path to target genome FASTA file
            e_value: E-value threshold for hits
            
        Returns:
            Dictionary mapping Pfam ID to list of hits
        """
        print("\n=== Pathway Coverage Analysis ===")
        
        # Use the extracted HMM file
        hmm_file = self.output_dir / "pathway_profiles.hmm"
        
        if not hmm_file.exists():
            print(f"Error: HMM file not found: {hmm_file}")
            return {}
        
        
        # Run search using pyhmmer
        all_hits = self.search_target_sequences_pyhmmer(str(hmm_file), target_fasta, e_value)
        
        if not all_hits:
            return {}
        
        # Organize hits by Pfam ID
        coverage = {pfam_id: [] for pfam_id in pfam_ids}
        
        for hit in all_hits:
            target_name = hit.name.decode()       # sequence name
            target_acc  = hit.accession.decode() if hit.accession else None
            
            # hmm_name comes from the query, not the Hit
            hmm_name = hmm.name.decode()
            
            # Try to match HMM name to Pfam ID
            matching_pfam = None
            for pfam_id in pfam_ids:
                if pfam_id in hmm_name or hmm_name == pfam_id:
                    matching_pfam = pfam_id
                    break
            
            if matching_pfam:
                coverage[matching_pfam].append({
                    'target': target_name,
                    'target_acc': target_acc,
                    'e_value': hit.evalue,
                    'score': hit.score,
                    'bias': hit.bias,
                    'hmm_name': hmm_name
                })

        # Generate report
        print(f"\n{'Pfam ID':<12} {'Status':<15} {'Hits':<6} {'Best E-value':<15}")
        print("-" * 60)
        
        present_count = 0
        for pfam_id in sorted(pfam_ids):
            hits = coverage[pfam_id]
            
            if hits:
                best_eval = min(hit['e_value'] for hit in hits)
                status = "PRESENT"
                present_count += 1
                print(f"{pfam_id:<12} {status:<15} {len(hits):<6} {best_eval:<15.2e}")
            else:
                status = "MISSING"
                print(f"{pfam_id:<12} {status:<15} {0:<6} {'N/A':<15}")
        
        print(f"\nSummary: {present_count}/{len(pfam_ids)} pathway components found in target genome")
        
        return coverage
    
    def save_results(self, coverage: Dict[str, List], output_file: str = "pathway_coverage_report.csv"):
        """
        Save pathway coverage analysis results
        
        Args:
            coverage: Results from analyze_pathway_coverage
            output_file: Output CSV file name
        """
        output_path = self.output_dir / output_file
        
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Pfam_ID', 'Status', 'Hit_Count', 'Best_E_value', 'Best_Score', 'Target_Proteins', 'HMM_Names'])
            
            for pfam_id in sorted(coverage.keys()):
                hits = coverage[pfam_id]
                
                if hits:
                    best_hit = min(hits, key=lambda h: h['e_value'])
                    target_proteins = ';'.join([hit['target'] for hit in hits])
                    hmm_names = ';'.join(list(set([hit['hmm_name'] for hit in hits])))
                    
                    writer.writerow([
                        pfam_id,
                        'PRESENT',
                        len(hits),
                        f"{best_hit['e_value']:.2e}",
                        f"{best_hit['score']:.2f}",
                        target_proteins,
                        hmm_names
                    ])
                else:
                    writer.writerow([pfam_id, 'MISSING', 0, 'N/A', 'N/A', '', ''])
        
        print(f"✓ Detailed results saved to: {output_path}")
    
    def run_complete_analysis(self, csv_file: str, target_fasta: str = None, 
                             extract_mode: str = "combined", e_value: float = 1e-5) -> bool:
        """
        Complete automated analysis workflow: hmmfetch + pyhmmer
        
        Args:
            csv_file: Path to your CSV results file
            target_fasta: Path to target genome FASTA (optional)
            extract_mode: 'combined' or 'individual'
            e_value: E-value threshold for searches
            
        Returns:
            True if successful
        """
        print("=== Pfam HMM Analysis: hmmfetch + pyhmmer ===")
        print(f"Input CSV: {csv_file}")
        print(f"Pfam DB: {self.pfam_db_path}")
        print(f"Output directory: {self.output_dir}")
        print(f"Extraction mode: {extract_mode}")
        if target_fasta:
            print(f"Target genome: {target_fasta}")
        print()
        
        # Step 1: Extract Pfam IDs from CSV
        pfam_ids = self.extract_pfam_ids_from_csv(csv_file)
        if not pfam_ids:
            print("No Pfam IDs found in CSV file")
            return False
        
        # Step 2: Extract HMM profiles using hmmfetch
        if extract_mode == "combined":
            success = self.extract_hmm_profiles_hmmfetch(pfam_ids)
        elif extract_mode == "individual":
            results = self.extract_individual_profiles_hmmfetch(pfam_ids)
            success = any(results.values())
        else:
            print(f"Unknown extraction mode: {extract_mode}")
            return False
        
        if not success:
            print("No HMM profiles were successfully extracted")
            return False
        
        # Step 3: Optional target genome analysis using pyhmmer
        if target_fasta and Path(target_fasta).exists():
            coverage = self.analyze_pathway_coverage(pfam_ids, target_fasta, e_value)
            
            if coverage:
                self.save_results(coverage)
        
        print(f"\n✓ Analysis complete! Files saved to: {self.output_dir}")
        return True
    
    def create_hmmsearch_script(self, target_genome_file: str = "target_genome.fasta"):
        """
        Create a script to run hmmsearch on target genome (alternative to pyhmmer)
        
        Args:
            target_genome_file: Path to target genome FASTA file
        """
        script_content = f"""#!/bin/bash

# Automated hmmsearch script for pathway analysis
# Generated by PfamHMMExtractor

PFAM_PROFILES="{self.output_dir}/pathway_profiles.hmm"
TARGET_GENOME="{target_genome_file}"
OUTPUT_DIR="{self.output_dir}/hmmsearch_results"

echo "Running hmmsearch analysis..."
echo "Profiles: $PFAM_PROFILES"
echo "Target: $TARGET_GENOME"
echo "Output: $OUTPUT_DIR"

mkdir -p $OUTPUT_DIR

# Run hmmsearch
hmmsearch --tblout $OUTPUT_DIR/hits.tbl \\
          --domtblout $OUTPUT_DIR/domain_hits.tbl \\
          -E 1e-5 \\
          $PFAM_PROFILES \\
          $TARGET_GENOME

echo "Analysis complete! Results in $OUTPUT_DIR/"
echo "Key files:"
echo "  - hits.tbl: Summary of all hits"
echo "  - domain_hits.tbl: Detailed domain matches"
"""
        
        script_file = self.output_dir / "run_hmmsearch.sh"
        with open(script_file, 'w') as f:
            f.write(script_content)
        
        # Make executable
        os.chmod(script_file, 0o755)
        
        print(f"Created hmmsearch script: {script_file}")
        print(f"Usage: ./{script_file}")


def main():
    """
    Command line interface
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Extract Pfam HMM profiles with hmmfetch, analyze with pyhmmer")
    parser.add_argument("csv_file", help="CSV file with Pfam IDs (from your UniProt script)")
    parser.add_argument("--pfam-db", default="Pfam/Pfam-A.hmm", 
                       help="Path to Pfam-A.hmm database (default: Pfam/Pfam-A.hmm)")
    parser.add_argument("--target", help="Target genome FASTA file for coverage analysis")
    parser.add_argument("--mode", choices=['combined', 'individual'], default='combined',
                       help="Extract as one combined file or individual files")
    parser.add_argument("--evalue", type=float, default=1e-5,
                       help="E-value threshold for hmmsearch (default: 1e-5)")
    parser.add_argument("--hmmsearch", action='store_true',
                       help="Run hmmsearch using pyhmmer on extracted profiles")
    parser.add_argument("--hmmsearch-target", help="Target genome for hmmsearch (if different from --target)")
    
    args = parser.parse_args()
    
    # Initialize extractor
    extractor = PfamHMMExtractor(args.pfam_db)
    
    # Run analysis
    success = extractor.run_complete_analysis(
        args.csv_file, 
        args.target, 
        args.mode, 
        args.evalue
    )
    
    # Optionally run hmmsearch separately
    if args.hmmsearch and success:
        target = args.hmmsearch_target or args.target or "target_genome.fasta"
        extractor.run_hmmsearch_pyhmmer(target_genome_file=target, e_value=args.evalue)
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())