import re
import time
from typing import Set, List, Dict, Optional, Tuple
from Bio.KEGG import REST


class KEGGModuleDiscovery:
    def __init__(self, delay: float = 0.1):
        """
        KEGG pathway module discovery and analysis tool
        
        Args:
            delay: seconds to wait between API calls
        """
        self.delay = delay
        self.module_cache = {}  # Cache module data to avoid repeated API calls
        
    def _standardize_pathway_id(self, pathway_id: str) -> str:
        """Convert pathway ID to standard KEGG format"""
        if pathway_id.startswith('ko') or pathway_id.startswith('map'):
            return pathway_id
        else:
            return f"ko{pathway_id}"
    
    def discover_pathway_modules(self, pathway_id: str) -> List[Dict[str, str]]:
        """
        Automatically discover all modules within a pathway
        
        Args:
            pathway_id: KEGG pathway ID (e.g., 'ko00720', 'map00720', '00720')
            
        Returns:
            List of dictionaries with module info: [{'id': 'M00376', 'name': '...', 'description': '...'}, ...]
        """
        pathway_id = self._standardize_pathway_id(pathway_id)
        
        try:
            print(f"Discovering modules in pathway {pathway_id}...")
            pathway_data = REST.kegg_get(pathway_id).read()
            
            # Extract all module references from pathway
            # Pattern matches MD:M##### or just M##### in KEGG data
            module_pattern = r'(?:MD:)?(M\d{5})'
            module_ids = list(set(re.findall(module_pattern, pathway_data)))
            
            print(f"Found {len(module_ids)} modules: {module_ids}")
            
            # Get detailed information for each module
            modules_info = []
            for module_id in module_ids:
                module_info = self._get_module_info(module_id)
                if module_info:
                    modules_info.append(module_info)
                time.sleep(self.delay)
            
            return modules_info
            
        except Exception as e:
            print(f"Error discovering modules in {pathway_id}: {e}")
            return []
    
    def _get_module_info(self, module_id: str) -> Optional[Dict[str, str]]:
        """
        Get detailed information about a specific module
        
        Args:
            module_id: KEGG module ID (e.g., 'M00376')
            
        Returns:
            Dictionary with module information
        """
        if module_id in self.module_cache:
            return self.module_cache[module_id]
        
        try:
            module_data = REST.kegg_get(module_id).read()
            
            info = {
                'id': module_id,
                'name': '',
                'definition': '',
                'class': '',
                'pathway': '',
                'reaction_count': 0,
                'orthology_count': 0
            }
            
            lines = module_data.split('\n')
            current_section = None
            
            for line in lines:
                line = line.strip()
                
                if line.startswith('NAME'):
                    info['name'] = line.replace('NAME', '').strip()
                elif line.startswith('DEFINITION'):
                    info['definition'] = line.replace('DEFINITION', '').strip()
                elif line.startswith('CLASS'):
                    info['class'] = line.replace('CLASS', '').strip()
                elif line.startswith('PATHWAY'):
                    current_section = 'pathway'
                elif current_section == 'pathway' and line and not line.startswith(' '):
                    current_section = None
                elif current_section == 'pathway' and line.startswith(' '):
                    info['pathway'] += line.strip() + '; '
                elif line.startswith('REACTION'):
                    # Count reactions in module
                    reaction_count = len(re.findall(r'R\d{5}', module_data))
                    info['reaction_count'] = reaction_count
                elif line.startswith('ORTHOLOGY'):
                    # Count KO numbers in module  
                    ko_count = len(re.findall(r'K\d{5}', module_data))
                    info['orthology_count'] = ko_count
            
            # Clean up pathway field
            info['pathway'] = info['pathway'].rstrip('; ')
            
            # Cache the result
            self.module_cache[module_id] = info
            return info
            
        except Exception as e:
            print(f"Error getting info for module {module_id}: {e}")
            return None
    
    def get_module_ec_numbers(self, module_id: str) -> Set[str]:
        """
        Extract EC numbers from a specific KEGG module
        
        Args:
            module_id: KEGG module ID (e.g., 'M00376')
            
        Returns:
            Set of EC numbers found in the module
        """
        try:
            print(f"Extracting EC numbers from module {module_id}...")
            module_data = REST.kegg_get(module_id).read()
            
            # Extract EC numbers from module
            ec_pattern = r'(\d+\.\d+\.\d+\.\d+)'
            ec_numbers = set(re.findall(ec_pattern, module_data))
            
            print(f"  Found {len(ec_numbers)} EC numbers in {module_id}")
            return ec_numbers
            
        except Exception as e:
            print(f"Error extracting EC numbers from {module_id}: {e}")
            return set()
    
    def get_pathway_ec_numbers(self, pathway_id: str) -> Set[str]:
        """
        Extract EC numbers directly from pathway (without module breakdown)
        
        Args:
            pathway_id: KEGG pathway ID
            
        Returns:
            Set of EC numbers found in the pathway
        """
        pathway_id = self._standardize_pathway_id(pathway_id)
        
        try:
            print(f"Extracting EC numbers from pathway {pathway_id}...")
            pathway_data = REST.kegg_get(pathway_id).read()
            
            # Method 1: Direct EC extraction
            ec_pattern_direct = r'EC:(\d+\.\d+\.\d+\.\d+)'
            ec_direct = set(re.findall(ec_pattern_direct, pathway_data))
            
            # Method 2: Via KO numbers for more comprehensive results
            ko_pattern = r'K\d{5}'
            ko_numbers = set(re.findall(ko_pattern, pathway_data))
            
            ec_via_ko = set()
            print(f"  Found {len(ko_numbers)} KO numbers, extracting their EC numbers...")
            
            for i, ko in enumerate(ko_numbers):
                try:
                    ko_data = REST.kegg_get(ko).read()
                    ko_ecs = re.findall(r'EC:(\d+\.\d+\.\d+\.\d+)', ko_data)
                    ec_via_ko.update(ko_ecs)
                    
                    if i % 20 == 0:  # Progress update
                        print(f"    Processed {i+1}/{len(ko_numbers)} KO numbers...")
                    
                    time.sleep(self.delay)
                    
                except Exception as e:
                    print(f"    Error with KO {ko}: {e}")
                    continue
            
            # Combine both methods
            all_ecs = ec_direct.union(ec_via_ko)
            print(f"  Total EC numbers: {len(all_ecs)} (direct: {len(ec_direct)}, via KO: {len(ec_via_ko)})")
            
            return all_ecs
            
        except Exception as e:
            print(f"Error extracting EC numbers from {pathway_id}: {e}")
            return set()
    
    def extract_pathway_components(self, pathway_id: str, granularity: str = 'pathway') -> Dict:
        """
        Extract pathway components with different levels of granularity
        
        Args:
            pathway_id: KEGG pathway ID
            granularity: 'pathway', 'modules', or 'reactions'
            
        Returns:
            Dictionary with extracted components based on granularity level
        """
        pathway_id = self._standardize_pathway_id(pathway_id)
        
        if granularity == 'pathway':
            # Return all EC numbers from the entire pathway
            ec_numbers = self.get_pathway_ec_numbers(pathway_id)
            return {
                'pathway_id': pathway_id,
                'granularity': 'pathway',
                'ec_numbers': ec_numbers,
                'total_ec_count': len(ec_numbers)
            }
            
        elif granularity == 'modules':
            # Break down by individual modules
            modules = self.discover_pathway_modules(pathway_id)
            
            result = {
                'pathway_id': pathway_id,
                'granularity': 'modules',
                'modules': {},
                'total_modules': len(modules),
                'total_ec_count': 0
            }
            
            for module_info in modules:
                module_id = module_info['id']
                ec_numbers = self.get_module_ec_numbers(module_id)
                
                result['modules'][module_id] = {
                    'info': module_info,
                    'ec_numbers': ec_numbers,
                    'ec_count': len(ec_numbers)
                }
                
                result['total_ec_count'] += len(ec_numbers)
                time.sleep(self.delay)
            
            return result
            
        elif granularity == 'reactions':
            # Break down by individual reactions (implementation placeholder)
            print("Reaction-level granularity not yet implemented")
            return self.extract_pathway_components(pathway_id, 'pathway')
            
        else:
            raise ValueError("Granularity must be 'pathway', 'modules', or 'reactions'")
    
    def interactive_module_selection(self, pathway_id: str) -> Dict:
        """
        Interactive workflow for module selection
        Implements the user workflow you described
        """
        print(f"\n{'='*60}")
        print(f"PATHWAY ANALYSIS: {pathway_id}")
        print(f"{'='*60}")
        
        # Step 1: Discover modules
        modules = self.discover_pathway_modules(pathway_id)
        
        if not modules:
            print("No modules found in this pathway.")
            print("Extracting EC numbers from entire pathway...")
            return self.extract_pathway_components(pathway_id, 'pathway')
        
        # Step 2: Display options to user
        print(f"\nFound {len(modules)} modules in this pathway:")
        print(f"{'#':<3} {'Module ID':<10} {'Name':<50}")
        print("-" * 70)
        
        for i, module in enumerate(modules):
            name = module['name'][:47] + "..." if len(module['name']) > 50 else module['name']
            print(f"{i+1:<3} {module['id']:<10} {name:<50}")
        
        # Step 3: Get user choice
        while True:
            try:
                print(f"\nOptions:")
                print(f"  1-{len(modules)}: Select specific module")
                print(f"  'all': Extract from all modules separately")
                print(f"  'entire': Extract from entire pathway (ignore modules)")
                
                choice = input("Your choice: ").strip().lower()
                
                if choice == 'all':
                    return self.extract_pathway_components(pathway_id, 'modules')
                elif choice == 'entire':
                    return self.extract_pathway_components(pathway_id, 'pathway')
                elif choice.isdigit():
                    idx = int(choice) - 1
                    if 0 <= idx < len(modules):
                        selected_module = modules[idx]
                        module_id = selected_module['id']
                        ec_numbers = self.get_module_ec_numbers(module_id)
                        return {
                            'pathway_id': pathway_id,
                            'granularity': 'single_module',  # <- ADD THIS LINE
                            'selected_module': selected_module,
                            'ec_numbers': ec_numbers,
                            'ec_count': len(ec_numbers)}
                    else:
                        print(f"Please enter a number between 1 and {len(modules)}")
                else:
                    print("Invalid input. Please try again.")
                    
            except KeyboardInterrupt:
                print("\nOperation cancelled by user.")
                return {}
            except Exception as e:
                print(f"Error: {e}. Please try again.")
    
    def save_results(self, results: Dict, filename: str = None):
        """
        Save extraction results to files for downstream analysis
        
        Args:
            results: Results dictionary from extract_pathway_components
            filename: Base filename (will add appropriate suffix)
        """
        if not filename:
            pathway_id = results.get('pathway_id', 'pathway').replace('ko', '').replace('map', '')
            filename = f"kegg_extraction_{pathway_id}"
        
        granularity = results.get('granularity', 'unknown')
        
        if granularity == 'pathway':
            # Save single EC number file
            ec_file = f"{filename}_ec_numbers.txt"
            with open(ec_file, 'w') as f:
                f.write(f"# EC numbers from KEGG pathway: {results['pathway_id']}\n")
                f.write(f"# Total EC numbers: {results['total_ec_count']}\n")
                f.write(f"# Extraction method: pathway-level\n\n")
                
                for ec in sorted(results['ec_numbers']):
                    f.write(f"{ec}\n")
            
            print(f"Saved {results['total_ec_count']} EC numbers to {ec_file}")
            
        elif granularity == 'modules':
            # Save separate files for each module
            for module_id, module_data in results['modules'].items():
                module_file = f"{filename}_{module_id}_ec_numbers.txt"
                
                with open(module_file, 'w') as f:
                    f.write(f"# EC numbers from KEGG module: {module_id}\n")
                    f.write(f"# Module name: {module_data['info']['name']}\n")
                    f.write(f"# Module definition: {module_data['info']['definition']}\n")
                    f.write(f"# Total EC numbers: {module_data['ec_count']}\n\n")
                    
                    for ec in sorted(module_data['ec_numbers']):
                        f.write(f"{ec}\n")
                
                print(f"Saved {module_data['ec_count']} EC numbers to {module_file}")
            
            # Also save a summary file
            summary_file = f"{filename}_modules_summary.txt"
            with open(summary_file, 'w') as f:
                f.write(f"# Module summary for pathway: {results['pathway_id']}\n")
                f.write(f"# Total modules: {results['total_modules']}\n")
                f.write(f"# Total EC numbers: {results['total_ec_count']}\n\n")
                
                for module_id, module_data in results['modules'].items():
                    f.write(f"{module_id}: {module_data['info']['name']} ({module_data['ec_count']} ECs)\n")
            
            print(f"Saved module summary to {summary_file}")
            
        elif granularity == 'single_module':
            # Save single module results
            module_info = results['selected_module']
            module_id = module_info['id']
            
            ec_file = f"{filename}_{module_id}_ec_numbers.txt"
            with open(ec_file, 'w') as f:
                f.write(f"# EC numbers from KEGG module: {module_id}\n")
                f.write(f"# Module name: {module_info['name']}\n")
                f.write(f"# Total EC numbers: {results['ec_count']}\n\n")
                
                for ec in sorted(results['ec_numbers']):
                    f.write(f"{ec}\n")
            
            print(f"Saved {results['ec_count']} EC numbers from {module_id} to {ec_file}")

def main():
    """
    Example usage and command-line interface
    """
    import sys
    
    print("KEGG Pathway Module Discovery Tool")
    print("=" * 40)
    
    # Get pathway ID from user
    if len(sys.argv) > 1:
        pathway_id = sys.argv[1]
    else:
        pathway_id = input("Enter KEGG pathway ID (e.g., ko00720, map00010, 00720): ").strip()
    
    # Initialize the tool
    discoverer = KEGGModuleDiscovery(delay=0.1)
    
    # Run interactive module selection
    results = discoverer.interactive_module_selection(pathway_id)

    if results:
        # Save results
        discoverer.save_results(results)
        
        print(f"\n{'='*50}")
        print("EXTRACTION COMPLETE!")
        print(f"{'='*50}")
        
        if 'ec_numbers' in results:
            print(f"Found {len(results['ec_numbers'])} EC numbers")
            if len(results['ec_numbers']) <= 20:  # Show if not too many
                print("EC numbers found:")
                for ec in sorted(results['ec_numbers']):
                    print(f"  {ec}")
        
        print(f"\nFiles saved. You can now run your Pfam analysis script!")


if __name__ == "__main__":
    main()