#!/usr/bin/env python3

import os
import re
import time
from pathlib import Path
from typing import Dict, List, Set, Optional
from Bio.KEGG import REST
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Rectangle, Circle
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import requests
from io import BytesIO


class KEGGPathwayAnnotator:
    def __init__(self):
        """
        Download and annotate official KEGG pathway images
        """
        self.cache_dir = Path("kegg_cache")
        self.cache_dir.mkdir(exist_ok=True)
        
    def download_pathway_image(self, pathway_id: str, output_file: str = None):
        """
        Download the official KEGG pathway image
        
        Args:
            pathway_id: KEGG pathway ID (e.g., '00720')
            output_file: Output file name
            
        Returns:
            Path to downloaded image
        """
        if not output_file:
            output_file = f"kegg_pathway_{pathway_id}.png"
        
        # Format pathway ID for URL - use 'map' prefix for the image
        clean_id = pathway_id.replace("ko", "").replace("map", "").replace("ec", "")
        map_id = f'map{clean_id}'
        
        # Try different URL formats
        urls = [
            f"https://www.kegg.jp/kegg/pathway/{map_id}.png",
            f"https://www.kegg.jp/pathway/{map_id}.png",
            f"https://rest.kegg.jp/get/{map_id}/image"
        ]
        
        for url in urls:
            try:
                print(f"Trying to download from: {url}")
                response = requests.get(url, headers={'User-Agent': 'Mozilla/5.0'})
                
                if response.status_code == 200:
                    # Check if it's actually an image
                    if 'image' in response.headers.get('content-type', ''):
                        img = Image.open(BytesIO(response.content))
                        img.save(output_file)
                        print(f"✓ KEGG pathway image saved to: {output_file}")
                        return output_file
                    
            except Exception as e:
                print(f"  Failed: {e}")
                continue
        
        print(f"Could not download pathway image for {pathway_id}")
        return None
    
    def get_pathway_info(self, pathway_id: str) -> Dict:
        """
        Get pathway information including EC numbers
        
        Args:
            pathway_id: KEGG pathway ID
            
        Returns:
            Dictionary with pathway information
        """
        # Try with 'ko' prefix for reference pathway
        if not pathway_id.startswith('ko') and not pathway_id.startswith('map'):
            pathway_id = f'ko{pathway_id}'
        
        try:
            print(f"\nGetting pathway information for {pathway_id}...")
            pathway_data = REST.kegg_get(pathway_id).read()
            
            info = {
                'name': '',
                'description': '',
                'ec_numbers': set(),
                'ko_numbers': set(),
                'compounds': set(),
                'modules': []
            }
            
            lines = pathway_data.split('\n')
            current_section = None
            
            for line in lines:
                if line.startswith('NAME'):
                    info['name'] = line.replace('NAME', '').strip()
                elif line.startswith('DESCRIPTION'):
                    info['description'] = line.replace('DESCRIPTION', '').strip()
                elif line.startswith('CLASS'):
                    info['class'] = line.replace('CLASS', '').strip()
                elif line.startswith('MODULE'):
                    current_section = 'MODULE'
                elif line.startswith('ENZYME'):
                    current_section = 'ENZYME'
                elif line.startswith('ORTHOLOGY'):
                    current_section = 'ORTHOLOGY'
                elif line.startswith('COMPOUND'):
                    current_section = 'COMPOUND'
                elif line.startswith('REFERENCE') or line.startswith('///'):
                    current_section = None
                elif current_section and line.strip():
                    if current_section == 'ENZYME':
                        # Extract EC numbers
                        ec_numbers = re.findall(r'(\d+\.\d+\.\d+\.\d+)', line)
                        info['ec_numbers'].update(ec_numbers)
                    elif current_section == 'MODULE':
                        # Extract module IDs
                        modules = re.findall(r'(M\d{5})', line)
                        info['modules'].extend(modules)
                    elif current_section == 'ORTHOLOGY':
                        # Extract KO numbers
                        ko_numbers = re.findall(r'(K\d{5})', line)
                        info['ko_numbers'].update(ko_numbers)
                    elif current_section == 'COMPOUND':
                        # Extract compound IDs
                        compounds = re.findall(r'(C\d{5})', line)
                        info['compounds'].update(compounds)
            
            # Print summary
            print("\n" + "="*60)
            print(f"Pathway: {pathway_id}")
            print(f"Name: {info['name']}")
            if info['description']:
                print(f"Description: {info['description']}")
            print(f"EC numbers found: {len(info['ec_numbers'])}")
            print(f"KO numbers found: {len(info['ko_numbers'])}")
            print(f"Compounds found: {len(info['compounds'])}")
            print(f"Modules found: {len(info['modules'])}")
            
            # Show some EC numbers
            if info['ec_numbers']:
                ec_list = sorted(info['ec_numbers'])[:10]
                print(f"\nFirst EC numbers: {', '.join(ec_list)}")
                if len(info['ec_numbers']) > 10:
                    print(f"  ... and {len(info['ec_numbers'])-10} more")
            
            print("="*60)
            
            return info
            
        except Exception as e:
            print(f"Error getting pathway info: {e}")
            return {}
    
    def annotate_pathway_image(self, image_path: str, coverage_data: Dict = None, 
                              output_file: str = None):
        """
        Annotate the KEGG pathway image with coverage information
        
        Args:
            image_path: Path to KEGG pathway image
            coverage_data: Optional coverage data (EC -> present/absent)
            output_file: Output file name for annotated image
            
        Returns:
            Path to annotated image
        """
        if not output_file:
            output_file = image_path.replace('.png', '_annotated.png')
        
        try:
            # Open the image
            img = Image.open(image_path)
            
            # Create a figure with matplotlib for better annotation control
            fig, ax = plt.subplots(figsize=(20, 16))
            ax.imshow(img)
            ax.axis('off')
            
            # Add title
            ax.text(img.width/2, 30, 'KEGG Pathway Analysis', 
                   fontsize=20, weight='bold', ha='center',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            if coverage_data:
                # Add legend for coverage
                legend_elements = [
                    patches.Patch(color='green', label='Present in genome'),
                    patches.Patch(color='red', label='Missing from genome'),
                    patches.Patch(color='yellow', label='Partial coverage')
                ]
                ax.legend(handles=legend_elements, loc='upper right', fontsize=12)
                
                # Add coverage statistics
                total = len(coverage_data)
                present = sum(1 for v in coverage_data.values() if v == 'present')
                missing = sum(1 for v in coverage_data.values() if v == 'missing')
                partial = sum(1 for v in coverage_data.values() if v == 'partial')
                
                stats_text = f"Coverage Statistics:\n"
                stats_text += f"Total components: {total}\n"
                stats_text += f"Present: {present} ({present/total*100:.1f}%)\n"
                stats_text += f"Missing: {missing} ({missing/total*100:.1f}%)\n"
                stats_text += f"Partial: {partial} ({partial/total*100:.1f}%)"
                
                ax.text(10, img.height-10, stats_text,
                       fontsize=11, va='bottom',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
            
            # Save annotated image
            plt.tight_layout()
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            plt.close()
            
            print(f"✓ Annotated image saved to: {output_file}")
            return output_file
            
        except Exception as e:
            print(f"Error annotating image: {e}")
            return None
    
    def highlight_pathway_image(self, image_path: str, ec_numbers: Set[str], 
                               present_ecs: Set[str], output_file: str = None):
        """
        Create a highlighted version showing which ECs are present
        
        Args:
            image_path: Path to KEGG pathway image
            ec_numbers: All EC numbers in pathway
            present_ecs: EC numbers present in genome
            output_file: Output file name
            
        Returns:
            Path to highlighted image
        """
        if not output_file:
            output_file = image_path.replace('.png', '_highlighted.png')
        
        try:
            # Open image
            img = Image.open(image_path)
            
            # Convert to RGBA for transparency support
            img = img.convert('RGBA')
            
            # Create overlay for highlighting
            overlay = Image.new('RGBA', img.size, (255, 255, 255, 0))
            draw = ImageDraw.Draw(overlay)
            
            # Since we can't get exact EC positions from the image,
            # we'll add a summary panel
            panel_height = 150
            panel = Image.new('RGBA', (img.width, panel_height), (255, 255, 255, 240))
            panel_draw = ImageDraw.Draw(panel)
            
            # Add coverage summary to panel
            total = len(ec_numbers)
            present = len(present_ecs)
            missing = len(ec_numbers - present_ecs)
            
            try:
                # Try to use a better font if available
                font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 16)
                small_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 12)
            except:
                font = ImageFont.load_default()
                small_font = font
            
            # Title
            panel_draw.text((10, 10), f"Pathway Coverage Analysis", fill='black', font=font)
            
            # Statistics
            y_offset = 40
            panel_draw.text((10, y_offset), f"Total EC numbers: {total}", fill='black', font=small_font)
            panel_draw.text((10, y_offset+20), f"Present in genome: {present} ({present/total*100:.1f}%)", 
                          fill='green', font=small_font)
            panel_draw.text((10, y_offset+40), f"Missing from genome: {missing} ({missing/total*100:.1f}%)", 
                          fill='red', font=small_font)
            
            # Show some present ECs
            if present_ecs:
                present_list = sorted(present_ecs)[:8]
                ec_text = "Present: " + ", ".join(present_list)
                if len(present_ecs) > 8:
                    ec_text += f" ... (+{len(present_ecs)-8} more)"
                panel_draw.text((10, y_offset+70), ec_text, fill='darkgreen', font=small_font)
            
            # Show some missing ECs
            missing_ecs = ec_numbers - present_ecs
            if missing_ecs:
                missing_list = sorted(missing_ecs)[:8]
                ec_text = "Missing: " + ", ".join(missing_list)
                if len(missing_ecs) > 8:
                    ec_text += f" ... (+{len(missing_ecs)-8} more)"
                panel_draw.text((10, y_offset+90), ec_text, fill='darkred', font=small_font)
            
            # Combine original image with panel
            combined = Image.new('RGBA', (img.width, img.height + panel_height))
            combined.paste(img, (0, 0))
            combined.paste(panel, (0, img.height), panel)
            
            # Convert back to RGB for saving as PNG
            combined = combined.convert('RGB')
            combined.save(output_file)
            
            print(f"✓ Highlighted image saved to: {output_file}")
            return output_file
            
        except Exception as e:
            print(f"Error creating highlighted image: {e}")
            return None


def main():
    """
    Main function to download and annotate KEGG pathway 00720
    """
    print("KEGG Pathway Image Downloader and Annotator")
    print("="*60)
    
    # Initialize annotator
    annotator = KEGGPathwayAnnotator()
    
    # Pathway ID for Carbon fixation pathways in prokaryotes
    pathway_id = "00720"
    
    # Get pathway information
    info = annotator.get_pathway_info(pathway_id)
    
    # Download the official KEGG pathway image
    image_file = annotator.download_pathway_image(pathway_id)
    
    if image_file and info.get('ec_numbers'):
        # Example: Create highlighted version with mock coverage data
        # In real use, this would come from your pyhmmer results
        all_ecs = info['ec_numbers']
        
        # Simulate some ECs being present (in real use, from your genome analysis)
        import random
        present_ecs = set(random.sample(list(all_ecs), min(len(all_ecs)//2, len(all_ecs))))
        
        # Create highlighted version
        highlighted = annotator.highlight_pathway_image(
            image_file, 
            all_ecs,
            present_ecs,
            f"pathway_{pathway_id}_coverage.png"
        )
        
        print("\n" + "="*60)
        print("Analysis complete! Generated files:")
        print(f"  1. Original KEGG image: {image_file}")
        if highlighted:
            print(f"  2. Coverage analysis: {highlighted}")
        print("="*60)
    else:
        print("\nFailed to download pathway image or extract EC numbers")


if __name__ == "__main__":
    main()