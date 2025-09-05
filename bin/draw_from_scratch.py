#!/usr/bin/env python3

import os
import re
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio.KEGG import REST
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Rectangle, FancyBboxPatch, Circle, FancyArrowPatch
import numpy as np
from PIL import Image
import requests
from io import BytesIO


class SimpleKEGGDrawer:
    def __init__(self):
        """
        Simple KEGG pathway drawer
        """
        self.cache_dir = Path("kegg_cache")
        self.cache_dir.mkdir(exist_ok=True)
        
    def get_pathway_kgml(self, pathway_id: str) -> ET.Element:
        """
        Get KEGG pathway data in KGML (XML) format
        
        Args:
            pathway_id: KEGG pathway ID (e.g., '00720')
            
        Returns:
            XML Element tree root
        """
        # Standardize pathway ID - use 'ko' prefix for reference pathway
        if not pathway_id.startswith('ko') and not pathway_id.startswith('map'):
            pathway_id = f'ko{pathway_id}'
        
        cache_file = self.cache_dir / f"{pathway_id}.xml"
        
        if cache_file.exists():
            print(f"Loading cached KGML for {pathway_id}")
            tree = ET.parse(cache_file)
            return tree.getroot()
        
        try:
            print(f"Downloading KGML for {pathway_id}...")
            kgml_data = REST.kegg_get(pathway_id, option='kgml').read()
            
            # Save to cache
            with open(cache_file, 'w') as f:
                f.write(kgml_data)
            
            # Parse XML
            root = ET.fromstring(kgml_data)
            return root
            
        except Exception as e:
            print(f"Error getting KGML for {pathway_id}: {e}")
            return None
    
    def parse_pathway_components(self, kgml_root: ET.Element) -> Dict:
        """
        Parse KGML to extract pathway components for drawing
        
        Args:
            kgml_root: KGML XML root element
            
        Returns:
            Dictionary with pathway components
        """
        components = {
            'title': kgml_root.get('title', 'KEGG Pathway'),
            'entries': {},
            'reactions': [],
            'relations': []
        }
        
        print(f"Pathway: {components['title']}")
        
        # Parse entries (enzymes, compounds, maps, etc.)
        for entry in kgml_root.findall('.//entry'):
            entry_id = entry.get('id')
            entry_type = entry.get('type')
            entry_name = entry.get('name', '')
            
            # Parse graphics information
            graphics = entry.find('graphics')
            if graphics is not None:
                components['entries'][entry_id] = {
                    'type': entry_type,
                    'name': entry_name,
                    'x': float(graphics.get('x', 0)),
                    'y': float(graphics.get('y', 0)),
                    'width': float(graphics.get('width', 46)),
                    'height': float(graphics.get('height', 17)),
                    'label': graphics.get('name', ''),
                    'shape': graphics.get('type', 'rectangle'),
                    'bgcolor': graphics.get('bgcolor', '#FFFFFF'),
                    'fgcolor': graphics.get('fgcolor', '#000000')
                }
                
                # Extract EC numbers if it's an enzyme
                if entry_type == 'enzyme':
                    ec_pattern = r'(\d+\.\d+\.\d+\.\d+)'
                    ec_numbers = re.findall(ec_pattern, entry_name)
                    components['entries'][entry_id]['ec_numbers'] = ec_numbers
        
        # Parse reactions
        for reaction in kgml_root.findall('.//reaction'):
            reaction_data = {
                'id': reaction.get('id'),
                'name': reaction.get('name'),
                'type': reaction.get('type'),
                'substrates': [],
                'products': []
            }
            
            for substrate in reaction.findall('substrate'):
                substrate_id = substrate.get('id')
                substrate_name = substrate.get('name')
                reaction_data['substrates'].append({'id': substrate_id, 'name': substrate_name})
            
            for product in reaction.findall('product'):
                product_id = product.get('id')
                product_name = product.get('name')
                reaction_data['products'].append({'id': product_id, 'name': product_name})
            
            components['reactions'].append(reaction_data)
        
        # Parse relations (arrows between entries)
        for relation in kgml_root.findall('.//relation'):
            relation_data = {
                'entry1': relation.get('entry1'),
                'entry2': relation.get('entry2'),
                'type': relation.get('type')
            }
            
            # Get subtype for arrow style
            subtype = relation.find('subtype')
            if subtype is not None:
                relation_data['subtype'] = subtype.get('name')
            
            components['relations'].append(relation_data)
        
        print(f"Parsed {len(components['entries'])} entries, {len(components['reactions'])} reactions, {len(components['relations'])} relations")
        
        return components
    
    def draw_pathway(self, pathway_id: str, output_file: str = None):
        """
        Draw the KEGG pathway
        
        Args:
            pathway_id: KEGG pathway ID (e.g., '00720')
            output_file: Output file name
        """
        if not output_file:
            output_file = f"pathway_{pathway_id}.png"
        
        # Get and parse KGML
        kgml = self.get_pathway_kgml(pathway_id)
        if kgml is None:
            print("Failed to retrieve pathway data")
            return None
        
        components = self.parse_pathway_components(kgml)
        
        # Create figure
        fig, ax = plt.subplots(figsize=(20, 16))
        
        # Set up the plot with KEGG coordinate system
        ax.set_xlim(0, 1200)
        ax.set_ylim(0, 900)
        ax.invert_yaxis()  # KEGG uses top-left origin
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Define colors for different entry types
        type_colors = {
            'enzyme': '#B4E7CE',      # Light green
            'compound': '#FFE5E5',     # Light pink
            'map': '#E5E5FF',         # Light blue
            'gene': '#FFFACD',        # Light yellow
            'group': '#F0F0F0',       # Light gray
            'ortholog': '#D4F1F4',    # Light cyan
            'other': '#FFFFFF'        # White
        }
        
        # Draw entries
        for entry_id, entry in components['entries'].items():
            x = entry['x']
            y = entry['y']
            width = entry['width']
            height = entry['height']
            entry_type = entry['type']
            shape = entry['shape']
            
            # Get color
            if entry['bgcolor'] != '#FFFFFF':
                color = entry['bgcolor']
            else:
                color = type_colors.get(entry_type, type_colors['other'])
            
            # Draw based on shape
            if shape == 'circle':
                # Draw circle for compounds
                circle = Circle((x, y), min(width, height)/2,
                              facecolor=color,
                              edgecolor='black',
                              linewidth=1,
                              alpha=0.8)
                ax.add_patch(circle)
                
            elif shape == 'roundrectangle' or entry_type == 'map':
                # Draw rounded rectangle
                rect = FancyBboxPatch(
                    (x - width/2, y - height/2), width, height,
                    boxstyle="round,pad=3",
                    facecolor=color,
                    edgecolor='black',
                    linewidth=1.5 if entry_type == 'map' else 1,
                    alpha=0.8
                )
                ax.add_patch(rect)
                
            else:
                # Draw regular rectangle (default)
                rect = Rectangle(
                    (x - width/2, y - height/2), width, height,
                    facecolor=color,
                    edgecolor='black',
                    linewidth=1,
                    alpha=0.8
                )
                ax.add_patch(rect)
            
            # Add label
            label = entry['label']
            if label:
                # Truncate long labels
                if len(label) > 15 and entry_type != 'map':
                    label = label[:12] + '...'
                
                # Adjust font size based on entry type
                fontsize = 10 if entry_type == 'map' else 8
                weight = 'bold' if entry_type == 'map' else 'normal'
                
                ax.text(x, y, label,
                       ha='center', va='center',
                       fontsize=fontsize,
                       weight=weight,
                       wrap=True)
        
        # Draw relations (arrows)
        for relation in components['relations']:
            entry1_id = relation['entry1']
            entry2_id = relation['entry2']
            
            if entry1_id in components['entries'] and entry2_id in components['entries']:
                entry1 = components['entries'][entry1_id]
                entry2 = components['entries'][entry2_id]
                
                x1, y1 = entry1['x'], entry1['y']
                x2, y2 = entry2['x'], entry2['y']
                
                # Determine arrow style based on relation type
                if relation.get('subtype') == 'inhibition':
                    arrowstyle = '-|'
                    color = 'red'
                    alpha = 0.6
                elif relation.get('subtype') == 'activation':
                    arrowstyle = '->'
                    color = 'green'
                    alpha = 0.6
                else:
                    arrowstyle = '->'
                    color = 'gray'
                    alpha = 0.4
                
                # Draw arrow
                arrow = FancyArrowPatch(
                    (x1, y1), (x2, y2),
                    arrowstyle=arrowstyle,
                    color=color,
                    linewidth=1.5,
                    alpha=alpha,
                    connectionstyle="arc3,rad=0.1"
                )
                ax.add_patch(arrow)
        
        # Draw reactions (connect compounds and enzymes)
        for reaction in components['reactions']:
            # This is simplified - in reality you'd need to map reaction components to entries
            pass
        
        # Add title
        title = components['title']
        ax.text(600, 30, title, fontsize=16, weight='bold', ha='center')
        
        # Add legend
        legend_elements = []
        for entry_type, color in type_colors.items():
            if any(e['type'] == entry_type for e in components['entries'].values()):
                label = entry_type.capitalize()
                legend_elements.append(patches.Patch(facecolor=color, edgecolor='black', label=label))
        
        if legend_elements:
            ax.legend(handles=legend_elements, loc='upper left', fontsize=10)
        
        # Add statistics
        stats_text = f"Entries: {len(components['entries'])}\n"
        stats_text += f"Reactions: {len(components['reactions'])}\n"
        stats_text += f"Relations: {len(components['relations'])}"
        
        ax.text(0.98, 0.02, stats_text,
               transform=ax.transAxes,
               fontsize=10,
               verticalalignment='bottom',
               horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Save figure
        plt.tight_layout()
        plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"\nPathway diagram saved to: {output_file}")
        return output_file
    
    def download_kegg_image(self, pathway_id: str, output_file: str = None):
        """
        Download the official KEGG pathway image
        
        Args:
            pathway_id: KEGG pathway ID (e.g., '00720')
            output_file: Output file name
        """
        if not output_file:
            output_file = f"kegg_pathway_{pathway_id}.png"
        
        # Format pathway ID for URL
        if not pathway_id.startswith('map'):
            pathway_id = f'map{pathway_id.replace("ko", "").replace("map", "")}'
        
        url = f"https://www.kegg.jp/kegg/pathway/{pathway_id}.png"
        
        try:
            print(f"Downloading official KEGG image from {url}...")
            response = requests.get(url)
            
            if response.status_code == 200:
                img = Image.open(BytesIO(response.content))
                img.save(output_file)
                print(f"Official KEGG image saved to: {output_file}")
                return output_file
            else:
                print(f"Failed to download: HTTP {response.status_code}")
                return None
                
        except Exception as e:
            print(f"Error downloading: {e}")
            return None
    
    def get_pathway_info(self, pathway_id: str):
        """
        Get basic information about a pathway
        
        Args:
            pathway_id: KEGG pathway ID
        """
        # Standardize pathway ID
        if not pathway_id.startswith('ko') and not pathway_id.startswith('map'):
            pathway_id = f'ko{pathway_id}'
        
        try:
            print(f"\nGetting pathway information for {pathway_id}...")
            pathway_data = REST.kegg_get(pathway_id).read()
            
            # Parse basic info
            lines = pathway_data.split('\n')
            info = {}
            
            for line in lines:
                if line.startswith('NAME'):
                    info['name'] = line.replace('NAME', '').strip()
                elif line.startswith('DESCRIPTION'):
                    info['description'] = line.replace('DESCRIPTION', '').strip()
                elif line.startswith('CLASS'):
                    info['class'] = line.replace('CLASS', '').strip()
                elif line.startswith('MODULE'):
                    if 'modules' not in info:
                        info['modules'] = []
                    module = line.replace('MODULE', '').strip()
                    info['modules'].append(module)
            
            # Print info
            print("\n" + "="*60)
            print(f"Pathway: {pathway_id}")
            if 'name' in info:
                print(f"Name: {info['name']}")
            if 'description' in info:
                print(f"Description: {info['description']}")
            if 'class' in info:
                print(f"Class: {info['class']}")
            if 'modules' in info:
                print(f"Associated modules: {len(info['modules'])}")
                for module in info['modules'][:5]:  # Show first 5
                    print(f"  - {module}")
                if len(info['modules']) > 5:
                    print(f"  ... and {len(info['modules'])-5} more")
            print("="*60)
            
            return info
            
        except Exception as e:
            print(f"Error getting pathway info: {e}")
            return {}


def main():
    """
    Main function to draw pathway 00720
    """
    print("KEGG Pathway Drawer")
    print("="*60)
    
    # Initialize drawer
    drawer = SimpleKEGGDrawer()
    
    # Pathway to draw (Carbon fixation pathways in prokaryotes)
    pathway_id = "00720"
    
    # Get pathway information
    info = drawer.get_pathway_info(pathway_id)
    
    # Draw the pathway
    print(f"\nDrawing pathway {pathway_id}...")
    output_file = drawer.draw_pathway(pathway_id)
    
    # Also download the official KEGG image for comparison
    kegg_file = drawer.download_kegg_image(pathway_id)
    
    print("\n" + "="*60)
    print("Done! Generated files:")
    if output_file:
        print(f"  - Custom drawing: {output_file}")
    if kegg_file:
        print(f"  - Official KEGG image: {kegg_file}")
    print("="*60)


if __name__ == "__main__":
    main()