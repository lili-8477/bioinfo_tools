#!/usr/bin/env python3
"""
Generalized Gene Homolog Conversion Tool

Command-line tool for converting gene symbols between mouse and human
using GSEAPY Biomart with robust error handling and naming convention fallback.

Author: Generated with Claude Code
Date: 2025-08-21
"""

import pandas as pd
from typing import Dict, List, Union
from gseapy import Biomart
import logging
import argparse
import sys
import os

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class HomologConverter:
    """
    A class for converting gene symbols between mouse and human using GSEAPY Biomart.
    """
    
    def __init__(self, timeout: int = 300, verbose: bool = True):
        """
        Initialize the HomologConverter.
        
        Args:
            timeout: Query timeout in seconds
            verbose: Whether to print progress messages
        """
        self.timeout = timeout
        self.verbose = verbose
        self.biomart = None
        
    def test_biomart_connection(self) -> bool:
        """
        Test Biomart connection with a simple query.
        
        Returns:
            bool: True if connection successful, False otherwise
        """
        try:
            if self.verbose:
                print("Testing Biomart connection...")
                
            bm = Biomart()
            bm.query(
                dataset='mmusculus_gene_ensembl',
                attributes=['external_gene_name'],
                filters={'external_gene_name': ['Actb']}
            )
            
            if self.verbose:
                print("✓ Biomart connection successful")
            return True
            
        except Exception as e:
            if self.verbose:
                print(f"✗ Biomart connection failed: {e}")
            return False
    
    def get_dataset_config(self, direction: str) -> Dict[str, str]:
        """
        Get dataset configuration for the specified conversion direction.
        
        Args:
            direction: 'm2h' for mouse-to-human or 'h2m' for human-to-mouse
            
        Returns:
            Dict containing dataset and attribute configurations
            
        Raises:
            ValueError: If direction is not 'm2h' or 'h2m'
        """
        if direction == 'm2h':
            return {
                'source_dataset': 'mmusculus_gene_ensembl',
                'source_attr': 'external_gene_name',
                'target_attr': 'hsapiens_homolog_associated_gene_name',
                'test_gene': 'Actb',
                'source_species': 'mouse',
                'target_species': 'human'
            }
        elif direction == 'h2m':
            return {
                'source_dataset': 'hsapiens_gene_ensembl', 
                'source_attr': 'external_gene_name',
                'target_attr': 'mmusculus_homolog_associated_gene_name',
                'test_gene': 'ACTB',
                'source_species': 'human',
                'target_species': 'mouse'
            }
        else:
            raise ValueError("Direction must be 'm2h' (mouse-to-human) or 'h2m' (human-to-mouse)")
    
    
    def apply_naming_convention(self, gene: str, direction: str) -> str:
        """
        Apply species-specific naming conventions.
        
        Args:
            gene: Gene symbol
            direction: 'm2h' or 'h2m'
            
        Returns:
            Gene symbol with appropriate naming convention
        """
        if direction == 'm2h':
            # Mouse to human: convert to uppercase
            return gene.upper()
        elif direction == 'h2m':
            # Human to mouse: first letter uppercase, rest lowercase
            return gene.capitalize()
        else:
            return gene
    
    def convert_gene_orthologs(
        self, 
        genes: Union[str, List[str]], 
        direction: str = 'm2h',
        use_biomart: bool = True,
        fallback_to_convention: bool = True
    ) -> Dict[str, Union[str, Dict]]:
        """
        Convert gene symbols between mouse and human using GSEAPY Biomart.
        
        Args:
            genes: Single gene string or list of gene symbols to convert
            direction: 'm2h' for mouse-to-human or 'h2m' for human-to-mouse
            use_biomart: Whether to attempt Biomart query
            fallback_to_convention: Whether to use naming convention for unmapped genes
            
        Returns:
            Dictionary containing:
            - 'mapping': Dict of gene mappings {source_gene: target_gene}
            - 'statistics': Dict with conversion statistics
            - 'unmapped': List of genes that couldn't be mapped
            - 'source': Dict indicating mapping source for each gene
        """
        # Ensure genes is a list
        if isinstance(genes, str):
            genes = [genes]
        
        # Get configuration for the specified direction
        config = self.get_dataset_config(direction)
        
        # Initialize results
        mapping = {}
        unmapped = []
        mapping_source = {}
        biomart_mapping = {}
        
        if self.verbose:
            print(f"Converting {len(genes)} genes from {config['source_species']} to {config['target_species']}")
        
        # Try Biomart approach if requested
        if use_biomart:
            biomart_mapping = self._query_biomart(genes, config)
            mapping.update(biomart_mapping)
            
            # Mark biomart genes
            for gene in biomart_mapping:
                mapping_source[gene] = 'biomart'
        
        # Apply naming convention for unmapped genes
        if fallback_to_convention:
            for gene in genes:
                if gene not in mapping:
                    # Apply naming convention
                    converted = self.apply_naming_convention(gene, direction)
                    mapping[gene] = converted
                    mapping_source[gene] = 'naming_convention'
        
        # Identify unmapped genes
        unmapped = [gene for gene in genes if gene not in mapping]
        
        # Calculate statistics
        stats = {
            'total_input': len(genes),
            'biomart_mapped': len(biomart_mapping),
            'convention_mapped': len([g for g, s in mapping_source.items() if s == 'naming_convention']),
            'unmapped': len(unmapped),
            'success_rate': (len(mapping) / len(genes)) * 100 if genes else 0
        }
        
        if self.verbose:
            self._print_statistics(stats, config)
        
        return {
            'mapping': mapping,
            'statistics': stats,
            'unmapped': unmapped,
            'source': mapping_source
        }
    
    def _query_biomart(self, genes: List[str], config: Dict[str, str]) -> Dict[str, str]:
        """
        Query Biomart for gene ortholog mappings.
        
        Args:
            genes: List of genes to map
            config: Dataset configuration
            
        Returns:
            Dictionary of gene mappings from Biomart
        """
        biomart_mapping = {}
        
        try:
            if self.verbose:
                print("Attempting GSEAPY Biomart conversion...")
                
            bm = Biomart()
            
            # Test connection first
            if not self.test_biomart_connection():
                return {}
            
            # Query orthologs
            if self.verbose:
                print(f"Querying {config['source_species']}-{config['target_species']} orthologs...")
                
            result = bm.query(
                dataset=config['source_dataset'],
                attributes=[config['source_attr'], config['target_attr']]
            )
            
            if result.empty:
                if self.verbose:
                    print("Empty result from Biomart query")
                return {}
            
            if self.verbose:
                print(f"Retrieved {len(result)} ortholog mappings from Biomart")
            
            # Create mapping dictionary with data validation
            all_biomart_mapping = {}
            skipped_count = 0
            
            for _, row in result.iterrows():
                source_gene = row[config['source_attr']]
                target_gene = row[config['target_attr']]
                
                # Only add if both genes are valid
                if (pd.notna(source_gene) and pd.notna(target_gene) and
                    str(target_gene).strip() and str(source_gene).strip()):
                    all_biomart_mapping[str(source_gene)] = str(target_gene)
                else:
                    skipped_count += 1
            
            if self.verbose:
                print(f"Created mapping dictionary with {len(all_biomart_mapping)} entries")
                print(f"Skipped {skipped_count} invalid mappings")
            
            # Filter for input genes
            for gene in genes:
                if str(gene) in all_biomart_mapping:
                    biomart_mapping[gene] = all_biomart_mapping[str(gene)]
            
            if self.verbose:
                print(f"Successfully mapped {len(biomart_mapping)} out of {len(genes)} input genes")
                
                # Show examples
                if biomart_mapping:
                    print("\nExample mappings from Biomart:")
                    for source, target in list(biomart_mapping.items())[:5]:
                        print(f"  {source} -> {target}")
            
            return biomart_mapping
            
        except ImportError as e:
            if self.verbose:
                print(f"GSEAPY import error: {e}")
                print("Please install gseapy: pip install gseapy")
        except ConnectionError as e:
            if self.verbose:
                print(f"Network connection error: {e}")
        except Exception as e:
            if self.verbose:
                print(f"GSEAPY Biomart failed with error: {e}")
                print("Falling back to naming convention")
        
        return {}
    
    def _print_statistics(self, stats: Dict, config: Dict[str, str]):
        """Print conversion statistics."""
        print("\n=== CONVERSION STATISTICS ===")
        print(f"Direction: {config['source_species']} -> {config['target_species']}")
        print(f"Total input genes: {stats['total_input']}")
        print(f"Biomart mapped: {stats['biomart_mapped']}")
        print(f"Convention mapped: {stats['convention_mapped']}")
        print(f"Unmapped: {stats['unmapped']}")
        print(f"Success rate: {stats['success_rate']:.1f}%")


# Convenience functions
def convert_mouse_to_human(genes: Union[str, List[str]], **kwargs) -> Dict:
    """
    Convert mouse gene symbols to human homologs.
    
    Args:
        genes: Mouse gene symbols to convert
        **kwargs: Additional arguments passed to convert_gene_orthologs
        
    Returns:
        Conversion results dictionary
    """
    converter = HomologConverter(**kwargs)
    return converter.convert_gene_orthologs(genes, direction='m2h')


def convert_human_to_mouse(genes: Union[str, List[str]], **kwargs) -> Dict:
    """
    Convert human gene symbols to mouse homologs.
    
    Args:
        genes: Human gene symbols to convert
        **kwargs: Additional arguments passed to convert_gene_orthologs
        
    Returns:
        Conversion results dictionary
    """
    converter = HomologConverter(**kwargs)
    return converter.convert_gene_orthologs(genes, direction='h2m')


def create_parser():
    """Create and configure argument parser."""
    parser = argparse.ArgumentParser(
        description='Convert gene symbols between mouse and human using GSEAPY Biomart',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert mouse genes to human
  %(prog)s -g Actb Kif11 Tpx2 -d m2h
  
  # Convert human genes to mouse
  %(prog)s -g ACTB KIF11 TPX2 -d h2m
  
  # Read genes from file
  %(prog)s -f mouse_genes.txt -d m2h -o human_orthologs.txt
  
  # Use only Biomart (no fallback)
  %(prog)s -g Actb -d m2h --no-fallback
        """
    )
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '-g', '--genes',
        nargs='+',
        help='Gene symbols to convert (space-separated)'
    )
    input_group.add_argument(
        '-f', '--file',
        help='File containing gene symbols (one per line)'
    )
    
    parser.add_argument(
        '-d', '--direction',
        choices=['m2h', 'h2m'],
        required=True,
        help='Conversion direction: m2h (mouse-to-human) or h2m (human-to-mouse)'
    )
    
    parser.add_argument(
        '-o', '--output',
        help='Output file for results (default: stdout)'
    )
    
    parser.add_argument(
        '--no-biomart',
        action='store_true',
        help='Skip Biomart query, use only naming convention'
    )
    
    parser.add_argument(
        '--no-fallback',
        action='store_true',
        help='Skip naming convention fallback for unmapped genes'
    )
    
    parser.add_argument(
        '--format',
        choices=['tsv', 'json', 'simple'],
        default='tsv',
        help='Output format (default: tsv)'
    )
    
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='Suppress progress messages'
    )
    
    parser.add_argument(
        '--timeout',
        type=int,
        default=300,
        help='Biomart query timeout in seconds (default: 300)'
    )
    
    return parser


def read_genes_from_file(filepath: str) -> List[str]:
    """Read gene symbols from file, one per line."""
    try:
        with open(filepath, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
        return genes
    except FileNotFoundError:
        print(f"Error: File '{filepath}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file '{filepath}': {e}", file=sys.stderr)
        sys.exit(1)


def format_output(result: Dict, format_type: str) -> str:
    """Format conversion results for output."""
    if format_type == 'simple':
        lines = []
        for source, target in result['mapping'].items():
            lines.append(f"{source}\t{target}")
        return '\n'.join(lines)
    
    elif format_type == 'json':
        import json
        return json.dumps(result, indent=2)
    
    else:  # tsv
        lines = ["Source_Gene\tTarget_Gene\tMapping_Source"]
        for source, target in result['mapping'].items():
            source_type = result['source'].get(source, 'unknown')
            lines.append(f"{source}\t{target}\t{source_type}")
        return '\n'.join(lines)


def main():
    """Main command-line interface."""
    parser = create_parser()
    args = parser.parse_args()
    
    # Read input genes
    if args.genes:
        genes = args.genes
    else:
        genes = read_genes_from_file(args.file)
    
    if not genes:
        print("Error: No genes provided", file=sys.stderr)
        sys.exit(1)
    
    # Set up converter
    converter = HomologConverter(
        timeout=args.timeout,
        verbose=not args.quiet
    )
    
    # Perform conversion
    try:
        result = converter.convert_gene_orthologs(
            genes=genes,
            direction=args.direction,
            use_biomart=not args.no_biomart,
            fallback_to_convention=not args.no_fallback
        )
        
        # Format output
        output_text = format_output(result, args.format)
        
        # Write output
        if args.output:
            try:
                with open(args.output, 'w') as f:
                    f.write(output_text)
                if not args.quiet:
                    print(f"Results written to {args.output}")
            except Exception as e:
                print(f"Error writing to file '{args.output}': {e}", file=sys.stderr)
                sys.exit(1)
        else:
            print(output_text)
        
        # Exit with appropriate code
        unmapped_count = result['statistics']['unmapped']
        if unmapped_count > 0 and not args.quiet:
            print(f"Warning: {unmapped_count} genes could not be mapped", file=sys.stderr)
        
        sys.exit(0 if unmapped_count == 0 else 1)
        
    except Exception as e:
        print(f"Error during conversion: {e}", file=sys.stderr)
        sys.exit(1)


# Example usage and testing
if __name__ == "__main__":
    main()