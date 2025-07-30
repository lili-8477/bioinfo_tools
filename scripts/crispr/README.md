# CRISPR/RNA Analysis Tools

Tools for CRISPR guide RNA design and RNA sequence analysis.

## Tools

### PFS_scanner/
Python tool for finding target sequences in RNA with customizable pattern matching.

**Features:**
- Flexible pattern matching with wildcards (N=[ACGU], V=[ACG])
- Guide RNA sequence generation with reverse complement
- Interactive command-line interface
- RNA sequence validation

**Files:**
- **PFS_Scanner.py** - Main scanner script
- **test.ipynb** - Jupyter notebook with examples and testing
- **GFP_Seq.txt** - Sample sequence file for testing

**Usage:**
```bash
cd PFS_scanner
python PFS_Scanner.py
```

**Interactive prompts:**
1. Enter RNA sequence (5' → 3')
2. Enter target pattern (e.g., "NAANNGC")
3. Get guide RNA sequences

**Pattern Syntax:**
- N = any nucleotide [ACGU]
- V = any nucleotide except U [ACG]
- Specific nucleotides: A, C, G, U

**Example:**
```
Enter the RNA sequence (5' => 3'): GGAUCCAUGCUGCACGCC...
Enter the pattern: NAANNGC
```

## Dependencies

- **re** - Regular expression matching (standard library)
- **jupyter** - For running test notebook (optional)

## Output

The tool outputs guide RNA sequences (reverse complement of target regions minus PAM sequence) that can be used for CRISPR design.