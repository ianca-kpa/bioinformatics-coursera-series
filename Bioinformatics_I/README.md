
# Bioinformatics I — Finding Hidden Messages in DNA

A collection of Python implementations and Jupyter notebooks developed during the *Finding Hidden Messages in DNA* course of the Bioinformatics Specialization offered by the University of California San Diego on Coursera.

This collection explores fundamental algorithms for DNA sequence analysis, pattern identification, k-mer frequency analysis, and motif discovery.

## Repository Structure

### 01 — Pattern Matching

[Explore notebooks](./01-pattern-matching/)

Fundamental algorithms for identifying and comparing patterns in DNA sequences.

Topics include:

- Exact pattern matching
- Pattern counting and locating
- Reverse complements
- Hamming distance
- Approximate pattern matching

### 02 — Frequency Analysis

[Explore notebooks](./02-frequency-analysis/)

Computational approaches to identifying frequently occurring patterns and analyzing nucleotide distributions.

Topics include:

- K-mer frequency counting
- Frequency maps
- Clump finding
- Symbol arrays
- GC skew analysis

### 03 — Mismatches and Neighbors

[Explore notebooks](./03-mismatches-and-neighbors/)

Algorithms for analyzing DNA patterns while allowing nucleotide mismatches.

Topics include:

- Approximate pattern matching
- K-mer neighborhoods
- Frequent patterns with mismatches
- Reverse complement analysis

### 04 — Motif Search

[Explore notebooks](./04-motif-search/)

Algorithms for identifying shared sequence patterns across collections of DNA sequences.

Topics include:

- Motif enumeration
- Greedy motif search
- Randomized motif search
- Gibbs sampling
- Motif scoring

### 05 — Probability Tools

[Explore notebooks](./05-probability-tools/)

Probability-based methods used in sequence analysis and motif discovery.

Topics include:

- Probability calculations
- Profile-based sequence analysis
- Weighted random selection
- Probabilistic motif identification

### 06 — Motif Utilities

[Explore notebooks](./06-motif-utilities/)

Supporting functions for representing, evaluating, and comparing DNA motifs.

Topics include:

- Profile matrices
- Consensus sequences
- Motif scoring
- Pseudocounts
- Sequence probability calculations

## Technical Concepts

The notebooks explore computational concepts relevant to biological sequence analysis, including:

- Python programming
- String processing algorithms
- K-mer analysis
- Sequence similarity and mismatch calculations
- Frequency analysis
- Algorithmic motif discovery
- Probabilistic methods
- Biological sequence interpretation

## Tools and Environment

- Python 3
- Jupyter Notebook or JupyterLab
- Python standard library

Individual notebooks may require additional dependencies or associated input datasets.

## Running the Notebooks

Clone the repository:

```bash
git clone https://github.com/ianca-kpa/bioinformatics-coursera-series.git
cd bioinformatics-coursera-series
```

Open the notebooks using Jupyter Notebook, JupyterLab, or Visual Studio Code with the Jupyter extension.

The notebooks are organized by topic and can be explored individually.

Some exercises use external input files stored in the `dataset/` directory.

When executing these notebooks, ensure that relative file paths correspond to the current Jupyter working directory.

## Educational Context

These implementations originated from coursework and personal programming practice during the UC San Diego Bioinformatics Specialization.

They illustrate foundational computational methods for biological sequence analysis.

The notebooks are educational implementations. Additional code review, testing, and documentation improvements may be introduced as the repository is progressively refined.

## Related Repository

[Bioinformatics — Coursera Series](../)

For genome assembly and peptide sequencing algorithms, see [Bioinformatics II](../Bioinformatics_II/).
