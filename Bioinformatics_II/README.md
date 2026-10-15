
# Bioinformatics II — Genome Sequencing

Python implementations and Jupyter notebooks developed during the *Genome Sequencing* course of the Bioinformatics Specialization offered by the University of California San Diego on Coursera.

This collection explores computational approaches to genome assembly, peptide sequencing, and biological sequence analysis.

Several notebooks have been subsequently reviewed to improve code readability, input validation, automated testing, and documentation.

## Repository Structure

### 01 — Graph-Based Genome Assembly

[Explore notebooks](./01-graph-based-assembly/)

Algorithms for reconstructing DNA sequences using graph representations of overlapping k-mers.

Topics include:

- K-mer composition and overlap graphs
- De Bruijn graph construction
- Eulerian paths and cycles
- String reconstruction from k-mers
- Maximal non-branching paths
- Contig generation

Selected notebooks have been reviewed to improve graph traversal, edge conservation, input validation, and automated testing.

### 02 — Peptide Sequencing

[Explore notebooks](./02-peptide-sequencing/)

Computational approaches to peptide sequencing using theoretical mass spectra.

Topics include:

- Linear and cyclic peptide spectra
- Cyclopeptide sequencing
- Spectrum scoring
- Leaderboard algorithms
- Spectral convolution

### 03 — Translation and Coding

[Explore notebooks](./03-translation-coding/)

Exercises connecting nucleotide sequences to peptide sequences.

Topics include:

- DNA and RNA sequence translation
- Peptide encoding
- Reading frames and sequence interpretation

### Additional Exercises

[Explore exercises](./exercises/)

Supplementary exercises and computational practice.

## Selected Graph-Based Assembly Implementations

The following notebooks have undergone additional review and validation.

| Notebook | Main concepts |
|---|---|
| [De Bruijn Graph Construction](./01-graph-based-assembly/DeBruijn_Patterns_Bioinfo_II.ipynb) | Graph construction, k-mer overlaps, repeated edges |
| [Eulerian Path](./01-graph-based-assembly/Eulerian_Path_Bioinfo_II.ipynb) | Graph traversal, degree constraints, connectivity validation |
| [String Reconstruction](./01-graph-based-assembly/StringReconstruction_Bioinfo_II.ipynb) | Sequence assembly, Eulerian traversal, k-mer conservation |
| [Contig Generation](./01-graph-based-assembly/Find_contigs_Bioinfo_II.ipynb) | Non-branching paths, isolated cycles, contig reconstruction |
| [Maximal Non-Branching Paths](./01-graph-based-assembly/MaximalNonBranchingPaths_Bioinfo_II.ipynb) | Graph topology, repeated edges, isolated-cycle detection |

## Implementation and Validation

The notebooks contain educational Python implementations developed during the course.

Selected graph-based assembly notebooks have subsequently been improved with:

- Function documentation and clearer code structure
- Input validation and error handling
- Automated tests using Python assertions
- Verification of graph traversal and edge multiplicity
- Tests for branching paths, isolated cycles, and repeated edges
- K-mer and edge-frequency conservation checks

These improvements are documented within the corresponding notebooks.

## Tools and Environment

- Python 3
- Jupyter Notebook or JupyterLab
- Visual Studio Code with the Jupyter extension (optional)
- Python standard library, including `collections`

Some notebooks may require additional packages or access to their associated datasets.

## Running the Notebooks

The notebooks can be opened using Jupyter Notebook, JupyterLab, or the Jupyter extension in Visual Studio Code.

### Requirements

- Python 3
- Jupyter Notebook or JupyterLab
- Additional packages when required by individual notebooks

### Getting Started

Clone the repository:

```bash
git clone https://github.com/ianca-kpa/bioinformatics-coursera-series.git
cd bioinformatics-coursera-series
```

Open the repository in VS Code or launch Jupyter from your preferred Python environment.

For notebooks in `01-graph-based-assembly/`, some examples read files from the local `Dataset/` directory.

These relative paths assume that the notebook's working directory is:

`Bioinformatics_II/01-graph-based-assembly/`

If the notebook is executed from another directory, adjust the dataset path accordingly.

### Dataset Organization

Input datasets used by selected graph-based assembly notebooks are stored in:

`Bioinformatics_II/01-graph-based-assembly/Dataset/`

For example, the contig generation notebook reads:

`Dataset/dataset_30189_5.txt`

This relative path works when the notebook's working directory is `01-graph-based-assembly/`.

### Validation

Selected notebooks include assertion-based tests covering algorithm correctness, input validation, and conservation of graph edges or k-mer frequencies.

Run the notebook cells in order to reproduce the examples and validation results.

The tests use small, representative datasets and are intended to demonstrate expected algorithmic behavior. They do not constitute exhaustive validation of every possible graph structure.

## Educational Context

These implementations originated from coursework and personal programming practice during the UC San Diego Bioinformatics Specialization.

Subsequent revisions focus on code clarity, correctness, testing, and reproducibility.

Any external reference implementations retained for educational comparison are identified separately from the author's own work.

The notebooks are educational implementations rather than production-ready genome assembly pipelines.

## Related Repository

[Bioinformatics — Coursera Series](../)
