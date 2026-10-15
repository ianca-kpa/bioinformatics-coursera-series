
# Bioinformatics — Coursera Series

A collection of Python implementations and Jupyter notebooks developed during the Bioinformatics Specialization offered by the University of California San Diego on Coursera.

The repository documents my learning and practice in bioinformatics algorithms, DNA sequence analysis, genome assembly, and computational problem-solving.

Several notebooks have been subsequently reviewed to improve code readability, input validation, automated testing, and documentation.

## Courses

### Bioinformatics I — Finding Hidden Messages in DNA

Explore algorithms for identifying patterns and motifs in DNA sequences.

Topics include:

- Pattern matching and approximate matching
- Hamming distance and reverse complements
- K-mer frequency analysis and clump finding
- Motif search algorithms
- Profile matrices and consensus sequences
- Greedy, randomized, and Gibbs sampling approaches

[Explore Bioinformatics I](./Bioinformatics_I/)

### Bioinformatics II — Genome Sequencing

Explore computational approaches to genome reconstruction and peptide sequencing.

Topics include:

- K-mer composition and overlap graphs
- De Bruijn graph construction
- Eulerian paths and cycles
- String reconstruction from k-mers
- Maximal non-branching paths and contig generation
- Peptide sequencing and mass spectrum analysis
- DNA translation and peptide encoding

[Explore Bioinformatics II](./Bioinformatics_II/)

## Selected Implementations

The following notebooks have undergone additional review and validation:

| Algorithm | Main concepts |
|---|---|
| [De Bruijn Graph Construction](./Bioinformatics_II/01-graph-based-assembly/DeBruijn_Patterns_Bioinfo_II.ipynb) | Graph construction, k-mer overlaps, repeated edges |
| [Eulerian Path](./Bioinformatics_II/01-graph-based-assembly/Eulerian_Path_Bioinfo_II.ipynb) | Graph traversal, degree constraints, connectivity validation |
| [String Reconstruction](./Bioinformatics_II/01-graph-based-assembly/StringReconstruction_Bioinfo_II.ipynb) | Sequence assembly, Eulerian traversal, k-mer conservation |
| [Contig Generation](./Bioinformatics_II/01-graph-based-assembly/Find_contigs_Bioinfo_II.ipynb) | Non-branching paths, isolated cycles, contig reconstruction |
| [Maximal Non-Branching Paths](./Bioinformatics_II/01-graph-based-assembly/MaximalNonBranchingPaths_Bioinfo_II.ipynb) | Graph topology, edge multiplicity, isolated-cycle detection |

## Technical Approach

The notebooks demonstrate practical work with:

- Python programming and algorithm implementation
- Graph-based representations of biological sequences
- Sequence reconstruction and k-mer analysis
- Input validation and exception handling
- Automated tests using assertions
- Edge and k-mer frequency verification
- Reproducible computational exercises

## Tools

- Python 3
- Jupyter Notebook
- Python standard library, including `collections`

Additional dependencies, where applicable, are defined by individual notebooks.

## Educational Context

These implementations originated from programming exercises completed during the UC San Diego Bioinformatics Specialization.

The repository also documents subsequent improvements made to selected notebooks, including additional tests, documentation, and robustness checks.

Where external reference implementations are retained for educational comparison, they are identified separately from my own work.

## Author

**Ianca Karine Prudencio de Albuquerque**

- [LinkedIn](https://linkedin.com/in/ianca-kpa)
- [ORCID](https://orcid.org/0000-0003-2495-4821)