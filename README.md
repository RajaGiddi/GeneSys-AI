# GeneSys AI

A conversational interface to a deterministic bioinformatics toolkit.

GeneSys lets biologists, clinicians, and students run genomic analyses in natural
language without writing code. The design point is that the language model never
performs the analysis. It selects and parameterizes a function; the function runs
real bioinformatics code and returns its output. Every result is a verifiable
operation, not generated text.

The routing layer is hand-written. GeneSys was built before agent frameworks and
before the Model Context Protocol existed, and it solves the same problem MCP now
standardizes: giving a language model reliable, typed access to domain tools.

## Architecture

- **Dispatch layer** — parses intent from conversation, selects a function from the
  DNA Tool Kit, extracts and validates arguments, and returns structured results.
- **DNA Tool Kit** — the analysis functions below. Each is independently callable
  and testable without the model in the loop.
- **Ingestion** — FASTA/FASTQ upload for DNA and protein sequences; tabular data
  (CSV, PDF, dataframes) for general analysis.

## Implemented

**Sequence analysis**
- Complementary sequence computation
- mRNA transcription; mRNA-to-protein translation
- Nucleotide counting; GC content, including per-window GC for a specified *k*
- Hamming distance between equal-length sequences
- Open reading frame (ORF) identification with translation to protein

**Structural and comparative**
- Multiple sequence alignment
- Phylogenetic tree generation from MSA
- SNP detection
- Motif finding
- Restriction site detection, with suggested restriction enzymes
- Protein mass calculation from amino acid composition
- 3D protein structure visualization

**Tabular data**
- Dataframe creation, editing, merging, sorting, filtering, aggregation
- Type conversion and random sampling
- Descriptive statistics and hypothesis testing
- Interactive charts and plots

## Roadmap

Not yet implemented: BLAST (BLASTn/p/x, tBLASTn, tBLASTx, BLAST2Seq), primer design
automation, Gene Ontology enrichment, sequence editing, pivot tables, data
transposition, export/import across formats.

Longer term: variant pathogenicity prediction, CRISPR-Cas9 off-target
identification, variant functional annotation, drug-target interaction prediction,
allele frequency computation, molecular docking, protein-protein interaction
prediction, RNA-Seq differential expression, and more storage-efficient nucleotide
data structures.
