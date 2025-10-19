# RNA-seq Analysis Pipeline

![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.10.5-brightgreen.svg)
![Docker](https://img.shields.io/badge/docker-enabled-blue.svg)
![Conda](https://img.shields.io/badge/conda-enabled-green.svg)

## Introduction

**RNA-seq Analysis Pipeline** is a comprehensive bioinformatics pipeline for RNA sequencing data analysis. The pipeline processes paired-end or single-end Illumina RNA-seq data through quality control, trimming, alignment, quantification, and generates comprehensive reports. It supports both alignment-based and pseudo-alignment approaches for transcript quantification, making it suitable for various RNA-seq analysis needs.

The pipeline integrates multiple state-of-the-art tools and provides flexible options for different analysis strategies, from basic gene-level counting to advanced transcript-level quantification with isoform detection.

## Pipeline Overview

The pipeline performs the following main steps:

### Quality Control & Pre-processing
1. **Raw Read QC** ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) - Quality assessment of raw sequencing reads
2. **Adapter Trimming** ([`TrimGalore`](https://github.com/FelixKrueger/TrimGalore) or [`SEQTK`](https://github.com/lh3/seqtk)) - Remove adapters and low-quality sequences
3. **Post-trim QC** ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) - Quality assessment after trimming

### Reference Preparation
4. **Transcriptome Generation** ([`GFFREAD`](https://github.com/gpertea/gffread)) - Extract transcript sequences from genome and GTF
5. **Index Building** - Build indices for selected aligners and pseudo-aligners

### Alignment & Quantification
6. **RNA-seq Alignment** (Choose one):
   - [`HISAT2`](https://github.com/DaehwanKimLab/hisat2) - Fast splice-aware aligner (default)
   - [`STAR`](https://github.com/alexdobin/STAR) - Ultra-fast universal RNA-seq aligner

7. **Post-alignment Processing**:
   - [`SAMtools`](https://github.com/samtools/samtools) - BAM sorting and indexing
   - [`Picard MarkDuplicates`](https://broadinstitute.github.io/picard/) - Mark duplicate reads

8. **Quantification** (Multiple approaches):
   - **Gene-level**: [`featureCounts`](https://github.com/ShiLab-Bioinformatics/subread) - Count reads per gene
   - **Transcript-level**: [`SALMON`](https://github.com/COMBINE-lab/salmon) - Fast pseudo-alignment and quantification

### Reporting
9. **Comprehensive Report** ([`MultiQC`](http://multiqc.info/)) - Aggregate all QC metrics and results

## Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/docs/latest/install) (≥24.10.5)
- [Docker](https://docs.docker.com/engine/installation/) or [Conda](https://conda.io/miniconda.html)

### Test Run

Test the pipeline with provided test data:

```bash
nextflow run . -profile test,docker --input ./assets/samplesheet.csv --outdir test_results
```

### Full Analysis

#### 1. Prepare Samplesheet

Create a CSV samplesheet with your input data:

`samplesheet.csv`:
```csv
sample,fastq_1,fastq_2,strandedness,library,seq_center
SRR6357070_2,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz,reverse,illumina,N.A.
SRR6357071_2,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357071_1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357071_2.fastq.gz,reverse,illumina,N.A.
SRR6357072_2,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357072_1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357072_2.fastq.gz,reverse,illumina,N.A.
SRR6357076_1,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004/SRR6357076_2.fastq.gz,reverse,illumina,N.A.
```

**Column descriptions:**
- `sample`: Unique sample identifier
- `fastq_1`: Path to first read file (R1)
- `fastq_2`: Path to second read file (R2, for paired-end)
- `strandedness`: Library strandedness (`forward`, `reverse`, or `unstranded`)
- `library`: Sequencing library type (e.g., `illumina`)
- `seq_center`: Sequencing center


#### 2. Run the Pipeline

**Basic usage:**
```bash
nextflow run . \
    -profile docker,test \
    --input samplesheet.csv \
    --outdir results \
    --genome sacCer3
```

**Advanced usage with custom parameters:**
```bash
nextflow run . \
    -profile docker,arm \
    --input samplesheet.csv \
    --outdir results \
    --genome sacCer3 \
    --mode hisat2 \
    --mark_duplicates true \
    --trimmer trimgalore \
    --run_fastqc_at_start true \
    --run_fastqc_after_trim true
```

## Configuration Options

### Execution Profiles

- `docker` - Run with Docker containers
- `conda` - Run with Conda environments  
- `singularity` - Run with Singularity containers
- `arm` - For Apple Silicon Macs (use with docker: `-profile docker,arm`)
- `test` - Run with test data

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | - | Path to samplesheet CSV |
| `--outdir` | - | Output directory |
| `--genome` | - | Reference genome (e.g., `sacCer3`, `hg38`) |
| `--mode` | `hisat2` | RNA-seq mode (`salmon`, `hisat2` or `star`) |
| `--trimmer` | `trimgalore` | Read trimmer (`trimgalore` or `seqtk`) |
| `--mark_duplicates` | `true` | Mark duplicates|
| `--run_fastqc_at_start` | `true` | QC on raw reads |
| `--run_fastqc_after_trim` | `true` | QC on trimmed reads |

## Pipeline Output

The pipeline generates the following outputs in the specified `--outdir`:

### Main Results
```
results/
├── fastqc/                     # Raw read quality control
├── trimgalore/                 # Adapter trimming logs and stats
├── hisat2/ (or star/)          # Alignment results and logs
├── samtools/                   # Sorted BAM files and indices
├── picard/                     # Duplicate marking metrics
├── subread/                    # Gene-level count matrices
├── salmon/                     # Transcript-level quantification
│   ├── *.quant.sf             # Transcript abundance estimates
│   └── aux_info/              # Auxiliary quantification info
├── multiqc/                    # Comprehensive QC report
│   └── multiqc_report.html    # Main report (open this first!)
└── pipeline_info/             # Execution reports and metadata
```

### Key Output Files

| File | Description |
|------|-------------|
| `multiqc/multiqc_report.html` | **Main QC report** - comprehensive overview |
| `subread/combined_counts.txt` | Gene-level count matrix (all samples) |
| `salmon/*/quant.sf` | Transcript abundance per sample for mode salmon |
| `subread/id.featureCounts.tsv` | For mode hisat2 / star|
| `pipeline_info/execution_report.html` | Pipeline execution metrics |

### Quality Control Metrics

The MultiQC report includes:
- **Read quality** (FastQC): Base quality, GC content, adapter content
- **Trimming stats** (TrimGalore): Reads trimmed, adapter removal
- **Alignment stats** (HISAT2/STAR): Mapping rates, splice junctions
- **Duplicate rates** (Picard): PCR duplicate percentages  
- **Gene counting** (featureCounts): Assignment rates, feature types
- **Transcript quantification** (SALMON): Library type, fragment length

## Troubleshooting

### Common Issues

**SALMON quantification fails with 0 assigned fragments:**
- Ensure GTF and genome FASTA are from the same reference build
- Check that read strandedness is correctly specified
- Verify input FASTQ files are not corrupted

**MultiQC shows incorrect aligner (e.g., Bowtie2 instead of HISAT2):**
- This is a known MultiQC detection issue - the actual analysis uses the correct aligner
- Check `pipeline_info/` for accurate tool versions and execution details

**Memory or disk space errors:**
- Adjust resource limits in `conf/base.config`
- Use smaller test datasets first
- Monitor system resources during execution

**Test profile can not be executed due (path)/assets/samplesheet does not match regex:**
- Ensure that the project path does not include whitespaces
- Use parameter --input ./assets/samplesheet.csv  in addition

## Version History

### v1.0.0dev (Current)
- Initial implementation with dual quantification approach
- Support for HISAT2 and STAR aligners
- SALMON pseudo-alignment integration
- Comprehensive MultiQC reporting
- Docker and Conda execution support

## Development & Contributions

### Repository Structure
```
├── main.nf                    # Main pipeline script
├── workflows/rnasequencing.nf # Core workflow logic
├── modules/nf-core/          # Individual tool modules
├── conf/                     # Configuration files
├── assets/                   # Pipeline assets and examples
└── misc/                     # Custom configurations and test data
```

### Contributing
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Credits & Acknowledgments

**RNA-seq Analysis Pipeline** was originally written by **Mykyta Borodin** and **Bastian Eisenmann** as part of the Computational Workflows course.

### Built With
- [Nextflow](https://nextflow.io/) - Workflow management system
- [nf-core](https://nf-co.re/) - Framework and modules
- [Docker](https://docker.com/) - Containerization platform

### Tool Credits
We acknowledge the developers of all integrated tools:
- FastQC, TrimGalore, HISAT2, STAR, SALMON, SAMtools, Picard, featureCounts, MultiQC

## Citations

## Citations

If you use this RNA-seq Analysis Pipeline for your research, please cite the following:

### Pipeline
- **RNA-seq Analysis Pipeline** - Borodin, M. & Eisenmann, B. (2025). Computational Workflows Course Project.

### Core Tools
- **Nextflow**: Di Tommaso, P. et al. Nextflow enables reproducible computational workflows. _Nat Biotechnol_ 35, 316–319 (2017). [doi:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)
- **nf-core**: Ewels, P. et al. The nf-core framework for community-curated bioinformatics pipelines. _Nat Biotechnol_ 38, 276–278 (2020). [doi:10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)

### Analysis Tools
- **FastQC**: Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.
- **TrimGalore**: Krueger, F. (2015). Trim Galore: a wrapper tool around Cutadapt and FastQC.
- **HISAT2**: Kim, D. et al. HISAT: a fast spliced aligner with low memory requirements. _Nat Methods_ 12, 357–360 (2015).
- **STAR**: Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. _Bioinformatics_ 29, 15–21 (2013).
- **SALMON**: Patro, R. et al. Salmon provides fast and bias-aware quantification of transcript expression. _Nat Methods_ 14, 417–419 (2017).
- **SAMtools**: Li, H. et al. The Sequence Alignment/Map format and SAMtools. _Bioinformatics_ 25, 2078–2079 (2009).
- **Picard**: Broad Institute. Picard Toolkit. http://broadinstitute.github.io/picard/
- **featureCounts**: Liao, Y. et al. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. _Bioinformatics_ 30, 923–930 (2014).
- **MultiQC**: Ewels, P. et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. _Bioinformatics_ 32, 3047–3048 (2016).

An extensive list of references for all tools can be found in the [`CITATIONS.md`](CITATIONS.md) file.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**Happy RNA-seq analysis! 🧬📊**
