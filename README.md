
# GeneWeaver
GPU-Accelerated CRISPR Alignment Engine

## Week 1 - Day 1

Implemented:
- Pure Python CRISPR off-target alignment baseline
- Forward and reverse strand scanning
- PAM sequence validation
- Mismatch counting
- BioPython FASTA parsing
- DNA normalization
- Genome chunking
- CPU benchmark
- Automated tests

## Installation

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

## Run tests

```powershell
$env:PYTHONPATH = "$PWD\src"
python -m pytest -v
```

## Run CPU benchmark

```powershell
$env:PYTHONPATH = "$PWD\src"
python benchmarks/benchmark_cpu.py --bases 1000000
```

## Current milestone

Week 1 - Day 1: CPU Baseline and FASTA Data Pipeline.

## Week 1 - Day 2

Implemented:
- Multi-record FASTA processing
- DNA sequence normalization
- DNA validation
- Numerical DNA encoding using NumPy uint8
- DNA decoding
- Genome chunk coordinate tracking
- Automatic alignment overlap calculation
- Alignment-safe chunk preparation
- DNA encoding benchmark
- Automated data pipeline tests

DNA encoding:
- A = 0
- C = 1
- G = 2
- T = 3
- N = 4

Status: Implementation complete; awaiting local test and benchmark verification.


## Day 3: Streaming and Large-Scale Genome Processing

### Synthetic datasets
- 100 KB synthetic FASTA dataset
- 1 MB synthetic FASTA dataset
- 10 MB synthetic FASTA dataset
- Reproducible generation using fixed random seeds

### Streaming FASTA implementation
- Line-by-line FASTA reading
- Chunk-based sequence processing
- Configurable overlap
- Genome coordinate preservation
- Multi-record FASTA support
- Input validation

### Streaming benchmark (10 MB)

| Metric | BioPython parser | Streaming reader |
|---|---:|---:|
| Bases processed | 10,000,000 | 10,000,000 |
| Chunks | 10 | 10 |
| Execution time | 1.1633 s | 6.3681 s |
| Peak traced RAM | 113.71 MB | 3.00 MB |

The streaming implementation reduced peak traced Python memory
by approximately 97.4%, but took approximately 5.47 times longer
in this benchmark.

Correctness check: PASS.


### Streaming optimization results

The streaming reader was optimized by replacing repeated string
concatenation with list-based sequence buffering.

Optimized benchmark (10 MB synthetic FASTA):

| Metric | BioPython parser | Optimized streaming |
|---|---:|---:|
| Bases processed | 10,000,000 | 10,000,000 |
| Chunks | 10 | 10 |
| Execution time | 1.1139 s | 2.1577 s |
| Peak traced RAM | 113.71 MB | 4.55 MB |

Optimization comparison:
- Original streaming time: 6.3681 seconds
- Optimized streaming time: 2.1577 seconds
- Approximately 2.95x faster than original streaming
- Approximately 66.1% lower execution time
- Approximately 96.0% lower peak traced memory than BioPython

Correctness: PASS
Full regression suite: 34 tests passed.

Note: Benchmark timings are from individual runs and may vary
with system load. Peak traced RAM is not total process RAM.

### CPU alignment optimization

Optimized the standard NGG PAM scanning path using
Python's optimized string-search operation.

Additional improvements:
- Early termination during mismatch counting
- PAM-first candidate filtering
- Separate generic PAM scanning path
- Preserved reverse-strand scanning
- Preserved genomic coordinate calculations

#### Benchmark results

| Genome size | Original | Optimized | Speedup |
|---|---:|---:|---:|
| 100,000 bases | 0.150006 s | 0.014916 s | 10.06x |
| 1,000,000 bases | 1.967636 s | 0.169692 s | 11.60x |

Candidate hit counts matched:
- 100 KB: 3 hits
- 1 MB: 1 hit

Alignment tests: 6 passed.
Full regression suite: 34 tests passed.

Note: These are individual benchmark runs. Results may vary
with hardware and system load. This is a computational prototype,
not a clinically validated CRISPR prediction tool.
