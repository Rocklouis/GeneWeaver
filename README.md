
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
