
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
