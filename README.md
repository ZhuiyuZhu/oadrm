# OADRM: Organ-Aging Drug Repurposing Matrix

## Overview
A computational framework for organ-specific anti-aging drug screening, 
bridging UK Biobank organ aging clocks with LINCS L1000 transcriptomic perturbations.

## Pipeline
1. **Organ Signature Construction**: Extract protein aging signatures → map to Gene Symbols → GTEx v8 validation
2. **LINCS Screening**: Modified Connectivity Score (mCS) across 8 organs × 21,304 compounds × 83 cell lines
3. **Tissue-Matching Weighting**: Homologous vs. heterologous cell-line stratification
4. **Candidate Prioritization**: Four-dimensional scoring (efficacy, specificity, robustness, safety)
5. **Real-World Validation**: UK Biobank prescription records → organ function trajectories

## Key Results
- **8 organs** with validated transcriptional signatures
- **139 organ-specific candidates** after cytotoxicity filtering
- **Zero broad-spectrum compounds** passed specificity thresholds
- Top candidates: hyperoside (lung), quinapril (brain), methylene-blue (blood)

## Requirements
- Python 3.13+
- pandas, numpy, scipy, h5py, matplotlib

## Citation
Zhu et al., "Organ-Aging Drug Repurposing Matrix", [Journal], [Year]
