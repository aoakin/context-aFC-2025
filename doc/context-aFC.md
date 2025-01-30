Continuous modeling of genetic regulatory variation across cellular contexts
---
**Authors:**\
Ayomikun O. Akinrinade\
Pejman Mohammadi\
*Winter 2025*

## Jan26-2025
Task: PCA on median TPM

- Top two PC axes make up ~75% of the variation in expression
- Entirely differentiated by whether or not a tissue is blood or not
- PC1: Top 10 eigengenes: 9 mitochondrial genes [-]  and hemoglobin subchain beta (HBB) [+]
- PC2: HBB and HBA2,  S100A9 (inflammation + immune response), 7 mitochondrial genes [all, +]

### Experiment 1 - remove blood and pancreas; all genes
- 73% of expr. var. (PC1) primarily epxlained by mitochondrial gene expr.

### Experiment 2 - all tissues; remove mitochondrial genes
- PC1 (blood vs not blood), PC2 (pancreas vs not pancreas)

### Experiment 3 - remove blood and pancreas; remove mitochondrial gene
- High hit genes look slightly more interesting biologically?? (***TODO*** - check back on this)

**PC1** (in order of eigenvector magnitude):\
PRL/+, GH1/+, POMC/+, CGA/+, ACTB/-, LHB/+, NNAT/+, S100A9, FTL/-, DES/-

**PC2**:\
S100A9/-, KRT13/-, KRT4/-, S100A8/-, SPRR3/-, KRTS/-, KRT6A/-, CRNN/-, RHCG/-, MTATP6P1/+


## Jan27-2025
- 