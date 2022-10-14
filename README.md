
To run prodigy evaluations on a toy set of complexes you coud run this command:
```
snakemake --configfile configs/test.yaml  --use-conda --conda-frontend mamba --use-singularity -c 6  -f prodigy_static_evaluation
```
To run a free sasa calculation on a toy set of complexes you coud run this command:
```
snakemake --configfile configs/test.yaml  --use-conda --conda-frontend mamba --use-singularity -c 6  -f get_freesasa
```