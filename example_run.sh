snakemake  --profile simple  --configfile configs/test.yaml -f  calculate_openmm_ff1
snakemake  --profile simple  --configfile configs/test.yaml -k  test_run/pdb_proc/process/3BZD_seqresMatched_promod.pdb
#nice -n 19 snakemake --use-conda --conda-frontend mamba --configfile configs/test.yaml -c 80 --use-singularity -f  rosetta_static_evaluation
