git lfs pull
snakemake  --profile simple  --configfile configs/prodigy.yaml -f t4t   -c 1 
snakemake  --profile simple  --configfile configs/prodigy.yaml -f  prodigy_run/pdb_proc/process/1WEJ_seqresMatched.pdb -c 1
snakemake  --profile local  --configfile configs/prodigy.yaml -f  prodigy_test
snakemake  --profile simple  --configfile configs/test.yaml -f  calculate_openmm_ff1
snakemake  --profile simple  --configfile configs/test.yaml -k  test_run/pdb_proc/process/3BZD_seqresMatched_promod.pdb
snakemake -c 24 --use-conda --configfile configs/prodigy.yaml    --use-conda --use-singularity -f prodigy_run/pdb_proc/process/4CPA_seqresMatched_woHOH_promod.pdb
#nice -n 19 snakemake --use-conda --conda-frontend mamba --configfile configs/test.yaml -c 80 --use-singularity -f  rosetta_static_evaluation
