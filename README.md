## Workflow setup

1. Install `conda`:
```bash
   wget -P miniconda https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh &&
   chmod 755 ./miniconda/Miniconda3-latest-Linux-x86_64.sh &&
   ./miniconda/Miniconda3-latest-Linux-x86_64.sh
```

2. Add path of `miniconda` to `.bashrc` if not selected option to add automatically during installation:
```bash
   cd ~/ && pth=$(pwd) &&
   echo "PATH=$PATH:$pth/miniconda3/bin" >> ~/.bashrc
```

3. Install mamba frontend for conda. This helps to handle dependancies installation.
```bash
    conda install -c conda-forge mamba
```

4. Clone the repository. nd enter it. Note - it contains submodules.
```bash 
    git clone --recurse-submodules  git@github.com:GediminasA/binding_evaluator.git
    cd binding_evaluator
```
5. Create your `conda` environment:
 ```bash
    mamba env create -f envs/binding_evaluator.yaml 
 ```

6. Activate created environment:
```bash
    conda activate binding_evaluator
```


To run prodigy evaluations on a toy set of complexes you could run this command:
```
snakemake --configfile configs/test.yaml  --use-conda --conda-frontend mamba --use-singularity -c 6  -f prodigy_static_evaluation
```
To run a free sasa calculation on a toy set of complexes you coud run this command:
```
snakemake --configfile configs/test.yaml  --use-conda --conda-frontend mamba --use-singularity -c 6  -f get_freesasa
```