### This workflow given a set of mutations evaluates their impact on a complex

## Funding
This project has received funding from European Regional Development Fund (project No 13.1.1-LMT-K-718-05-0023) under grant agreement with the Research Council of Lithuania (LMTLT). Funded as European Union's measure in response to Cov-19 pandemic."

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

7. Checkout exemplary data git lfs:
```bash
    git lfs fetch
    git lfs install
    git lfs checkout
```

Make sure that uidmap is installed in the system, this coud ne done like:
 ```bash
 sudo apt-get install uidmap
 ```

The target mutations can by generated as:
```
snakemake --profile local --configfile configs/antibody.yaml  -k antibody_run/mutants_structure_generation/TEMPLATES/todoList
```

A workflow assesing muttions impact on protein binding can be invoked like this 

```
snakemake --profile local --configfile configs/antibody.yaml  -k get_summary_of_binding
```
