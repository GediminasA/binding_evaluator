### PROCESS CSV file with mutations ###
import pandas as pd 
df_muts = pd.read_csv(config["mutants"])
pdb_template_stems = list(set(df_muts.PDB))
seqs_4mutations = list(set([s.strip() for s in df_muts.Template]))
pair4mut = ["=".join(l) for l in  list(set(list(zip([s.strip() for s in df_muts.PDB],[s.strip() for s in df_muts.Template]))))]
mutrez = work_dir + "/rezults"

rule mutants_targets:
    input:
        pdb_templates = expand(
            work_dir+"/processed/{stem}.pdb",
            stem=pdb_template_stems),
        pdb_sequences = expand(
            work_dir+"/processed_info/{stem}.fasta",
            stem=pdb_template_stems),
        pdb_sequences_aligned_to_templates = expand(
            work_dir+"/mutants_sequence_generation/{stem}_vs_templates_blastp.txt",
            stem=pdb_template_stems),
        pdb_templates_map = expand(
            work_dir+"/mutants_sequence_generation/{pair}.csv", #map input seq to pdb map 
            pair = pair4mut),
        pdb_tempplate_mutation_data = expand(
            work_dir+"/mutants_sequence_generation/{pair}=template_mutation_data__EVOEF.csv",
            pair = pair4mut),
        ddg_results = mutrez + "/mutations_data_4EVOEFmodels_results.csv"
        

rule collect_data_befor_modelling:
    input:
        cleaned_haps = expand(
            work_dir+"/mutants_sequence_generation/{pair}=template_mutation_data__cleanhap.csv",
            pair = pair4mut),
        evoef_data = expand(
            work_dir+"/mutants_sequence_generation/{pair}=template_mutation_data__EVOEF.csv",
            pair = pair4mut)
    output:
        cleaned_haps = mutrez + "/mutations_data_cleanHap.csv",
        evoef_data = mutrez + "/mutations_data_4EVOEFmodels.csv"
    shell:
        """
            cat {input.cleaned_haps} > {output.cleaned_haps}
            cat {input.evoef_data} > {output.evoef_data}
        """


#### rules to prepare for EVOEF modelling

rule align_to_mutations_template:
    input:
        template = config["mutants_templates"],
        fasta = work_dir+"/processed_info/{stem}.fasta"
    output:
        coordinates = work_dir+"/mutants_sequence_generation/{stem}_vs_templates_blastp.txt"
    shell:
        """
        blastp -query {input.fasta} -subject {input.template} -outfmt '6 qseqid sseqid qstart qend sstart send length pident qseq sseq'  -subject_besthit &>  {output}
        """

rule get_template_pdb_map:
    input:
        coordinates = work_dir+"/mutants_sequence_generation/{pdb}_vs_templates_blastp.txt"
    output:
        work_dir+"/mutants_sequence_generation/{pdb}={seqtempl,[^=]+}.csv",
    script:
     "notebooks/get_template_pdb_map.py.py"
    # notebook:
    #   "notebooks/get_template_pdb_map.py.ipynb"

rule get_template_mutation_data:
    input:
        map_s_p = work_dir+"/mutants_sequence_generation/{pdb}={seqtempl}.csv",
        data = config["mutants"]
    output:
        work_dir+"/mutants_sequence_generation/{pdb,[^=]+}={seqtempl,[^=]+}=template_mutation_data__EVOEF.csv",
        work_dir+"/mutants_sequence_generation/{pdb,[^=]+}={seqtempl,[^=]+}=template_mutation_data__cleanhap.csv"
    # notebook:
    #    "notebooks/get_template_mutation_data.r.ipynb"   
    script:
        "notebooks/get_template_mutation_data.r.R"   

#### rules to run EVOEF modelling

#use checkpoint is used as a prtiori number ofstructure that must be modelled is not known, mostly due to variable sequence length available in PDBs
checkpoint create_idividual_tasks_4_modeling_with_EVOEF:
    input:
        mutrez + "/mutations_data_4EVOEFmodels.csv",
    output:
        directory(work_dir + "/mutants_structure_generation/EVOEF/todoList") # writes out small text file for each mutants to be generated
    notebook:
        "notebooks/create_idividual_tasks_4_modeling_with_EVOEF.py.ipynb"


def aggregate_EVOEF_results(wildcards):
    checkpoint_output = checkpoints.create_idividual_tasks_4_modeling_with_EVOEF.get(**wildcards).output[0]
    stems = glob_wildcards(os.path.join(checkpoint_output, "{i,[^\.].+}")).i #rehex prevents snakemake temps to be included
    print(stems)
    return(expand(work_dir + "/mutants_structure_scoring/EVOEF/scores/{stem}.ddg",stem=stems))
    #return(expand(work_dir + "/mutants_structure_scoring/EVOEF/scores/{stem}.sc",stem=stems))

rule run_evoEF1_modelling:
    input:
        container = "containers/evoef1.sif", # to make sure it is built
        pdb = work_dir+"/processed/{pdb}.pdb",
        mutant_file = work_dir + "/mutants_structure_generation/EVOEF/todoList/{pdb}={chain}={mutations}"
    params:
        pdb = os.path.abspath(work_dir+"/processed/{pdb}.pdb"),
        mutant_file = os.path.abspath(work_dir + "/mutants_structure_generation/EVOEF/todoList/{pdb}={chain}={mutations}"),
        wdir = work_dir + "/mutants_structure_generation/EVOEF/structures/{pdb}={chain}={mutations}",
        num_of_runs = 10,
        outmut = work_dir + "/mutants_structure_generation/EVOEF/structures/{pdb}={chain}={mutations}/{pdb}_Model_0001.pdb",
        outwt = work_dir + "/mutants_structure_generation/EVOEF/structures/{pdb}={chain}={mutations}/{pdb}_Model_0001_WT.pdb"

    output:
        mut = work_dir + "/mutants_structure_generation/EVOEF/structures/{pdb}={chain}={mutations}.pdb",
        wt = work_dir + "/mutants_structure_generation/EVOEF/structures/{pdb}={chain}={mutations}=WT.pdb"
    log:
        os.path.abspath(work_dir + "/mutants_structure_generation/EVOEF/structures/{pdb}={chain}={mutations}.log")
    singularity:
        "containers/evoef1.sif"
    shell:
        """
            startD=`pwd`
            echo {params.pdb}
            mkdir -p {params.wdir}
            cd {params.wdir}
            /EvoEF-master/EvoEF  --command=BuildMutant   --pdb {params.pdb} --num_of_runs {params.num_of_runs} --mutant_file {params.mutant_file} > {log} 
            cd $startD
            mv {params.outmut} {output.mut}
            mv {params.outwt} {output.wt}

        """


rule run_evoEF1_eval:
    input:
        mut = work_dir + "/mutants_structure_generation/EVOEF/structures/{pdb}={chain}={mutations}.pdb",
        wt = work_dir + "/mutants_structure_generation/EVOEF/structures/{pdb}={chain}={mutations}=WT.pdb",
        interacting_groups = work_dir+"/processed_info/{pdb,[^_]+}_interactigGroups.tsv"
    output:
        mut = work_dir + "/mutants_structure_scoring/EVOEF/scores/{pdb}={chain}={mutations}.sc",
        wt = work_dir + "/mutants_structure_scoring/EVOEF/scores/{pdb}={chain}={mutations}=WT.sc",
    params:
        mut = os.path.abspath(work_dir + "/mutants_structure_generation/EVOEF/structures/{pdb}={chain}={mutations}.pdb"),
        wt = os.path.abspath(work_dir + "/mutants_structure_generation/EVOEF/structures/{pdb}={chain}={mutations}=WT.pdb"),
        wdir = work_dir + "/mutants_structure_scoring/EVOEF/scores/{pdb}={chain}={mutations}",
        num_of_runs = 10,
        part1_part2 = get_interacting_chains4evoEF1
    singularity:
        "containers/evoef1.sif"
    shell:
        """
            /EvoEF-master/EvoEF --command=ComputeBinding --pdb  {params.mut} --split  {params.part1_part2} > {output.mut}  
            /EvoEF-master/EvoEF --command=ComputeBinding --pdb  {params.wt} --split  {params.part1_part2} > {output.wt}  

        """
    #log:
    #    os.path.abspath(work_dir + "/mutants_structure_generation/EVOEF/structures/{pdb}={chain}={mutations}.log")


rule run_evoEF1_eval_substract:
    input:
        mut = work_dir + "/mutants_structure_scoring/EVOEF/scores/{pdb}={chain}={mutations}.sc",
        wt = work_dir + "/mutants_structure_scoring/EVOEF/scores/{pdb}={chain}={mutations}=WT.sc",
    output:
        mut = work_dir + "/mutants_structure_scoring/EVOEF/scores/{pdb}={chain}={mutations}.ddg",
    script:
        "notebooks/substract_evoef1.r.R"
    #notebook:
    #    "notebooks/substract_evoef1.r.ipynb"


rule run_evoEF1_eval_summary:
    input:
        evoEF1_rez =  aggregate_EVOEF_results,
        muation_data = mutrez + "/mutations_data_4EVOEFmodels.csv",
    output:
        muation_data = mutrez + "/mutations_data_4EVOEFmodels_results.csv",
    #notebook:
    #    "notebooks/collect_evoef1_results.r.ipynb"
    script:
        "notebooks/collect_evoef1_results.r.R"




        
