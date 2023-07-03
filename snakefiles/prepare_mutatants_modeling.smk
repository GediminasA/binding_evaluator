### PROCESS CSV file with mutations ###
import pandas as pd 
df_muts = pd.read_csv(config["mutants"])
pdb_template_stems = list(set(df_muts.PDB))
seqs_4mutations = list(set([s.strip() for s in df_muts.Template]))
pair4mut = ["=".join(l) for l in  list(set(list(zip([s.strip() for s in df_muts.PDB],[s.strip() for s in df_muts.Template]))))]
mutrez = work_dir + "/rezults"
nconf = config["conformers_to_evaluate"]

rule mutants_targets_EVOEF:
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
            work_dir+"/mutants_sequence_generation/{pair}=template_mutation_data__4models.csv",
            pair = pair4mut),
        ddg_results = mutrez + "/mutations_data_4EVOEFmodels_results.csv"
        

rule collect_data_befor_modelling:
    input:
        cleaned_haps = expand(
            work_dir+"/mutants_sequence_generation/{pair}=template_mutation_data__cleanhap.csv",
            pair = pair4mut),
        evoef_data = expand(
            work_dir+"/mutants_sequence_generation/{pair}=template_mutation_data__4models.csv",
            pair = pair4mut)
    output:
        cleaned_haps = mutrez + "/mutations_data_cleanHap.csv",
        evoef_data = mutrez + "/mutations_data_4models.csv"
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
        work_dir+"/mutants_sequence_generation/{pdb,[^=]+}={seqtempl,[^=]+}=template_mutation_data__4models.csv",
        work_dir+"/mutants_sequence_generation/{pdb,[^=]+}={seqtempl,[^=]+}=template_mutation_data__cleanhap.csv"
    # notebook:
    #   "notebooks/get_template_mutation_data.r.ipynb"   
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
        muation_data = mutrez + "/mutations_data_4models.csv",
    output:
        muation_data = mutrez + "/mutations_data_4EVOEFmodels_results.csv",
    #notebook:
    #    "notebooks/collect_evoef1_results.r.ipynb"
    script:
        "notebooks/collect_evoef1_results.r.R"


#### rules to prepare for template based modelling modelling
        
checkpoint create_idividual_tasks_4_modeling_with_sequence: #for tools that can model full sequence
    input:
        mutrez + "/mutations_data_4models.csv",
    output:
        directory(work_dir + "/mutants_structure_generation/TEMPLATES/todoList") # writes out small text file for each mutants to be generated
    notebook:
        "notebooks/create_idividual_tasks_4_modeling_based_on_templates.py.ipynb"

rule get_template_sequence:
    input:
        coordinates = work_dir+"/mutants_sequence_generation/{pdb}_vs_templates_blastp.txt",
        pdbfasta = work_dir+"/processed_info/{pdb}.fasta",
        pdbtemplates = config["mutants_templates"],
        todo = work_dir + "/mutants_structure_generation/TEMPLATES/todoList/{pdb}={chain}={mutations}",
    log:
        work_dir + "/mutants_structure_generation/TEMPLATES/sequences/{pdb}={chain}={mutations}.log"
    output:
        todo = work_dir + "/mutants_structure_generation/TEMPLATES/sequences/{pdb}={chain}={mutations}.fasta" # this contains all template sequence
    # notebook:
    #    "notebooks/get_template_sequence.py.ipynb"
    script:
        "notebooks/get_template_sequence.py.py"


def aggregate_TEMPLATE_seqs(wildcards):
    checkpoint_output = checkpoints.create_idividual_tasks_4_modeling_with_sequence.get(**wildcards).output[0]
    stems = glob_wildcards(os.path.join(checkpoint_output, "{i,[^\.].+}")).i #rehex prevents snakemake temps to be included
    return(expand(work_dir + "/mutants_structure_generation/TEMPLATES/sequences/{stem}.fasta",stem=stems))

def aggregate_TEMPLATE_promod_models(wildcards):
    checkpoint_output = checkpoints.create_idividual_tasks_4_modeling_with_sequence.get(**wildcards).output[0]
    stems = glob_wildcards(os.path.join(checkpoint_output, "{i,[^\.].+}")).i #rehex prevents snakemake temps to be included
    return(expand(work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{stem}.pdb",stem=stems))

def aggregate_TEMPLATE_promod_models_ready_4_eval(wildcards):
    checkpoint_output = checkpoints.create_idividual_tasks_4_modeling_with_sequence.get(**wildcards).output[0]
    stems = glob_wildcards(os.path.join(checkpoint_output, "{i,[^\.].+}")).i #rehex prevents snakemake temps to be included
    return(expand(work_dir + "/evaluation/starting_structures/{stem}_1.pdb",stem=stems))

def aggregate_TEMPLATE_promod_models_ready_4_eval_with_conformers(wildcards):
    checkpoint_output = checkpoints.create_idividual_tasks_4_modeling_with_sequence.get(**wildcards).output[0]
    stems = glob_wildcards(os.path.join(checkpoint_output, "{i,[^\.].+}")).i #rehex prevents snakemake temps to be included
    return(expand(work_dir + "/evaluation/starting_structures/{stem}_{id}.pdb",stem=stems, id = range(1, nconf+2)))

rule model_mutants_promod:
    input:
        "containers/promod.sif",
        structure = work_dir+"/processed/{pdb}.pdb",
        sequence = work_dir + "/mutants_structure_generation/TEMPLATES/sequences/{pdb}={chain}={mutations}.fasta" 
    output:
        model = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb}={chain}={mutations}.pdb"
    log:
        work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb}={chain}={mutations}.log"
    container: "containers/promod.sif"
    threads: 8
    shell:
        """
        export OPENMM_CPU_THREADS={threads}
        PYTHONPATH=covid-lt covid-lt/bin/promod-model --simulate --trim --template {input.structure} --sequences {input.sequence} 1> {output.model} 2> {log} 
        """ 


rule copy_for_evaluation_static:
    input:
        model = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb}={chain}={mutations}.pdb"
    output:
        model = work_dir + "/evaluation/starting_structures/{pdb}={chain}={mutations}_1.pdb"
    shell:
        """
            cp {input} {output}
        """

    

rule generate_conformations_CABS_part1:
    input:
        "{stem}.pdb"        
    output:
        ["{stem}_CABS/output_pdbs/model_"+str(id)+".pdb" for id in range(0,nconf)]
    params:
        lwdir = "{stem}_CABS",
        nconf = nconf
    conda:
        "../envs/CABS.yaml"
    shell:
        """
            echo {params.lwdir}
            mkdir -p  {params.lwdir}
            CABSflex --dssp mkdssp -o M -a 2 -y {params.nconf} -i {input} -w {params.lwdir}
        """
            #mkdir -p {prams.lwdir}
            #CABSflex --dssp mkdssp -o M -a 2 -y 10 -i {input} -w {params.lwdir} 

rule generate_conformations_CABS_part2:
    input:
        ["{stem}_CABS/output_pdbs/model_"+str(id)+"_FULL.pdb" for id in range(0,nconf)]
    output:
        ["{stem}_"+str(id)+"_CABS.pdb" for id in range(2,nconf+2)] #second number is reserved for WT
    run:
        n = len(input)
        print(n)
        for i in range(0,n):
            inf = input[i]
            outputf = output[i]
            shell("cp {inf} {outputf}")


rule copy_for_evaluation_dynamic:
    input:
        model = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{stem}_{id}_CABS.pdb"
    output:
        model = work_dir + "/evaluation/starting_structures/{stem}_{id,[^1]\d?|1\d+}.pdb" # one is reserved for static
    shell:
        """
            cp {input} {output}
        """

rule get_full_atoms_structure_from_CA:
    input:
        "{stem}.pdb" 
    output:
        "{stem}_FULL.pdb"
    conda:
        "../envs/modeller.yaml"
    shell:
        """
            python2.7 external/ca2all.py  -i {input} -o {output}
        """



#### aggregates for evaluation collection and testing


def aggregate_TEMPLATE_promod_models_prodigy(wildcards):
    checkpoint_output = checkpoints.create_idividual_tasks_4_modeling_with_sequence.get(**wildcards).output[0]
    stems = glob_wildcards(os.path.join(checkpoint_output, "{i,[^\.].+}")).i #rehex prevents snakemake temps to be included
    return(expand(
        
        work_dir + "/evaluation/scores/{stem}_{id}_prodigy.tsv",
        stem=stems, id=[1]))

def aggregate_TEMPLATE_promod_models_prodigy_with_conformers(wildcards):
    checkpoint_output = checkpoints.create_idividual_tasks_4_modeling_with_sequence.get(**wildcards).output[0]
    stems = glob_wildcards(os.path.join(checkpoint_output, "{i,[^\.].+}")).i #rehex prevents snakemake temps to be included
    return(expand(
        
        work_dir + "/evaluation/scores/{stem}_{id}_prodigy.tsv",
        stem=stems, id=range(1,12)))

def aggregate_TEMPLATE_promod_models_evoef1(wildcards):
    checkpoint_output = checkpoints.create_idividual_tasks_4_modeling_with_sequence.get(**wildcards).output[0]
    stems = glob_wildcards(os.path.join(checkpoint_output, "{i,[^\.].+}")).i #rehex prevents snakemake temps to be included
    return(expand(
        
        work_dir + "/evaluation/scores/{stem}_1_evoef1.tsv",
        stem=stems))

def aggregate_TEMPLATE_promod_models_evoef1_with_conformers(wildcards):
    checkpoint_output = checkpoints.create_idividual_tasks_4_modeling_with_sequence.get(**wildcards).output[0]
    stems = glob_wildcards(os.path.join(checkpoint_output, "{i,[^\.].+}")).i #rehex prevents snakemake temps to be included
    return(expand(
        
        work_dir + "/evaluation/scores/{stem}_{id}_evoef1.tsv",
        stem=stems, id=range(1,12)))


rule mutants_targets_templates:
    input:
        #todolist = work_dir + "/mutants_structure_generation/TEMPLATES/todoList",
        #target = aggregate_TEMPLATE_promod_models_prodigy
        # modelling_templates = aggregate_TEMPLATE_seqs,
        promod_models = aggregate_TEMPLATE_promod_models,
        promod_models_copied_4_eval = aggregate_TEMPLATE_promod_models_ready_4_eval,
        promod_models_copied_4_eval_with_conformers = aggregate_TEMPLATE_promod_models_ready_4_eval_with_conformers,
        promod_models_prodigy_evals = aggregate_TEMPLATE_promod_models_prodigy,
        promod_models_prodigy_evals_with_conformers = aggregate_TEMPLATE_promod_models_prodigy_with_conformers,
        promod_models_evoef1_evals = aggregate_TEMPLATE_promod_models_evoef1,
        promod_models_evoef1_evals_with_conformers = aggregate_TEMPLATE_promod_models_evoef1_with_conformers

rule get_summary_of_binding:
    input:
        promod_models_prodigy_evals = aggregate_TEMPLATE_promod_models_prodigy_with_conformers,
        promod_models_evoef1_evals = aggregate_TEMPLATE_promod_models_evoef1_with_conformers
    output:
        ddg_results_on_promod_full = mutrez + "/promod_models_results_full.csv",
        ddg_results_on_promod_main = mutrez + "/promod_models_results_main.csv",
    notebook:
        "notebooks/collect_results_4promod.r.ipynb"

