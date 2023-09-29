### PROCESS CSV file with mutations ###
import pandas as pd 
import tempfile
df_muts = pd.read_csv(config["mutants"])
df_muts = df_muts.dropna()
pdb_template_stems = list(set(df_muts.PDB))
seqs_4mutations = list(set([str(s).strip() for s in df_muts.Template]))
pair4mut = ["=".join(l) for l in  list(set(list(zip([str(s).strip() for s in df_muts.PDB],[str(s).strip() for s in df_muts.Template]))))]
mutrez = work_dir + "/rezults"
nconf = config["conformers_to_evaluate"]

rule mutants_targets_EVOEF:
    input:
        # pdb_templates = expand(
        #     work_dir+"/processed/{stem}.pdb",
        #     stem=pdb_template_stems),
        pdb_test = expand(
            work_dir+"/initial_cleanup/{stem}_noambi.pdb",
            stem=pdb_template_stems),
        #prodigy_test/initial_cleanup/1S1Q_noambi.pdb
        # pdb_sequences = expand(
        #     work_dir+"/processed_info/{stem}.fasta",
        #     stem=pdb_template_stems),
        # pdb_sequences_aligned_to_templates = expand(
        #     work_dir+"/mutants_sequence_generation/{stem}_vs_templates_blastp.txt",
        #     stem=pdb_template_stems),
        # pdb_templates_map = expand(
        #     work_dir+"/mutants_sequence_generation/{pair}.csv", #map input seq to pdb map 
        #     pair = pair4mut),
        # pdb_tempplate_mutation_data = expand(
        #     work_dir+"/mutants_sequence_generation/{pair}=template_mutation_data__4models.csv",
        #     pair = pair4mut),
        # ddg_results = mutrez + "/mutations_data_4EVOEFmodels_results.csv"
        

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
    # notebook:
    #     "notebooks/collect_data_before_modelling.r.ipynb"   
    script:
        "notebooks/collect_data_before_modelling.r.R"   


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
        #work_dir+"/mutants_sequence_generation/{pdb}={seqtempl,[^=]+}.csv",
        work_dir+"/mutants_sequence_generation/{pdb,[^=]+}={seqtempl,[^=]+}.csv",
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
    notebook:
        "notebooks/get_template_mutation_data.r.ipynb"   
    #script:
    #     "notebooks/get_template_mutation_data.r.R"   

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
    container:
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
    container:
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
        todo = work_dir + "/mutants_structure_generation/TEMPLATES/sequences/{pdb}={chain}={mutations,[^_]+}.fasta" # this contains all template sequence
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
        model = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb}={chain}={mutations,[^_]+}_before_faspr.pdb"
    log:
        work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb}={chain}={mutations}.log"
    container:
        "containers/promod.sif"
    threads: 8
    shell:
        """
        export OPENMM_CPU_THREADS={threads}
        PYTHONPATH=covid-lt covid-lt/bin/promod-model --simulate --trim --template {input.structure} --sequences {input.sequence} > {output.model} 2> {log}
        """

rule model_mutants_promod_faspr:
    input:
        structure = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb}={chain}={mutations,[^_]+}_before_faspr.pdb",
        container = "containers/faspr.sif"
    output:
        work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb}={chain}={mutations,[^_]+}.pdb"
    log:
        work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb}={chain}={mutations}.log"
    container:
        "containers/faspr.sif"
    shell:
        "FASPR -i {input.structure} -o {output} 2>&1 | cat >> {log}"

rule model_mutants_faspr_denovo:
    input:
        structure = work_dir + "/pdb_proc/pristine/{pdb}.pdb",
        groups = work_dir + "/processed_info/{pdb}_interactigGroups.tsv",
        container = "containers/faspr.sif"
    output:
        work_dir + "/mutants_structure_generation/TEMPLATES/faspr_models/{pdb}={chain}={mutations,[^_]+}_before_promod.pdb"
    container:
        "containers/faspr.sif"
    shell:
        """
        rm -f {output}

        MUTATIONS=
        if [ "{wildcards.mutations}" != nan ]
        then
            MUTATIONS=$(echo {wildcards.mutations} \
                | tr + ' ' \
                | xargs -n 1 echo \
                | awk '{{ print substr($0, 0, 1) "{wildcards.chain}" substr($0, 2) }}' \
                | xargs -i echo --replace {})
        fi

        TMPFILE=$(mktemp --suffix .pdb)
        covid-lt-new/bin/pdb_select --first-model --chain $(cat {input.groups} | cut -f 1,2 | sed 's/[\t,]//g') {input.structure} \
            | PYTHONPATH=covid-lt-new covid-lt-new/bin/pdb_resolve_alternate_locations > $TMPFILE
        PYTHONPATH=covid-lt-new covid-lt-new/bin/pdb_atom2fasta --with-initial-gaps $TMPFILE \
            | covid-lt-new/bin/fasta2pdb_seqres \
            | covid-lt-new/bin/pdb_mutate_seqres $MUTATIONS --trim-gaps \
            | covid-lt-new/bin/pdb_seqres2fasta \
            | grep -v '^>' \
            | grep -o . \
            | grep -v X \
            | xargs echo \
            | sed 's/ //g' \
            | FASPR -i $TMPFILE -s /dev/stdin -o {output} || true
        test -e {output} || echo -n > {output}
        rm $TMPFILE
        """

# This rule is needed to add hydrogens to FASPR-optimized structures as FASPR does not do that itself
rule model_mutants_promod_after_faspr:
    input:
        structure = work_dir + "/mutants_structure_generation/TEMPLATES/faspr_models/{pdb}={chain}={mutations,[^_]+}_before_promod.pdb",
        container = "containers/promod.sif"
    output:
        work_dir + "/mutants_structure_generation/TEMPLATES/faspr_models/{pdb}={chain}={mutations,[^_]+}.pdb"
    container:
        "containers/promod.sif"
    shell:
        "PYTHONPATH=covid-lt-new covid-lt-new/bin/promod-fix-pdb --do-not-fill-gaps --simulate {input.structure} > {output}"

rule copy_for_evaluation_static:
    input:
        #model = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb}={chain}={mutations}_DeepRefine.pdb"
        model = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{pdb}={chain}={mutations}.pdb"
    output:
        model = work_dir + "/evaluation/starting_structures/{pdb}={chain}={mutations,[^_]+}_1.pdb"
    shell:
        """
            cp {input} {output}
        """

    

rule generate_conformations_CABS_part1:
    input:
        "{stem}.pdb"        
    output:
        #["{stem}_CABS/output_pdbs/model_"+str(id)+".pdb" for id in range(0,nconf)]
        ["{stem}_CABS/output_pdbs/model_"+str(id)+"_FULL.pdb" for id in range(0,nconf)]
    params:
        lwdir = "{stem}_CABS",
        nconf = nconf
    conda:
        "../envs/CABS.yaml"
    shell:
        """
            echo {params.lwdir}
            mkdir -p  {params.lwdir}
            CABSflex --dssp mkdssp -o M -i {input} -w {params.lwdir}
            #CABSflex --dssp mkdssp -o M -a 2 -y {params.nconf} -i {input} -w {params.lwdir}
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


rule generate_conformations_Galaxy_part1:
    input:
        "{stem}.pdb",        
        "containers/GalaxyRefineComplex.sif"
    output:
        ["{stem}_GalaxyRefineComplex/galaxyref/model/model_"+str(id)+".pdb" for id in range(1,nconf+1)]
    params:
        lwdir = "{stem}_GalaxyRefineComplex",
        nconf = nconf,
    container: "containers/GalaxyRefineComplex.sif"
    threads: 16
    shell:
        """
            mkdir -p  {params.lwdir}
            startd=`pwd`
            cd {params.lwdir}
            export GALAXY_HOME=/opt/GalaxyRefineComplex
            export NSLOTS={threads}
            /opt/GalaxyRefineComplex/bin/GalaxyRefineComplex.ubuntu1604 --protocol1 -p $startd/{input[0]}  -t galaxyref 
        """

rule generate_conformations_GALAXY_part2:
    input:
        ["{stem}_GalaxyRefineComplex/galaxyref/model/model_"+str(id)+".pdb" for id in range(1,nconf+1)]
    output:
        ["{stem}_"+str(id)+"_GalaxyRefineComplex.pdb" for id in range(2,nconf+2)] #second number is reserved for WT
    run:
        n = len(input)
        print(n)
        for i in range(0,n):
            inf = input[i]
            outputf = output[i]
            shell("cp {inf} {outputf}")

rule copy_for_evaluation_dynamic:
    input:
        model = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{stem}_{id}_GalaxyRefineComplex.pdb"
        #model = work_dir + "/mutants_structure_generation/TEMPLATES/promod_models/{stem}_{id}_CABS_DeepRefine.pdb"
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
def get_tmp_dir(wildcards):
 out = tempfile.TemporaryDirectory(dir = "/dev/shm").name
 return(out)

rule DeepRefine:
    input:
        "{stem}/{name,[^/]+}.pdb"
    output:
        file = "{stem}/{name,[^/]+}_DeepRefine.pdb",
        score = "{stem}/{name,[^/]+}_DeepRefine_plddt.tsv"
    threads: 4
    params:
        tmp =get_tmp_dir,
        prog_path = os.getcwd() + "/DeepRefine",
        cur_path = os.getcwd()
    conda:
        "../envs/DeepRefine.yaml"
    shell:
        """
            startdir=`pwd`
            indir={params.tmp}/in/input/
            infile={params.tmp}/in/input/input.pdb
            outdir={params.tmp}/out/
            outfile={params.tmp}/out/input/input_refined.pdb
            outscore={params.tmp}/out/input/input_refined_plddt.csv
            
            mkdir -p $indir
            mkdir -p $outdir
            cp  {input} $infile

            DR_DIR={params.prog_path}
            
            ckpt_dir="$DR_DIR"/project/checkpoints/EGR_All_Atom_Models
            ckpt_name=LitPSR_EGR_AllAtomModel1_Seed42.ckpt
            atom_selection_type=all_atom
            seed=42
            nn_type=EGR
            graph_return_format=dgl
            
            export PYTHONPATH=$DR_DIR
            cd $DR_DIR/project
            python3  lit_model_predict.py --device_type gpu --num_devices 1 --num_compute_nodes 1 --num_workers 1 --batch_size 1 --input_dataset_dir {params.tmp}/in --output_dir {params.tmp}/out  --ckpt_dir "$ckpt_dir" --ckpt_name "$ckpt_name" --atom_selection_type "$atom_selection_type" --seed "$seed" --nn_type "$nn_type" --graph_return_format "$graph_return_format" --perform_pos_refinement 
            
            cp $outfile $startdir/{output.file}
            cp $outscore $startdir/{output.score}
            rm -r {params.tmp} 
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
        stem=stems, id=range(1, nconf+2)))

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
        stem=stems, id = range(1, nconf+2)))


rule mutants_targets_templates:
    input:
        target = aggregate_TEMPLATE_promod_models_prodigy,
        # modelling_templates = aggregate_TEMPLATE_seqs,
        #promod_models = aggregate_TEMPLATE_promod_models,
        #promod_models_copied_4_eval = aggregate_TEMPLATE_promod_models_ready_4_eval,
        #promod_models_copied_4_eval_with_conformers = aggregate_TEMPLATE_promod_models_ready_4_eval_with_conformers,
        #promod_models_prodigy_evals = aggregate_TEMPLATE_promod_models_prodigy,
        #promod_models_prodigy_evals_with_conformers = aggregate_TEMPLATE_promod_models_prodigy_with_conformers,
        #promod_models_evoef1_evals = aggregate_TEMPLATE_promod_models_evoef1,
        #promod_models_evoef1_evals_with_conformers = aggregate_TEMPLATE_promod_models_evoef1_with_conformers

rule collect_ddG_binding:
    input:
        promod_models_prodigy_evals = aggregate_TEMPLATE_promod_models_prodigy_with_conformers,
        promod_models_evoef1_evals = aggregate_TEMPLATE_promod_models_evoef1_with_conformers,
        ddg = work_dir + "/rezults/mutation_ddg_predictions.csv"
    output:
        ddg_results_on_promod_full = mutrez + "/promod_models_results_full.csv",
        ddg_results_on_promod_main = mutrez + "/promod_models_results_main.csv",
    notebook:
        "notebooks/collect_results_4promod.r.ipynb"

rule get_summary_of_binding:
    input:
        mutdef = mutrez + "/mutations_data_4models.csv",
        ddg_results_on_promod_main = mutrez + "/promod_models_results_main.csv",
    output:
        rez = mutrez + "/sequence_variants_per_ddG.csv",
    notebook:
        "notebooks/match_ddg_with_rezults.r.ipynb"

