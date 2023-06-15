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
            work_dir+"/mutants_sequence_generation/{pair}.csv",
            pair = pair4mut),
        pdb_tempplate_mutation_data = expand(
            work_dir+"/mutants_sequence_generation/{pair}=template_mutation_data__EVOEF.csv",
            pair = pair4mut),
        pdb_tempplate_mutation_data_hap = expand(
            work_dir+"/mutants_sequence_generation/{pair}=template_mutation_data__cleanhap.csv",
            pair = pair4mut)
        

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




