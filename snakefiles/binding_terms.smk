def aggregate_mutation_stems(wildcards):
    checkpoint_output = checkpoints.create_idividual_tasks_4_modeling_with_sequence.get(**wildcards).output[0]
    stems = glob_wildcards(os.path.join(checkpoint_output, "{i,[^\.].....=[A-Z].+}")).i # regex excludes temporary files and WT structures
    return expand(work_dir + "/mutants_structure_scoring/CADscore/scores/{stem}.sc", stem=stems) + \
        expand(work_dir + "/mutants_structure_scoring/DSSP/scores/complex/{stem}.sum", stem=stems) + \
        expand(work_dir + "/mutants_structure_scoring/DSSP/scores/part/{stem}.sum", stem=stems) + \
        expand(work_dir + "/mutants_structure_scoring/EvoEF1/scores/{stem}.diff", stem=stems) + \
        expand(work_dir + "/mutants_structure_scoring/OpenMM/scores/{stem}.diff", stem=stems) + \
        expand(work_dir + "/mutants_structure_scoring/PROVEAN/scores/{stem}.sum", stem=stems)


def aggregate_mutation_stems4dgg(wildcards):
    checkpoint_output = checkpoints.create_idividual_tasks_4_modeling_with_sequence.get(**wildcards).output[0]
    stems = glob_wildcards(os.path.join(checkpoint_output, "{i,[^\.].....=[A-Z].+}")).i # regex excludes temporary files and WT structures
    return(expand(work_dir + "/mutants_structure_scoring/all_terms/{stem}.csv", stem = stems))

rule collect_binding_terms4ddg:
    input:
        cad = work_dir + "/mutants_structure_scoring/CADscore/scores/{stem}.sc",
        dssp_complex = work_dir + "/mutants_structure_scoring/DSSP/scores/complex/{stem}.sum",
        dssp_part = work_dir + "/mutants_structure_scoring/DSSP/scores/part/{stem}.sum",
        evoef = work_dir + "/mutants_structure_scoring/EvoEF1/scores/{stem}.diff",
        openmm = work_dir + "/mutants_structure_scoring/OpenMM/scores/{stem}.diff",
        provean = work_dir + "/mutants_structure_scoring/PROVEAN/scores/{stem}.sc",
        #aggregate_mutation_stems
    output:
        work_dir + "/mutants_structure_scoring/all_terms/{stem}.csv"
    # notebook:
    #     "notebooks/collect_terms4ddg.r.ipynb"
    script:
       "notebooks/collect_terms4ddg.r.R"

rule merge_colected_binding_terms:
    input:
        aggregate_mutation_stems4dgg
    output:
        table = work_dir + "/rezults/mutation_terms_4ddg.csv"
    notebook:
        "notebooks/merge_ddgrez.r.ipynb"
    # script:
    #     "notebooks/merge_ddgrez.r.R"

rule collect_binding_terms:
    input:
        aggregate_mutation_stems
    output:
        work_dir + "/rezults/mutation_terms.csv"
    shell:
        """
        (
            (
                echo label,CADscore,dS,SA_com,SA_part,ddG_EvoEF
                ls -1 {work_dir}/mutants_structure_scoring/OpenMM/scores/*.sc \
                    | head -n 1 \
                    | xargs cut -f 1
                echo CS
            ) \
                | xargs echo \
                | sed 's/ /,/g'
            ls -1 {work_dir}/mutants_structure_scoring/CADscore/scores/*.sc \
                | xargs -i basename {{}} .sc \
                | grep -v '=nan$' \
                | while read LABEL
                  do
                    (
                        echo $LABEL
                        awk '{{print $5 " " $7 - $6}}' < {work_dir}/mutants_structure_scoring/CADscore/scores/$LABEL.sc
                        cat {work_dir}/mutants_structure_scoring/DSSP/scores/complex/$LABEL.sum
                        cat {work_dir}/mutants_structure_scoring/DSSP/scores/part/$LABEL.sum
                        cat {work_dir}/mutants_structure_scoring/EvoEF1/scores/$LABEL.diff
                        cut -f 2 {work_dir}/mutants_structure_scoring/OpenMM/scores/$LABEL.diff
                        cat {work_dir}/mutants_structure_scoring/PROVEAN/scores/$LABEL.sum
                    ) \
                        | xargs echo \
                        | sed 's/ /,/g'
                  done
        ) > {output}
        """

rule predict_ddG:
    input:
        table = work_dir + "/rezults/mutation_terms_4ddg.csv",
        model = "covid-lt-new/binding-evaluator-model-our.RData",
        container = "containers/r-cran.sif"
    output:
        work_dir + "/rezults/mutation_terms_predicted.csv"
    container:
        "containers/r-cran.sif"
    shell:
        """
        covid-lt-new/bin/random-forest {input.table} --input-format csv --input-model {input.model} > {output}
        """

rule extract_results_ddg:
    input:
        table = work_dir + "/rezults/mutation_terms_4ddg.csv",
        ddg = work_dir + "/rezults/mutation_terms_predicted.csv"
    output:
        ddg = work_dir + "/rezults/mutation_ddg_predictions.csv"
    # notebook:
    #     "notebooks/get_mutation_ddg_predictions.r.ipynb"
    script:
        "notebooks/get_mutation_ddg_predictions.r.R"

