def aggregate_mutation_stems(wildcards):
    checkpoint_output = checkpoints.create_idividual_tasks_4_modeling_with_sequence.get(**wildcards).output[0]
    stems = glob_wildcards(os.path.join(checkpoint_output, "{i,[^\.].....=[A-Z].+}")).i # regex excludes temporary files and WT structures
    return expand(work_dir + "/mutants_structure_scoring/CADscore/scores/{stem}.sc", stem=stems) + \
        expand(work_dir + "/mutants_structure_scoring/DSSP/scores/complex/{stem}.sum", stem=stems) + \
        expand(work_dir + "/mutants_structure_scoring/DSSP/scores/part/{stem}.sum", stem=stems) + \
        expand(work_dir + "/mutants_structure_scoring/EvoEF1/scores/{stem}.diff", stem=stems) + \
        expand(work_dir + "/mutants_structure_scoring/OpenMM/scores/{stem}.diff", stem=stems) + \
        expand(work_dir + "/mutants_structure_scoring/PROVEAN/scores/{stem}.sc", stem=stems)

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
                        tail -n 1 {work_dir}/mutants_structure_scoring/PROVEAN/scores/$LABEL.sc | cut -f 2
                    ) | xargs echo \
                      | sed 's/ /,/g'
                  done
        ) > {output}
        """

rule predict_ddG:
    input:
        table = work_dir + "/rezults/mutation_terms.csv",
        container = "containers/r-cran.sif"
    output:
        work_dir + "/rezults/mutation_terms_predicted.csv"
    container:
        "containers/r-cran.sif"
    shell:
        """
        covid-lt-new/bin/random-forest {input.table} --input-format csv --input-model covid-lt-new/binding-evaluator-model.RData > {output}
        """
