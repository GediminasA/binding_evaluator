def aggregate_mutation_stems(wildcards):
    checkpoint_output = checkpoints.create_idividual_tasks_4_modeling_with_sequence.get(**wildcards).output[0]
    stems = glob_wildcards(os.path.join(checkpoint_output, "{i,[^\.].....=[A-Z].+}")).i # regex excludes temporary files and WT structures
    return expand(work_dir + "/mutants_structure_scoring/DSSP/scores/complex/{stem}.sum", stem=stems) + \
        expand(work_dir + "/mutants_structure_scoring/DSSP/scores/part/{stem}.sum", stem=stems)

rule collect_binding_terms:
    input:
        aggregate_mutation_stems
