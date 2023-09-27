

### EVALUATION RULES ###

#
rule evaluate_freesasa:
    input: 
        structure = "{directory}/{stem}_{part}_{frame}.pdb",
        container = "containers/freesas.sif"
    output: 
        tsv = temp("{directory}/{stem,[^_]+}_{part,[\d]+}_{frame,[^_]+}_freesasa_initial.tsv"),
    container:
        "containers/freesas.sif"
    shell:
        """
            freesasa  --depth=structure  --output={output}   {input.structure} 
        """
rule clean_freesasa_output:
    input:
        tsv = temp("{directory}/{stem}_{frame}_{part}_freesasa_initial.tsv"),
    output:
        tsv = temp("{directory}/{stem,[^_]+}_{frame,[\d]+}_{part,[^_]+}_freesasa.tsv"),
    params:
        complexid = "{stem}",
        part = "{part}",
        frameid = "{frame}"
    notebook:
        "notebooks/clean_freesasa_output.r.ipynb"

