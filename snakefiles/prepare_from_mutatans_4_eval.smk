#aggregate_TEMPLATE_promod_models_ready_4_eval


rule run_prodigy_mutants:
    input:
        complex_str = work_dir + "/evaluation/starting_structures/{pdb}={chain}={mutations}_{frame}.pdb",
        interacting_groups = work_dir+"/processed_info/{pdb,[^_]+}_interactigGroups.tsv"
    output:
        complex_prodigy = work_dir + "/evaluation/scores/{pdb}={chain}={mutations,[^_]+}_{frame,[^_]+}_prodigy_ini.out"
    params:
        id= "{pdb}={chain}={mutations}_{frame}",
    log:
        work_dir + "/evaluation/scores/{pdb}={chain}={mutations}_{frame}_prodigy_ini.log"
    run:
        part1,part2 = get_interacting_chains(input.interacting_groups)
        part1 = ",".join(part1)
        part2 = ",".join(part2)
        part14out = "".join(part1)
        part24out = "".join(part2)
        shell(
            """
                set +e
                out=`prodigy  {input.complex_str} --selection {part1} {part2} `
                exitcode=$?
                if [ $exitcode -eq 1 ]
                then
                    echo crashed {params.id} >>{log}
                    echo NA   > {output}
                else
                    echo $out &>> {output}
                    echo $out >> {log} 
                fi
            """
        )

rule parse_prodigy_mutants:
    input:
        complex_prodigy = work_dir + "/evaluation/scores/{pdb}={chain}={mutations}_{frame}_prodigy_ini.out"
    output:
        complex_prodigy = work_dir + "/evaluation/scores/{pdb}={chain}={mutations}_{frame,[^_]+}_prodigy.tsv"
    params:
        id= "{pdb}={chain}={mutations}_{frame}",
    # notebook:
    #     "notebooks/parse_prodigy.r.ipynb"
    script:
        "notebooks/parse_prodigy.r.R"


rule run_evoEF1_eval_mutants:
    input:
        mut = work_dir + "/evaluation/starting_structures/{pdb}={chain}={mutations}_{frame}.pdb",
        interacting_groups = work_dir+"/processed_info/{pdb,[^_]+}_interactigGroups.tsv",
        container = "containers/evoef1.sif"
    output:
        mut = work_dir + "/evaluation/scores/{pdb}={chain}={mutations}_{frame,[^_]+}_evoef1.tsv"
    params:
        mut = os.path.abspath(work_dir + "/evaluation/starting_structures/{pdb}={chain}={mutations}_{frame}.pdb"),
        part1_part2 = get_interacting_chains4evoEF1
    container:
        "containers/evoef1.sif"
    shell:
        """
            /EvoEF-master/EvoEF --command=ComputeBinding --pdb  {params.mut} --split  {params.part1_part2} > {output.mut}  
        """