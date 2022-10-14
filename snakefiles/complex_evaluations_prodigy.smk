
rule run_prodigy:
    input:
        complex_str = work_dir+"/processed/{stem}.pdb",
        interacting_groups = work_dir+"/processed_info/{stem,[^_]+}_interactigGroups.tsv"
    output:
        complex_prodigy = work_dir + "/static/{stem,[^_]+}_prodigy_ini.out"
    params:
        id= "{stem}",
    log:
        work_dir + "/static/{stem}_prodigy_ini.log"
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

rule parse_prodigy:
    input:
        complex_prodigy = work_dir + "/static/{stem}_prodigy_ini.out"
    output:
        complex_prodigy = work_dir + "/static/{stem}_prodigy.tsv"
    params:
        id= "{stem}",
    notebook:
        "notebooks/parse_prodigy.r.ipynb"

