
rule run_prodigy_pdb2pqr:
    input:
        complex_str = work_dir + "/promod/{id}.pdb"
    output:
        complex_prodigy = work_dir + "/promod_prodigy/{id}.sc"
    params:
        id= "{id}",
        main_chain=config["main_chain"]
    log:
        work_dir + "/promod_prodigy/{id}.log"
    container: None 
    run:
        part1,part2 = get_interacting_chains(input[0],params.main_chain)
        part2 = ",".join(part2)
        shell(
            """
                set +e
                out=`prodigy  {input} --selection {part1} {part2} -q `
                exitcode=$?
                if [ $exitcode -eq 1 ]
                then
                    echo crashed {params.id} >>{log}
                    echo NA  {part1} {part2} {params.id} > {output}
                else
                    echo $out | cut -f 2 -d " " | sed "s/$/ {part1} {part2} {params.id}/g" &>> {output}
                    echo $out >> {log} 
                fi
            """
        )
