
rule run_prodigy:
    input:
        complex_str = work_dir+"/processed/{stem}.pdb"
    output:
        complex_prodigy = work_dir + "/prodigy_static/{stem,[^_]+}.sc"
    params:
        id= "{stem}",
        main_chain=config["main_chain"]
    log:
        work_dir + "/prodigy_static/{stem}.log"
    container: None 
    run:
        part1,part2 = get_interacting_chains(input.complex_str,params.main_chain)
        part1 = ",".join(part1)
        part2 = ",".join(part2)
        part14out = "".join(part1)
        part24out = "".join(part2)
        shell(
            """
                set +e
                out=`prodigy  {input} --selection {part1} {part2} -q `
                exitcode=$?
                if [ $exitcode -eq 1 ]
                then
                    echo crashed {params.id} >>{log}
                    echo NA  {part14out} {part24out} {params.id} > {output}
                else
                    echo $out | cut -f 2 -d " " | sed "s/$/ {part14out} {part24out} {params.id}/g" &>> {output}
                    echo $out >> {log} 
                fi
            """
        )
