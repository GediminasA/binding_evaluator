import os
rule repack_docking:
    input:
        initial_complex = work_dir+"/processed/{stem}.pdb"
    output:
        repacked = work_dir+"/rosetta_on_static/prepack/{stem}_0001.pdb"
    params:
        bin = config["rosetta"]["folder"]+"/main/source/bin/" + "docking_prepack_protocol" + config["rosetta"]["binaries_suffix"],
        db = config["rosetta"]["folder"] + "/main/database",
        main_chain = config["main_chain"],
        
    log:
        work_dir+"/rosetta_on_static/prepack/{stem}_prepack.log"
    run:
        part1,part2 = get_interacting_chains(input.initial_complex,params.main_chain)
        part1 = "".join(part1)
        part2 = "".join(part2)
        interactors = part1+"_"+part2
        output_basename = os.path.basename(output.repacked)
        output_dirname = os.path.dirname(output.repacked)
        shell("""
        mkdir -p {output_dirname}
        {params.bin} \
        -in:file:s {input.initial_complex} \
        -out:prefix {output_dirname}/ \
        -ex1 -ex2aro -docking:dock_rtmin \
        -partners {interactors}  &> {log}
        """)


rule run_localdocking:
    input:
        complex = work_dir+"/rosetta_on_static/prepack/{stem}_0001.pdb"
    output:
        complex_str = work_dir+"/rosetta_on_static/localdock/{stem,[^_]+}_{id}.out"
    log:
        complex_log = work_dir+"/rosetta_on_static/localdock/{stem}_{id}.log"
    params:
        bin = config["rosetta"]["folder"]+"/main/source/bin/" + "docking_protocol" + config["rosetta"]["binaries_suffix"],
        main_chain = config["main_chain"],
    run:
        part1,part2 = get_interacting_chains(input.complex,params.main_chain)
        part1 = "".join(part1)
        part2 = "".join(part2)
        interactors = part1+"_"+part2
        shell("""
    {params.bin} \
    -s {input.complex} \
    -out:file:silent {output.complex_str} \
    -nstruct 1 \
    -partners {interactors}   -ex1  -ex2aro  -dock_pert 3 8   \
    &> {log.complex_log}""")

rule run_localdocking_refinament:
    input:
        complex = work_dir+"/rosetta_on_static/localdock/{stem}_{id}.out",
        complexpdb = work_dir+"/rosetta_on_static/prepack/{stem}_0001.pdb"
    output:
        complex_str = work_dir+"/rosetta_on_static/localdockref/{stem,[^_]+}_{id}.out"
    log:
        complex_log = work_dir+"/rosetta_on_static/localdockref/{stem}_{id}.log"
    params:
        bin = config["rosetta"]["folder"]+"/main/source/bin/" + "docking_protocol" + config["rosetta"]["binaries_suffix"],
        main_chain = config["main_chain"],
    run:
        part1,part2 = get_interacting_chains(input.complexpdb,params.main_chain)
        part1 = "".join(part1)
        part2 = "".join(part2)
        interactors = part1+"_"+part2
        shell("""
    {params.bin} \
    -in:file:silent {input.complex} \
    -out:file:silent {output.complex_str} \
    -nstruct 1 \
    -partners {interactors}   -ex1  -ex2aro   -use_input_sc  -docking_local_refine    \
    &> {log.complex_log}""")

rule run_localdocking_refinament_eval:
    input:
        complex = work_dir+"/rosetta_on_static/localdockref/{stem}_{id}.out",
        complexpdb = work_dir+"/rosetta_on_static/prepack/{stem}_0001.pdb"
    output:
        complex_score = work_dir+"/rosetta_on_static/localdockref_eval/{stem,[^_]+}_{id}.sc"
    log:
        complex_log = work_dir+"/rosetta_on_static/localdockref_eval/{stem}_{id}.log"
    params:
        bin = config["rosetta"]["folder"]+"/main/source/bin/" + "InterfaceAnalyzer" + config["rosetta"]["binaries_suffix"],
        main_chain = config["main_chain"],
    run:
        part1,part2 = get_interacting_chains(input.complexpdb,params.main_chain)
        part1 = "".join(part1)
        part2 = "".join(part2)
        interactors = part1+"_"+part2
        shell("""
    {params.bin} \
    -in:file:silent {input.complex} \
    -out:file:score_only {output.complex_score} \
    -nstruct 1 \
    -interface {interactors}  -pack_separated     \
    &> {log.complex_log}""")


local_docking_runs_ids = list(range(1,config['rosetta']['local_docking_runs'],1))
rule rosetta_static_summary:
    input:
        #expand(work_dir+"/rosetta_on_static/{stem}_0001.pdb",stem=pdb_stems)
        expand(work_dir+"/rosetta_on_static/localdockref_eval/{stem}_{id}.sc",stem=pdb_stems,id=local_docking_runs_ids)
    params:
        id = "{stem}"
    output:
        work_dir+"/scores/rosetta_static_{stem,[^_]+}_alldata.csv",
        work_dir+"/scores/rosetta_static_{stem,[^_]+}_summary.csv",
    notebook:
        "notebooks/summarise_rosetta_on_static.r.ipynb"

rule run_extract_pdb:
    input:
        complex_str = work_dir+"/rosetta_on_static/localdockref/{stem}_{id}.out",
    output:
        complex_str = work_dir+"/rosetta_on_static/localdockref/{stem,[^_]+}_{id}.pdb",
    log:
        complex_log = work_dir+"/rosetta_on_static/localdockref/{stem}_{id}_get_pdb.log"
    params:
        bin = config["rosetta"]["folder"]+"/main/source/bin/" + "InterfaceAnalyzer" + config["rosetta"]["binaries_suffix"],
        db = config["rosetta"]["folder"] + "/main/database",
        complex_str_dir = work_dir+"/rosetta_on_static/localdockref/{stem}_{id}_pdb" 
    threads: 1
    shell:"""
    rm -rf {params.complex_str_dir}
    mkdir -p {params.complex_str_dir}
{params.bin} \
-in:file:silent {input.complex_str} \
-out:prefix {params.complex_str_dir}/ &> {log.complex_log}
mv {params.complex_str_dir}/*.pdb  {output.complex_str}
# """


rule run_prodigy_after_rosetta:
    input:
        complex_str = work_dir+"/rosetta_on_static/localdockref/{stem,[^_]+}_{id}.pdb",
    output:
        complex_prodigy = work_dir+"/rosetta_on_static/localdockref/{stem,[^_]+}_{id}_prodigy.sc",
    params:
        id= "{id}",
        stem = "{stem}",
        main_chain=config["main_chain"]
    log:
        complex_prodigy = work_dir+"/rosetta_on_static/localdockref/{stem}_{id}_prodigy.log",
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
                    echo NA  {part1} {part2} {params.stem} {params.id} > {output}
                else
                    echo $out | cut -f 2 -d " " | sed "s/$/ {part1} {part2} {params.stem} {params.id}/g" &>> {output}
                    echo $out >> {log} 
                fi
            """
        )

rule rosetta_static_summary_prodogy:
    input:
        expand(work_dir+"/rosetta_on_static/localdockref/{stem}_{id}_prodigy.sc",stem=pdb_stems,id=local_docking_runs_ids)
    params:
        id = "{stem}"
    output:
        work_dir+"/scores/rosetta_static_prodigy_{stem,[^_]+}_alldata.csv",
        work_dir+"/scores/rosetta_static_prodigy_{stem,[^_]+}_summary.csv",
    notebook:
        "notebooks/summarise_rosetta_prodigy_on_static.r.ipynb"