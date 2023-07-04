#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import json
from Bio import SeqIO
import copy


# In[ ]:


df_todo = pd.read_csv(snakemake.input.todo,sep=";",names=['mutations', 'PDB','chain','template'])
df_todo


# In[ ]:


fastaf = open(snakemake.input.pdbfasta,"r")
pdb_seqs = {}
for record in SeqIO.parse(fastaf, "fasta"):
    chain = record.id.split("|")[1]
    pdb_seqs[chain] = str(record.seq)


# In[ ]:


fastaf = open(snakemake.input.pdbtemplates,"r")
templates_seqs = {}
for record in SeqIO.parse(fastaf, "fasta"):
    templates_seqs[record.id] = str(record.seq)


# In[ ]:


if snakemake.wildcards.chain != "nan":
    df = pd.read_csv(snakemake.input.coordinates, sep = "\t", names = "qseqid sseqid qstart qend sstart send length pident qseq sseq".split())
    df = df.query("sseqid == @df_todo.template[0]")
    df = df.sort_values(by=["length"], ascending= False)
    data = dict(df.iloc[0])


# In[ ]:


if snakemake.wildcards.chain != "nan":
    gap_ct_q = 0
    gap_ct_s = 0
    spos_seqs = []
    qpos_seqs = []

    for i in range(0, len(data['qseq'])):
        qs = data['qseq'][i]
        ss = data['sseq'][i]
        if qs == '-':
            gap_ct_q += 1
        if ss == '-':
            gap_ct_s += 1
        qpos_aln = (i + 1) - gap_ct_q
        spos_aln = (i + 1) - gap_ct_s
        qpos_seq = data["qstart"] + qpos_aln -1
        spos_seq = data["sstart"] + spos_aln -1
        qpos_seqs.append(qpos_seq)
        spos_seqs.append(spos_seq)



# In[ ]:


if snakemake.wildcards.chain != "nan":
    df = pd.DataFrame.from_dict({"PDB":qpos_seqs,"TEMPLATE":spos_seqs})
    df["CHAIN"] = data["qseqid"].split("|")[1]


# In[ ]:


# prepare dicts
if snakemake.wildcards.chain != "nan":
    PDB_TEMPLATE = dict(zip(df.PDB,df.TEMPLATE))
    TEMPLATE_PDB = dict(zip(df.TEMPLATE, df.PDB))


# In[ ]:


if snakemake.wildcards.chain != "nan":
    pdb_seq = pdb_seqs[df_todo.chain[0]]
    template_seq = templates_seqs[df_todo.template[0]]
    pdb_seq_new = list(pdb_seq)


# In[ ]:


if snakemake.wildcards.chain != "nan":
    log = open(snakemake.log[0],"w")
    seqsd = f">PDB|{df_todo.chain[0]}\n{pdb_seq}\n>{df_todo.template[0]}\n{template_seq}\n"
    foundmis = False
    for m in str(df_todo.mutations[0]).split(","):  
        if m.find("ins") == -1 and m != 'nan' :
            wts = m[0]
            muts = m[len(m)-1]
            pospdb = int(m[1:len(m)-1])
            postempl = PDB_TEMPLATE[pospdb]
            pdb_symb = pdb_seq[pospdb-1]
            templ_symb = template_seq[postempl-1]
            pdb_seq_new[pospdb-1] = muts
            warnm = ""
            if wts != pdb_symb:
                warnm += "!MISMATCH_IN_PDB!"
                foundmis = True
            if wts != templ_symb:
                warnm += " !MISMATCH_IN_TEMPLATE!"
                foundmismis = True
            log.write(f"WT:{wts} MUT:{muts} POSINPDB:{pospdb} POSINTEMPL:{postempl} PDBSYMB:{pdb_symb} TEMPLSYMB:{templ_symb} {warnm}\n")
    if foundmis:
        log.write(seqsd)
    log.close()


# In[ ]:


if snakemake.wildcards.chain != "nan":
    pdb_seq_new4out = ''.join(pdb_seq_new).replace('-','')
    with open(snakemake.output[0],"w") as of:
        chain = df_todo.chain[0]
        of.write(f">PDB|{chain}\n{pdb_seq_new4out}\n")
        for chn in pdb_seqs:
            if chn != chain:
                seqold = pdb_seqs[chn]
                of.write(f">PDB|{chn}\n{seqold}\n")
else:
    with open(snakemake.output[0],"w") as of:
        for chn in pdb_seqs:
            seqold = pdb_seqs[chn]
            of.write(f">PDB|{chn}\n{seqold}\n")

