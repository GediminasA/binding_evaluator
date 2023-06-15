#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import json


# In[ ]:


df = pd.read_csv(snakemake.input.coordinates, sep = "\t", names = "qseqid sseqid qstart qend sstart send length pident qseq sseq".split())
df = df.query("sseqid == @snakemake.wildcards.seqtempl")
df = df.sort_values(by=["length"], ascending= False)
data = dict(df.iloc[0])
data


# In[ ]:


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


df = pd.DataFrame.from_dict({"PDB":qpos_seqs,"TEMPLATE":spos_seqs})
df["CHAIN"] = data["qseqid"].split("|")[1]
df.to_csv(snakemake.output[0], index = False)

