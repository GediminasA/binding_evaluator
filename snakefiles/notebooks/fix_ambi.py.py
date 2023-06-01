#!/usr/bin/env python
# coding: utf-8

# In[2]:


from Bio import PDB, SeqIO, pairwise2
import pandas as pd
sys.path.append(snakemake.params.path)
from pdbio.pdbfile import PDBFile
import sys, warnings
import sys
from Bio.Data.IUPACData import protein_letters_3to1_extended
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.pairwise2 import format_alignment



# In[3]:


kext = list(protein_letters_3to1_extended.keys())
kmain = list(protein_letters_3to1.keys())
spec3a = list(set(kext)- set(kmain))
spec1a = [protein_letters_3to1_extended[a] for a in spec3a]
spec3a
map1to3 = {}
map3to1 = {}
for k in protein_letters_3to1_extended:
    map1to3[protein_letters_3to1_extended[k]] = k.upper()
    map3to1[k.upper()] = protein_letters_3to1_extended[k]

    


# In[4]:


df = pd.read_csv(snakemake.input.swiss, sep="\t",
                 header=0,names=["E","Pid","START","END","CHAIN","SID","QSEQ","SSEQ"],)
df.head(n=100)


# In[5]:


# chack what needs to be masked
map2n = {}
ct = 0
for ch in df.CHAIN:
    ct += 1
    chid = ch.replace("chain","")
    lid = 0
    for l in df.QSEQ[ct-1]:
        lid += 1
        ql = df.QSEQ[ct-1][lid-1]
        sl = df.SSEQ[ct-1][lid-1]
        if ql in spec1a:
            kk = (chid,lid,ql)
            if not kk in map2n.keys():
                map2n[kk] = sl
map2n    


# In[6]:


pdb = PDBFile(snakemake.input[1])


# In[7]:


# for chain in pdb:
#     seqres_seq = chain.sequence_seqres()
#     atom_seq = chain.sequence_atom()
#     for a in pairwise2.align.globalxx(seqres_seq, atom_seq,one_alignment_only=True):
#         print(format_alignment(*a))
#         print("SEQRES",seqres_seq)
#         print("ATOMRE",atom_seq)
    
    


# In[8]:


SEQRES_updated = ""
for chain in pdb:
    seqres_seq = chain.sequence_seqres()
    ct = 0
    ctchain = 0
    
    out_laines = []
    out = []
    for s in seqres_seq:
        ct += 1
        ctchain += 1
        var3L = map1to3[s]
        #print(s,chain.name,ctchain,var3L)
        if var3L.find("X") > -1:
            id4sub = (chain.name,ctchain,s)
            #print(id4sub)
            #print(var3L,s,ct)
            #print(id4sub,"id")
            if id4sub in map2n.keys():
                su = map2n[id4sub]
                var3L = map1to3[su]
            else:
                raise  Exception(f"unknown proper value for residue {s}, chain {chain.name}, number {ct}")            
        out.append(var3L)
        if ct  == 13:
            outl = " ".join(out)
            out_laines.append(outl)
            out=[]
            ct = 0
    outl = " ".join(out)
    out_laines.append(outl)
    out_laines_with_prefixes = []
    ct = 0
    nch = len(seqres_seq)
    ctl = 0
    for l in out_laines:
        if len(l) > 0:
            ctl += 1
            ll = "SEQRES "+str(ctl).rjust(3)+chain.name.rjust(2)+" "+str(nch).rjust(4)+"  "+l
            SEQRES_updated += ll+"\n"
#print(SEQRES_updated)


# In[9]:


map3D = {}
for chain in pdb:
    seqres_seq = chain.sequence_seqres()
    atom_seq = chain.sequence_atom()
    for a in pairwise2.align.globalxx(seqres_seq, atom_seq,one_alignment_only=True):
        outs = format_alignment(*a)
        ctA = 0
        ctB = 0 
        for i,lA in enumerate(a.seqA):
            if lA != "-":
                ctA += 1
            lB = a.seqB[i]
            if lB != "-":
                ctB += 1
            lA3 = map1to3[lA]
            var3L = ""
            if lA3.find("X") > -1:
                id4sub = (chain.name,ctA,lA)
                if id4sub in map2n.keys():
                    su = map2n[id4sub]
                    var3L = map1to3[su]
                else:
                    raise  Exception(f"unknown proper value for residue {s}, chain {chain.name}, number {ct}")
            if lB != "-":
                lB3 = map1to3[lB]
                if lB3.find("X") > -1:
                    id4sub = ctB
                    map3D[(chain.name,id4sub,lB3)] = var3L
                    if var3L == "":
                        raise Exception(f"Inconsistencies between seqres and sequence records, residue {lB3}, chain {chain.name}, number {id4sub}")
            
print(map3D)
    


# In[10]:


def parse_atom_line(line):
    """
    Parse an ATOM/HETATM line from a PDB file and return a dictionary of the fields.

    Parameters
    ----------
    line : str
        The ATOM/HETATM line to parse.

    Returns
    -------
    dict
        A dictionary containing the parsed fields, with the keys being the field names.
    """
    fields = {
        'record_name': line[0:6].strip(),
        'atom_number': int(line[6:11]),
        'atom_name': line[12:16].strip(),
        'alt_loc': line[16],
        'res_name': line[17:20].strip(),
        'chain_id': line[21],
        'res_seq': int(line[22:26]),
        'i_code': line[26],
        'x': float(line[30:38]),
        'y': float(line[38:46]),
        'z': float(line[46:54]),
        'occupancy': float(line[54:60]),
        'temp_factor': float(line[60:66]),
        'element': line[76:78].strip(),
        'charge': line[78:80].strip()
    }
    return fields






# In[11]:


def generate_atom_line(fields):
    """
    Generate an ATOM/HETATM line for a PDB file from a dictionary of the fields.

    Parameters
    ----------
    fields : dict
        A dictionary containing the fields for the ATOM/HETATM line, with the keys being the field names.

    Returns
    -------
    str
        The generated ATOM/HETATM line.
    """
    format_string = "{:<6}{:>5} {:^4}{:1}{:<3} {:1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:>2}"
    if len(fields['atom_name']) == 1:
        format_string = "{:<6}{:>5} {:^4}{:1}{:<3} {:1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:>2}"
    if len(fields['atom_name']) == 3:
        format_string = "{:<6}{:>5} {:>4}{:1}{:<3} {:1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:>2}"
    args = (
        fields['record_name'],
        fields['atom_number'],
        fields['atom_name'],
        fields['alt_loc'],
        fields['res_name'],
        fields['chain_id'],
        fields['res_seq'],
        fields['i_code'],
        fields['x'],
        fields['y'],
        fields['z'],
        fields['occupancy'],
        fields['temp_factor'],
        fields['element'].strip(),
        fields['charge'].strip()
    )
    return format_string.format(*args)


# # Example usage
# atom_line_1 = 'ATOM     52  OE2BGLU A   3      10.301   2.252  11.187  0.50  5.38           O  '
# fields_1 = parse_atom_line(atom_line_1)
# print(fields_1)
# regenerated_line_1 = generate_atom_line(fields_1)
# print(atom_line_1)
# print(regenerated_line_1)
# print(atom_line_1 == regenerated_line_1)

# # #atom_line_2 = 'ATOM      2  CA  ALA A   1       8.710  19.126   0.464  1.00 36.68           C  '
# # #fields_2 = parse_atom_line(atom_line_2)
# # #regenerated_line_2 = generate_atom_line(fields_2)
# # #print(atom_line_2 == regenerated_line_2)


# In[12]:


# get atom content
Afields = []
for l in pdb.content:
    if l[0:5].find("ATOM") > -1:
        Afields.append(parse_atom_line(l))
    


# In[13]:


#modify atom content
resct = 0
chain = "NA"
resini = "NaN"
for i in range(0,len(Afields)):
    fields = Afields[i]
    chain_id = fields["chain_id"]
    res_name = fields["res_name"]
    if chain_id != chain:
        resct = 0
        chain = chain_id
    if resini != res_name:
        resct += 1
        resini = res_name

    if res_name.find("X") > -1:
        k = (chain_id, resct,res_name)
        res_name = map3D[k]
        fields["res_name"] = res_name
        


# In[14]:


# regenerate lines
ATOM_lines = []
for i in range(0,len(Afields)):
    fields = Afields[i]
    l = generate_atom_line(fields)
    ATOM_lines.append(l)
    if l.find("GLX") > -1:
        print(l)


# In[15]:


#regenerate all
out = ""
ATOMSPUT = False
SEQRESPUT = False
#SEQRES_updated
for l in pdb.content:
    if l[0:5].find("ATOM") > -1:
        if not ATOMSPUT:
            out += "\n".join(ATOM_lines)+"\n"
            ATOMSPUT = True
    elif l[0:7].find("SEQRES ") > -1:
        if not SEQRESPUT:
            out += SEQRES_updated
            SEQRESPUT = True
    else:
        out += l
with open(snakemake.output[0],"w") as outF:
    outF.write(out)

