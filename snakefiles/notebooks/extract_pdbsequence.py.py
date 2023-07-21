#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import warnings
warnings.filterwarnings('ignore')
from Bio import PDB, SeqIO, pairwise2
import pandas as pd
sys.path.append(snakemake.params.path)
from pdbio.pdbfile import PDBFile
import sys, warnings
import sys
from Bio.Data.IUPACData import protein_letters_3to1_extended
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.pairwise2 import format_alignment
from Bio.Seq import MutableSeq
from Bio.Align import substitution_matrices
logf = open(snakemake.log[0],"w")


# In[ ]:


#prepare mappings
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

#kmain


# In[ ]:


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
    #format_string = "{:<6}{:>5}  {:<3}{:1}{:<3} {:1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:>2}"
    
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



# for l in PDBFile(snakemake.input[0]).content:
#     if l[0:5].find("ATOM") > -1 or l[0:5].find("HETA") > -1:  
#         f = parse_atom_line(l)
#         l=l.replace("\n","")
#         l2 = generate_atom_line(f)
#         if l != l2:
#                 print(f"|{l}|")
#                 print(f"|{l2}|")
#                 print("--")


# In[ ]:


# read in PDB
pdb = PDBFile(snakemake.input[0])
# get atom content
Afields = []
for l in pdb.content:
    startlets = l[0:5]
    if startlets.find("ATOM") > -1 or startlets.find("HETA") > -1 :#HETA
        Afields.append(parse_atom_line(l))
    


# In[ ]:


#collect hetero residues
hets = []
standard3As = [s.upper() for s in kext]
for a in Afields:
    if a["record_name"] != "ATOM":
        hets.append(a["res_name"])
hets = list(set([i.upper() for i in hets]))
hets = list(set(hets).difference(set(standard3As)))
logf.write(f"DTECTED HETS: {hets}\n")
#print(f"DTECTED HETS: {hets}\n")
# remove AA seymbols from hets # in some structures like 4CPA some AA are trated as HETERO residues (GLY in 4CPA case)


# In[ ]:


matrix = substitution_matrices.load("BLOSUM62")
SEQRES_updated = ""
hets4unchain = [] # a tuple of resname and chain - for removal of heteroatoms
for chain in pdb:
    seqres_seq, seqreq_seq_3L = chain.sequence_seqres(return_3L=True)
    #print(seqres_seq, seqreq_seq_3L)
    atom_seq = chain.sequence_atom()
    
    #check how many to keep at the ends
    align = pairwise2.align.globalds(seqres_seq, atom_seq, matrix, -11, -1, penalize_end_gaps=(False, False),one_alignment_only=True)[0]
    sequence = MutableSeq(align[0])
    structure = MutableSeq(align[1])
#     print(f"CHAIN {chain.name}")
#     print(seqres_seq)
#     print(atom_seq)
#     print("---")
#     print(sequence)
#     print("---")
#     print(structure)
    #get position at the begining at seqres protruding no more than 1 residue at the N- terminus
    lenaln = len(sequence)
    posST = 0 # position structure
    posSE = 0 # position sequence
    startpos = 1
    for ali in range(0,lenaln):
        if sequence[ali]=='X':
            sequence[ali]='-'
        if structure[ali]=='X':
            structure[ali]='-'
    for ali in range(0,lenaln):
        if sequence[ali] != '-':
            posSE += 1
        if structure[ali] != '-':
            posST += 1
        if posST == 1:
            startpos = posSE
            if startpos - 1  >= 1:
                startpos = startpos - 1
    seqlen_woG = len(str(sequence).replace('-',''))
    strlen_woG = len(str(structure).replace('-',''))



    # get the position at the end of seqres protruding no more than one residue at the C- terminus
    endpos = seqlen_woG
    posST = 0 # position structure
    posSE = 0 # position sequence
    for ali in range(0,lenaln):
        if sequence[ali] != '-':
            posSE += 1
        if structure[ali] != '-':
            posST += 1
        if posST == strlen_woG:
            endpos = posSE
            if posSE + 1 <= seqlen_woG :
                endpos = endpos + 1
            break
#     print(posSE,posST, strlen_woG, endpos)
        
                
                
    #check for hetero in sequences
    ct = 0
    ctseqres = 0
    ctall = 0
    ctaasymb = 0
    out_laines = []
    out = []
    for s in seqres_seq:
        ctall += 1
        var3L = seqreq_seq_3L[ctall-1].upper() 
        if s != 'X' and (not (var3L in hets)):
            ctaasymb += 1
            if ctaasymb >= startpos and ctaasymb <= endpos :
                ctseqres += 1
                ct += 1
                out.append(var3L)
                #print(var3L, ct)
                if ct  == 13:
                    outl = " ".join(out)
                    out_laines.append(outl)
                    #print("AAA",outl)
                    out=[]
                    ct = 0
        else:
            hets4unchain.append((var3L,chain.name))
            
            logf.write(f"removing {var3L} residues in chain {chain.name}\n") 
            #print(f"removing {var3L} residues in chain {chain.name}\n", s)
        #print(out_laines)
        
    outl = " ".join(out)
    out_laines.append(outl)
    out_laines_with_prefixes = []
    ct = 0
    ctl = 0
    if abs(ctseqres-strlen_woG) > 2:
                    logf.write(f"WARNING updated SEQRES length {ctseqres} and STRUCTURE length {strlen_woG} for chain {chain.name}\n") 
    for l in out_laines:
        if len(l) > 0:
            ctl += 1
            ll = "SEQRES "+str(ctl).rjust(3)+chain.name.rjust(2)+" "+str(ctseqres).rjust(4)+"  "+l
            SEQRES_updated += ll+"\n"
hets4unchain = list(set(hets4unchain))


# In[ ]:


#modify atom content
for i in range(0,len(Afields)):
    fields = Afields[i]
    chain_id = fields["chain_id"]
    res_name = fields["res_name"]
    #mark hetatoms that shoulc be removed - that was part of seqres
    if (res_name, chain_id) in hets4unchain:
        fields['res_name'] = 'BAAAD'
    #remove chain assaignments for  all hetatoms
    if fields["record_name"] == 'HETATM':
        fields["chain_id"] = ' ' 
ATOM_lines = []
for i in range(0,len(Afields)):
    fields = Afields[i]
    l = generate_atom_line(fields)
    if fields['res_name'] != 'BAAAD':
        ATOM_lines.append(l)


# In[ ]:


#regenerate all
out = ""
ATOMSPUT = False
SEQRESPUT = False
#SEQRES_updated
for l in pdb.content:
    if l[0:5].find("ATOM") > -1 or l[0:5].find("HETA") > -1:
        if not ATOMSPUT:
            out += "\n".join(ATOM_lines)+"\n"
            ATOMSPUT = True
    elif l[0:7].find("SEQRES ") > -1:
        if not SEQRESPUT:
            out += SEQRES_updated
            SEQRESPUT = True
    else:
        out += l
with open(snakemake.output.pdb,"w") as outF:
    outF.write(out)


# In[ ]:


pdb = PDBFile(snakemake.output.pdb)
out=[]
for chain in pdb:
    idch = ">chain"+chain.name
    out.append(idch)
    atom_seq = chain.sequence_atom()
    seqres_seq = chain.sequence_seqres()
    if not seqres_seq == None:
        out.append(seqres_seq)
    else:
        out.append(atom_seq)
        
with open(snakemake.output.fasta,"w") as fo:
    fo.write("\n".join(out))

