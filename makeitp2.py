import math
from pathlib import Path
import re
import os

def get_dis(pos1, pos2) :
    x1,y1,z1=[float(pos1[i]) for i in range(3)]
    x2,y2,z2=[float(pos2[i]) for i in range(3)]
    dis=((x1-x2)**2 + (y1-y2)**2 +(z1-z2)**2)**0.5
    return dis
#print(get_dis([11.729,10.896,0.180],[11.701,10.851,0.095]))

class gro_atoms(object):
    def __init__(self,line:str):
        self.index=[(0,8),(8,15),(15,20),(20,28),(28,36),(36,44)]
        self.allinfo=[line[s:e] for s,e in self.index]
        self.residue_name=self.allinfo[0]
        self.name=self.allinfo[1]
        self.index=int(self.allinfo[2])
        self.pos=[10*float(self.allinfo[i]) for i in [3,4,5]]
        #乘10将nm转为A
#test1=gro_atoms("  247MET      N    1  11.729  10.896   0.180")
#test2=gro_atoms("  247MET     H1    2  11.701  10.851   0.095")
#print(get_dis(test1.pos,test2.pos))
#ans=1.00169855745129

topol=Path("topol.top")
if topol.is_file():
    print("Find the topol.top of system")
else:
    print("No topol.top found\nPlease check your input")
    exit()
topol_itp=[item for item in map(str,list(Path('.').glob('topol*itp')))]
posre_itp=[item for item in map(str,list(Path('.').glob('posre_Protein*itp')))]
if len(topol_itp) < 1 or len(posre_itp) < 1:
    print("Sorry, protein's topol files don't exist in local dir.\nMaybe all file is in the topol.top, which cannot be handled.\nPlease check your input files.")
    exit()
elif len(topol_itp) != len(posre_itp):
    print(f"{len(topol_itp)} chain while {len(posre_itp)} porse\n There is something wrong")
    exit()
else:
    print(f"{len(topol_itp)} pair top file of protein in current dir")
lig_name=[item for item in map(str,list(Path('.').glob('*prm')))]
lig_name=lig_name[0].split(".")[0].split("/")[-1].upper()
lig_name_lower=lig_name[0].split(".")[0].split("/")[-1]
lig_top=[item for item in map(str,list(Path('/home/chengyj/git_handbook/makeitp/test/5q0n_9L4c').glob(f'{lig_name_lower}*')))]
if len(lig_top) != 5:
    print("This script needs 5 ligand topol file,\n which are .prm, .itp, .top, _ini.pdb and .gro\n There is not enough file in current dir")
    exit()
else:
    print(f"{len(lig_top)} of ligand topol files in current dir")


topol_order=[]
with open ('topol.top') as topol_all_file:
    for line in topol_all_file:
        if "topol_Protein_chain" in line:
            topol_order.append(line.split()[1][1:-1]) 

end=[]
for item in topol_order:
    lines_index=[""]*2
    with open (item) as topol_chain_file: 
        for line in topol_chain_file:
            if line == "[ bonds ]\n":
                data_re=[lines_index[1].split()[i] for i in [0,2,3]]
                data_re.append(re.sub('topol_Protein','posre',item))
                end.append(data_re)
                break  
            lines_index[0],lines_index[1] = line,lines_index[0]  

pro_atoms_of_gro=[]
lig_atoms_of_gro=[]
with open('solv_ions.gro') as file_of_gro:
    for line in file_of_gro:
        if line[0]==' ':
            if 'SOL' in line:
                break
            if lig_name in line:
                lig_atoms_of_gro.append(gro_atoms(line))
            else:
                pro_atoms_of_gro.append(gro_atoms(line))



split_chain=[(0,""),]
for item in end:
    split_chain.append((split_chain[-1][0]+int(item[0]),item[-1]))

restain_atoms=list(pro_atoms_of_gro)
for lig_atoms in lig_atoms_of_gro:
    del_count=0
    for res_atom_index in range(len(restain_atoms)):
        #去除12A以内及H原子
        if get_dis(lig_atoms.pos,restain_atoms[res_atom_index-del_count].pos) < 12 or 'H' in restain_atoms[res_atom_index-del_count].name:
            restain_atoms.pop(res_atom_index-del_count)
            del_count+=1

count=0
for index_chain in range(len(split_chain)-1):
    with open(split_chain[index_chain+1][1],'w') as prose_itp:
        prose_itp.write("[ position_restraints ]\n")
        for atom in restain_atoms:
            if atom.index > split_chain[index_chain+1][0]:
                break
            if split_chain[index_chain][0] < atom.index <= split_chain[index_chain+1][0]:
                prose_itp.write(f"{str(atom.index-split_chain[index_chain][0]).rjust(5,' ')}    1       1000       1000       1000\n")
                count+=1
if count == len(restain_atoms):
    print("All itp files have been written")
else:
    print("There is something wrong, some restain atoms are missing")

for chain_topol in topol_order:
    os.system(f"sed -i \"s/posre_Protein_chain/posre_chain/\" {chain_topol}")
print("The topol of every chain is modified.\nNow everything has finished")
