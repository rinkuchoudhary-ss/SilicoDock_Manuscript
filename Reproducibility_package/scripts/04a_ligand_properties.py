import os, pandas as pd, numpy as np
from math import sqrt
from scipy.stats import spearmanr
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski
RDLogger.DisableLog('rdApp.*')
R=r"C:\Users\ADMIN\Desktop\SilicoScientia\AI_Docking_tool"; B=R+r"\Manuscript_Revision\benchmarking data"
def props(m):
    ri=m.GetRingInfo(); sizes=[len(r) for r in ri.AtomRings()]
    # fused ring systems
    rings=[set(r) for r in ri.AtomRings()]; systems=[]
    for r in rings:
        hit=[s for s in systems if len(s&r)>=2]
        for s in hit: systems.remove(s); r=r|s
        systems.append(r)
    return dict(ha=m.GetNumHeavyAtoms(), rotb=Lipinski.NumRotatableBonds(m), rings=len(sizes), maxring=max(sizes) if sizes else 0,
                macro=int(any(s>=12 for s in sizes)), ringsys=len(systems), fsp3=rdMolDescriptors.CalcFractionCSP3(m), stereo=len(Chem.FindMolChiralCenters(m,includeUnassigned=True)))
rows=[]
cr=pd.read_csv(B+r"\casf_new_pose_rmsd.csv"); cr['c']=cr.entry_id.str.lower()
for c in cr.c.unique():
    m=Chem.MolFromMolFile(f"{R}\\RigiDiff\\data\\All_P-L\\{c}\\{c}_ligand.sdf",removeHs=True)
    if m is None: m=Chem.MolFromMol2File(f"{R}\\RigiDiff\\data\\All_P-L\\{c}\\{c}_ligand.mol2",removeHs=True)
    rows.append(dict(set='CASF',c=c,**(props(m) if m else {})))
pb=pd.read_csv(B+r"\posebuster_pose_rmsd.csv")
for c in pb.entry_id.unique():
    m=Chem.MolFromMolFile(f"{B}\\posebusters_results\\posebusters_results\\{c}\\pose_01.sdf",removeHs=True)
    rows.append(dict(set='PB',c=c,**(props(m) if m else {})))
P=pd.DataFrame(rows); print(P.groupby('set').size().to_dict(), 'missing props', P.ha.isna().sum())
P.to_csv(r"C:\Users\ADMIN\AppData\Local\Temp\rdk\ligprops.csv",index=False)
print(P.groupby('set')[['ha','rotb','rings','maxring','macro','ringsys']].describe().T.loc[(slice(None),['min','50%','max']),:].round(1).to_string())
print('macrocycles:',P.groupby('set').macro.sum().to_dict())
