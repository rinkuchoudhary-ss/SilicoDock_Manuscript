import numpy as np, pandas as pd, os, warnings
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*'); warnings.filterwarnings('ignore')
R=r"C:\Users\ADMIN\Desktop\SilicoScientia\AI_Docking_tool"; A=R+r"\RigiDiff\data\All_P-L"
B=R+r"\Manuscript_Revision\benchmarking data"
ids=[str(x).lower() for x in pd.read_csv(B+r"\casf_new_pose_rmsd.csv").entry_id.unique()]
def prot_xyz(path):
    xyz=[]
    for l in open(path,errors='ignore'):
        if l.startswith('ATOM') and l[76:78].strip()!='H' and not l[12:16].strip().startswith('H'):
            try: xyz.append((float(l[30:38]),float(l[38:46]),float(l[46:54])))
            except ValueError: pass
    return np.array(xyz)
rows=[]
for i,c in enumerate(ids):
    m=Chem.MolFromMolFile(f"{A}\\{c}\\{c}_ligand.sdf",removeHs=True)
    if m is None or m.GetNumConformers()==0: continue
    lx=m.GetConformer().GetPositions(); px=prot_xyz(f"{A}\\{c}\\{c}_protein.pdb")
    dmin=np.min(np.linalg.norm(px[:,None,:]-lx[None,:,:],axis=2)) if len(px)*len(lx)<6e7 else np.nan
    # internal: nonbonded pairs with topological distance >=3
    tm=Chem.GetDistanceMatrix(m); D=np.linalg.norm(lx[:,None,:]-lx[None,:,:],axis=2); iu=np.triu_indices(len(lx),1); sel=tm[iu]>=3
    intmin=D[iu][sel].min() if sel.any() else np.nan
    # bond deviation crystal vs ETKDG+MMFF conformer
    dev=np.nan
    try:
        mh=Chem.AddHs(m); ps=AllChem.ETKDGv3(); ps.randomSeed=42
        if AllChem.EmbedMolecule(mh,ps)==0:
            try: AllChem.MMFFOptimizeMolecule(mh,maxIters=500)
            except Exception: pass
            cm=Chem.RemoveHs(mh).GetConformer().GetPositions(); ds=[]
            for b in m.GetBonds():
                a1,a2=b.GetBeginAtomIdx(),b.GetEndAtomIdx(); l1=np.linalg.norm(lx[a1]-lx[a2]); l2=np.linalg.norm(cm[a1]-cm[a2]); ds.append(abs(l1-l2)/l2)
            dev=max(ds) if ds else 0
    except Exception: pass
    rows.append(dict(c=c,prot_min=dmin,int_min=intmin,bond_dev=dev))
    if i%60==0: print(i,flush=True)
d=pd.DataFrame(rows); d.to_csv(r"C:\Users\ADMIN\AppData\Local\Temp\rdk\crystal_anchors.csv",index=False)
print(len(d))
for col,lab in (('prot_min','closest ligand-protein heavy-atom distance in the crystal complex (A)'),('int_min','closest internal nonbonded heavy-atom distance in the crystal ligand (A)'),('bond_dev','largest relative bond-length difference, crystal vs ETKDG+MMFF conformer')):
    x=d[col].dropna(); print(lab,'n',len(x),'min %.2f p1 %.2f p5 %.2f median %.2f p95 %.2f p99 %.2f max %.2f'%(x.min(),x.quantile(.01),x.quantile(.05),x.median(),x.quantile(.95),x.quantile(.99),x.max()))
print('share of crystals with prot_min >= 1.8:',(d.prot_min>=1.8).mean(),' >=2.0:',(d.prot_min>=2.0).mean())
print('share int_min >= 1.0:',(d.int_min>=1.0).mean())
print('share bond_dev <= 0.15:',(d.bond_dev.dropna()<=0.15).mean(),' <=0.10:',(d.bond_dev.dropna()<=0.10).mean(),' <=0.25:',(d.bond_dev.dropna()<=0.25).mean())
