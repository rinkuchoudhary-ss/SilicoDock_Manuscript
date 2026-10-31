import json, os, pickle, numpy as np, pandas as pd
from collections import defaultdict
from Bio.Align import PairwiseAligner, substitution_matrices
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
R = r"C:\Users\ADMIN\Desktop\SilicoScientia\AI_Docking_tool"; A = R + r"\RigiDiff\data\All_P-L"; B = R + r"\Manuscript_Revision\benchmarking data"
OUT = r"C:\Users\ADMIN\AppData\Local\Temp\rdk"
T3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','HID':'H','HIE':'H','HIP':'H','CYX':'C','MSE':'M'}
def pocket(prot_path, lig_xyz, cutoff=6.0):
    res = defaultdict(list); names = {}
    for l in open(prot_path, errors='ignore'):
        if not l.startswith('ATOM'): continue
        if l[76:78].strip() == 'H' or l[12:16].strip().startswith('H'): continue
        try: xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        except ValueError: continue
        k = (l[21], int(l[22:26]) if l[22:26].strip().lstrip('-').isdigit() else 0, l[26]); res[k].append(xyz); names[k] = T3.get(l[17:20].strip(), 'X')
    keep = []
    for k, a in res.items():
        a = np.array(a)
        if np.min(np.linalg.norm(a[:, None, :] - lig_xyz[None, :, :], axis=2)) <= cutoff: keep.append(k)
    keep.sort()
    return ''.join(names[k] for k in keep)
split = json.load(open(R + r"\RigiDiff\data\test_split.json")); train = [i.lower() for i in split['train']]
cache = OUT + r"\pocket_train.pkl"
if os.path.exists(cache): ref = pickle.load(open(cache, 'rb'))
else:
    ref = {}
    for n, i in enumerate(train):
        d = f"{A}\\{i}"
        m = Chem.MolFromMolFile(f"{d}\\{i}_ligand.sdf", removeHs=True)
        if m is None or m.GetNumConformers() == 0: continue
        try: s = pocket(f"{d}\\{i}_protein.pdb", m.GetConformer().GetPositions())
        except Exception: continue
        if len(s) >= 5: ref[i] = s
        if n % 3000 == 0: print('ref', n, flush=True)
    pickle.dump(ref, open(cache, 'wb'))
print('training pockets', len(ref), 'median length', np.median([len(s) for s in ref.values()]), flush=True)
ids = list(ref); seqs = [ref[i] for i in ids]; K = 2
inv = defaultdict(list)
for k, s in enumerate(seqs):
    for km in {s[j:j+K] for j in range(len(s) - K + 1)}: inv[km].append(k)
inv = {k: np.array(v, dtype=np.int32) for k, v in inv.items()}
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load('BLOSUM62'); al.open_gap_score = -10; al.extend_gap_score = -1; al.mode = 'global'
def best(q, selfid):
    cnt = np.zeros(len(seqs), dtype=np.int32)
    for km in {q[j:j+K] for j in range(len(q) - K + 1)}:
        a = inv.get(km)
        if a is not None: cnt[a] += 1
    order = np.argsort(-cnt)[:60]; b = 0.0
    for k in order:
        if ids[k] == selfid: continue
        a = al.align(q, seqs[k])[0]; b = max(b, a.counts().identities / len(q))
    return b
rows = []
cr = pd.read_csv(B + r"\casf_new_pose_rmsd.csv"); casf = sorted({x.lower() for x in cr.entry_id.unique()})
for n, c in enumerate(casf):
    m = Chem.MolFromMolFile(f"{A}\\{c}\\{c}_ligand.sdf", removeHs=True); q = pocket(f"{A}\\{c}\\{c}_protein.pdb", m.GetConformer().GetPositions())
    rows.append(dict(id=c, set='CASF', pocket_len=len(q), pocket_id=best(q, c)))
pb = pd.read_csv(B + r"\posebuster_pose_rmsd.csv"); PBD = B + r"\posebusters_results\posebusters_results"
for n, c in enumerate(pb.entry_id.unique()):
    m = Chem.MolFromMolFile(f"{PBD}\\{c}\\pose_01.sdf", removeHs=True); q = pocket(f"{PBD}\\{c}\\{c}_protein_prepared.pdb", m.GetConformer().GetPositions())
    rows.append(dict(id=c, set='PoseBusters', pocket_len=len(q), pocket_id=best(q, c)))
    if n % 100 == 0: print('pb', n, flush=True)
df = pd.DataFrame(rows); df.to_csv(OUT + r"\pocket_sim.csv", index=False)
print(df.groupby('set').pocket_id.describe().round(3).to_string())
