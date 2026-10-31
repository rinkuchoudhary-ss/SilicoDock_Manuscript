import json, os, sys, pickle, numpy as np, pandas as pd
from collections import defaultdict
from Bio.Align import PairwiseAligner, substitution_matrices
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog('rdApp.*')
R = r"C:\Users\ADMIN\Desktop\SilicoScientia\AI_Docking_tool"
A = R + r"\RigiDiff\data\All_P-L"
B = R + r"\Manuscript_Revision\benchmarking data"
OUT = r"C:\Users\ADMIN\AppData\Local\Temp\rdk"
T3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','HID':'H','HIE':'H','HIP':'H','CYX':'C','MSE':'M'}

def read_pdb(path):
    """returns {chain: (seq, coords of all atoms)}"""
    seqs = defaultdict(list); xyz = defaultdict(list); seen = set()
    with open(path, errors='ignore') as f:
        for l in f:
            if not l.startswith('ATOM'): continue
            ch = l[21]; rn = l[17:20].strip(); key = (ch, l[22:27])
            try: xyz[ch].append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            except ValueError: continue
            if l[12:16].strip() == 'CA' and key not in seen:
                seen.add(key); seqs[ch].append(T3.get(rn, 'X'))
    return {c: (''.join(s), np.array(xyz[c])) for c, s in seqs.items() if len(s) >= 20}

def load_lig(path_sdf, path_mol2=None):
    m = Chem.MolFromMolFile(path_sdf, sanitize=True, removeHs=True) if os.path.exists(path_sdf) else None
    if m is None and path_mol2 and os.path.exists(path_mol2):
        m = Chem.MolFromMol2File(path_mol2, sanitize=True, removeHs=True)
    return m

def fp_scaf(m):
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)
    try:
        sc = MurckoScaffold.GetScaffoldForMol(m); ex = Chem.MolToSmiles(sc) if sc.GetNumAtoms() else ''
        gs = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(sc)) if sc.GetNumAtoms() else ''
    except Exception: ex = gs = ''
    return fp, ex, gs

split = json.load(open(R + r"\RigiDiff\data\test_split.json"))
train = [i.lower() for i in split['train']]; trainset = set(train)
casf = pd.read_csv(B + r"\casf_new_pose_rmsd.csv"); casf['id'] = casf.entry_id.str.lower()
pb = pd.read_csv(B + r"\posebuster_pose_rmsd.csv")
casf_ids = sorted(casf.id.unique()); pb_ids = sorted(pb.entry_id.unique())

cache = OUT + r"\ref_train.pkl"
if os.path.exists(cache):
    ref = pickle.load(open(cache, 'rb'))
else:
    ref = {'seq': [], 'seq_owner': [], 'ids': [], 'fp': [], 'ex': [], 'gs': []}
    for n, i in enumerate(train):
        d = f"{A}\\{i}"
        try:
            for ch, (s, _) in read_pdb(f"{d}\\{i}_protein.pdb").items():
                ref['seq'].append(s); ref['seq_owner'].append(i)
        except Exception: pass
        m = load_lig(f"{d}\\{i}_ligand.sdf", f"{d}\\{i}_ligand.mol2")
        if m is not None:
            fp, ex, gs = fp_scaf(m); ref['ids'].append(i); ref['fp'].append(fp); ref['ex'].append(ex); ref['gs'].append(gs)
        if n % 2000 == 0: print('ref', n, flush=True)
    pickle.dump(ref, open(cache, 'wb'))
print('ref chains', len(ref['seq']), 'ligands', len(ref['ids']), flush=True)

# unique sequence -> owners
useq = defaultdict(set)
for s, o in zip(ref['seq'], ref['seq_owner']): useq[s].add(o)
useqs = list(useq); K = 3
inv = defaultdict(list)
for k, s in enumerate(useqs):
    for km in {s[j:j+K] for j in range(len(s) - K + 1)}: inv[km].append(k)
inv = {k: np.array(v, dtype=np.int32) for k, v in inv.items()}
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load('BLOSUM62'); al.open_gap_score = -10; al.extend_gap_score = -0.5; al.mode = 'global'

def best_identity(q, selfid):
    cnt = np.zeros(len(useqs), dtype=np.int32)
    for km in {q[j:j+K] for j in range(len(q) - K + 1)}:
        a = inv.get(km)
        if a is not None: cnt[a] += 1
    order = np.argsort(-cnt)[:25]; best = 0.0
    for k in order:
        if cnt[k] == 0: break
        if useq[useqs[k]] <= {selfid}: continue          # only the query itself has this sequence
        a = al.align(q, useqs[k])[0]
        ident = a.counts().identities / len(q)
        best = max(best, ident)
    return best

def analyse(cid, prot, lig):
    ch = read_pdb(prot); lx = lig.GetConformer().GetPositions() if lig.GetNumConformers() else None
    contact = []
    for c, (s, xyz) in ch.items():
        if lx is None: contact.append(c); continue
        dmin = np.min(np.linalg.norm(xyz[:, None, :] - lx[None, :, :], axis=2)) if len(xyz) * len(lx) < 4e7 else 99
        if dmin <= 6.0: contact.append(c)
    if not contact: contact = list(ch)
    ids = [best_identity(ch[c][0], cid) for c in contact if c in ch]
    fp, ex, gs = fp_scaf(lig)
    sims = np.array(DataStructs.BulkTanimotoSimilarity(fp, ref['fp']))
    for k, o in enumerate(ref['ids']):
        if o == cid: sims[k] = -1
    exs = ex != '' and any(e == ex for o, e in zip(ref['ids'], ref['ex']) if o != cid)
    gss = gs != '' and any(e == gs for o, e in zip(ref['ids'], ref['gs']) if o != cid)
    return dict(id=cid, seq_id=max(ids) if ids else np.nan, tan=float(sims.max()), scaf_exact=bool(exs), scaf_generic=bool(gss))

rows = []
for n, cid in enumerate(casf_ids):
    d = f"{A}\\{cid}"; lig = load_lig(f"{d}\\{cid}_ligand.sdf", f"{d}\\{cid}_ligand.mol2")
    if lig is None: rows.append(dict(id=cid, set='CASF')); continue
    r = analyse(cid, f"{d}\\{cid}_protein.pdb", lig); r['set'] = 'CASF'; rows.append(r)
    if n % 50 == 0: print('casf', n, flush=True)
PBD = B + r"\posebusters_results\posebusters_results"
for n, cid in enumerate(pb_ids):
    d = f"{PBD}\\{cid}"; prot = f"{d}\\{cid}_protein_prepared.pdb"
    lig = None
    for f in ('pose_01.sdf', 'best_pose.sdf'):
        if os.path.exists(f"{d}\\{f}"):
            lig = Chem.MolFromMolFile(f"{d}\\{f}", removeHs=True)
            if lig is not None: break
    if lig is None or not os.path.exists(prot): rows.append(dict(id=cid, set='PoseBusters')); continue
    r = analyse(cid, prot, lig); r['set'] = 'PoseBusters'; rows.append(r)
    if n % 50 == 0: print('pb', n, flush=True)
df = pd.DataFrame(rows)
# success
pb['id'] = pb.entry_id
best = pd.concat([casf.groupby('id').rmsd_ang.min().rename('best').reset_index().assign(set='CASF'),
                  pb.groupby('id').rmsd_ang.min().rename('best').reset_index().assign(set='PoseBusters')])
df = df.merge(best, on=['id', 'set'], how='left'); df['ok'] = df.best <= 2.0
df['grp'] = df.id.map({**{i: 'train' for i in split['train']}, **{i: 'val' for i in split['val']}, **{i: 'test' for i in split['test']}}).where(df.set == 'CASF')
df.to_csv(OUT + r"\sim_split.csv", index=False)
print(df.groupby('set').agg(n=('id', 'size'), have=('seq_id', lambda s: s.notna().sum()), ok=('ok', 'sum')))
