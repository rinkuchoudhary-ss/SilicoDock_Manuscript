import pandas as pd, numpy as np, os, glob, json
from math import sqrt
B=r"C:\Users\ADMIN\Desktop\SilicoScientia\AI_Docking_tool\Manuscript_Revision\benchmarking data"
def wil(k,n,z=1.96):
    p=k/n; d=1+z*z/n; ce=(p+z*z/(2*n))/d; h=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return 100*(ce-h),100*(ce+h)
A=pd.read_csv(r"C:\Users\ADMIN\AppData\Local\Temp\rdk\crystal_anchors.csv")
anch={}
for col,thr,fn in (('prot_min',1.8,lambda x:(x>=1.8).mean()),('int_min',1.0,lambda x:(x>=1.0).mean()),('bond_dev',0.15,lambda x:(x<=0.15).mean())):
    x=A[col].dropna(); s=100 if col!='bond_dev' else 100
    f=(1 if col!='bond_dev' else 100)
    anch[col]=dict(n=len(x),mn=x.min()*f,p1=x.quantile(.01)*f,med=x.median()*f,p95=x.quantile(.95)*f,mx=x.max()*f,share=100*fn(x))
sens={}
for name,root,rm in (('PB',B+r"\posebusters_results\posebusters_results",B+r"\posebuster_pose_rmsd.csv"),('CASF',B+r"\casf_results_new\casf_results",B+r"\casf_new_pose_rmsd.csv")):
    fs=glob.glob(root+r"\*\docking_summary.csv")
    d=pd.concat([pd.read_csv(f).assign(c=os.path.basename(os.path.dirname(f))) for f in fs]); r=pd.read_csv(rm)
    d['key']=d.c.str.lower()+'|'+d.filename; r['key']=r.entry_id.str.lower()+'|'+r.pose_file
    m=d.merge(r[['key','rmsd_ang']].rename(columns={'rmsd_ang':'rm'}),on='key',how='left'); assert m.rm.notna().all()
    m['chem_ok']=m.chemistry.eq('OK'); m['ok']=m.rm<=2; n_c=m.c.nunique(); out={}
    for lab,mask in (('No receptor-distance or chemistry filter',m.min_dist_ang>=0),('Minimum receptor distance \u2265 1.5 \u00c5',m.min_dist_ang>=1.5),('Minimum receptor distance \u2265 1.8 \u00c5 (used)',m.min_dist_ang>=1.8),('Minimum receptor distance \u2265 2.0 \u00c5',m.min_dist_ang>=2.0),('Minimum receptor distance \u2265 2.2 \u00c5',m.min_dist_ang>=2.2),('Minimum receptor distance \u2265 2.5 \u00c5',m.min_dist_ang>=2.5),('Minimum receptor distance \u2265 3.0 \u00c5',m.min_dist_ang>=3.0),('No bond-length flag',m.chem_ok),('No bond-length flag and \u2265 1.8 \u00c5',m.chem_ok&(m.min_dist_ang>=1.8))):
        s=m[mask]; k=int(s.groupby('c').ok.max().sum()); lo,hi=wil(k,n_c); out[lab]=dict(ret=100*mask.mean(),has=int(s.c.nunique()),n=n_c,sp=100*k/n_c,lo=lo,hi=hi)
    sens[name]=out; sens[name+'_flag']=100*(~m.chem_ok).mean()
json.dump(dict(anch=anch,sens=sens),open(r"C:\Users\ADMIN\AppData\Local\Temp\rdk\th_numbers.json",'w'),indent=1)
print(json.dumps(anch,indent=0)); print(sens['PB_flag'],sens['CASF_flag'])
