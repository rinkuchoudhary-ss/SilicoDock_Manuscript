import pandas as pd, numpy as np
from math import sqrt
from scipy.stats import spearmanr
B=r"C:\Users\ADMIN\Desktop\SilicoScientia\AI_Docking_tool\Manuscript_Revision\benchmarking data"
P=pd.read_csv(r"C:\Users\ADMIN\AppData\Local\Temp\rdk\ligprops.csv")
def wil(k,n,z=1.96):
    p=k/n; d=1+z*z/n; ce=(p+z*z/(2*n))/d; h=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return 100*(ce-h),100*(ce+h)
out=[]
for s,name,rm in (('CASF','CASF',B+r"\casf_new_pose_rmsd.csv"),('PB','PoseBusters',B+r"\posebuster_pose_rmsd.csv")):
    r=pd.read_csv(rm); r['c']=r.entry_id.str.lower() if s=='CASF' else r.entry_id; r['ok']=r.rmsd_ang<=2
    g=r.groupby('c').agg(best=('rmsd_ang','min'),frac=('ok','mean'),med=('rmsd_ang','median')); g['sp']=(g.best<=2).astype(int)
    e=pd.read_csv(r"C:\Users\ADMIN\AppData\Local\Temp\rdk\e2e_"+('CASF' if s=='CASF' else 'PoseBusters')+".csv"); e['ok']=e.rmsd_ang<=2
    g['t1']=e[e['rank']==1].set_index('complex').ok.reindex(g.index).fillna(False).astype(int)
    d=P[P.set==s].set_index('c').join(g); d['set']=s; out.append(d.reset_index())
D=pd.concat(out); D.to_csv(r"C:\Users\ADMIN\AppData\Local\Temp\rdk\flex_all.csv",index=False)
def table(d,col,edges,labels):
    d=d.assign(b=pd.cut(d[col],edges,labels=labels,right=True)); rows=[]
    for b,x in d.groupby('b',observed=True):
        n=len(x);k=int(x.sp.sum());lo,hi=wil(k,n);k1=int(x.t1.sum());lo1,hi1=wil(k1,n)
        rows.append(f'   {str(b):10s} n={n:3d}  sampling {100*k/n:5.1f}% [{lo:4.1f}-{hi:4.1f}]  per-sample {100*x.frac.mean():5.1f}%  top1(SMiDeN) {100*k1/n:5.1f}% [{lo1:4.1f}-{hi1:4.1f}]')
    return '\n'.join(rows)
for s,d in (('CASF',D[D.set=='CASF']),('PB',D[D.set=='PB'])):
    print('\n=====',s,len(d))
    print(' rotatable bonds'); print(table(d,'rotb',[-1,2,5,8,100],['0-2','3-5','6-8','>=9']))
    print(' heavy atoms'); print(table(d,'ha',[0,15,25,35,100],['<=15','16-25','26-35','>35']))
    print(' rings'); print(table(d,'rings',[-1,2,3,4,100],['0-2','3','4','>=5']))
    print(' ring systems'); print(table(d,'ringsys',[-1,1,2,3,100],['0-1','2','3','>=4']))
    print(' macrocycle (ring>=12)'); 
    for v,x in d.groupby('macro'):
        n=len(x);k=int(x.sp.sum());lo,hi=wil(k,n);print(f'   macro={v} n={n} sampling {k}/{n}={100*k/n:.1f}% [{lo:.1f}-{hi:.1f}] top1 {int(x.t1.sum())}/{n}')
    for col in ('rotb','ha','rings','ringsys','maxring'):
        r1,p1=spearmanr(d[col],d.best); r2,p2=spearmanr(d[col],d.frac); print(f'   Spearman {col:8s} vs best RMSD rho={r1:+.2f} p={p1:.2g} | vs fraction of samples correct rho={r2:+.2f} p={p2:.2g}')
