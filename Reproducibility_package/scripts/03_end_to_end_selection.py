import pandas as pd, numpy as np, re
from math import sqrt
B=r"C:\Users\ADMIN\Desktop\SilicoScientia\AI_Docking_tool\Manuscript_Revision\benchmarking data"; D=r"C:\Users\ADMIN\Desktop\SilicoScientia\AI_Docking_tool\Manuscript_Revision\benchmarking_smiden_data"
def wil(k,n,z=1.96):
    p=k/n; d=1+z*z/n; ce=(p+z*z/(2*n))/d; h=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return 100*(ce-h),100*(ce+h)
num=lambda s:int(re.search(r'(\d+)',s).group(1))
res={}
for name,rm,pp,rep,n_all in (('PoseBusters',B+r"\posebuster_pose_rmsd.csv",D+r"\per_pose_ranking_posebuster.csv",D+r"\ranking_report_posebuster.csv",428),('CASF',B+r"\casf_new_pose_rmsd.csv",D+r"\per_pose_ranking.csv",D+r"\ranking_report.csv",285)):
    r=pd.read_csv(rm); r['complex']=r.entry_id.str.lower() if name=='CASF' else r.entry_id; r['pn']=r.pose_file.map(num); r=r[['complex','pn','rmsd_ang']]
    s=pd.read_csv(pp); s['pn']=s.pose.map(num)
    if name=='CASF': s['complex']=s.complex.str.lower()
    m=s.merge(r,on=['complex','pn'],how='left'); print('\n==',name,'rows',len(s),'unmatched',m.rmsd_ang.isna().sum(),'complexes scored',s.complex.nunique(),'of',n_all)
    if name=='PoseBusters':
        rp=pd.read_csv(rep); t=m[m['rank']==1].merge(rp[['complex','top1_rmsd','best_rmsd_available']],on='complex')
        d=(t.rmsd_ang-t.top1_rmsd).abs(); print('mapping check vs report top1_rmsd: max abs diff',round(d.max(),3),'n equal(<0.01)',int((d<0.01).sum()),'of',len(t))
    top=m[m['rank']==1]; n=len(top); k=int((top.rmsd_ang<=2).sum()); lo,hi=wil(k,n)
    print(f'top-1 <=2A scored: {k}/{n} = {100*k/n:.1f}% [{lo:.1f}-{hi:.1f}]'); lo,hi=wil(k,n_all); print(f'  unscored as failures: {k}/{n_all} = {100*k/n_all:.1f}% [{lo:.1f}-{hi:.1f}]')
    bp=m.groupby('complex').rmsd_ang.min(); print(f'  oracle (scored poses): {int((bp<=2).sum())}/{len(bp)} = {100*(bp<=2).mean():.1f}%')
    fr=m.assign(ok=m.rmsd_ang<=2).groupby('complex').ok.mean(); print(f'  random-pick expectation: {100*fr.mean():.1f}%')
    for t_ in (2,3): kk=int((m[m['rank']<=t_].groupby('complex').rmsd_ang.min()<=2).sum()); print(f'  top-{t_}: {kk}/{n} = {100*kk/n:.1f}%')
    # on the full-denominator sampling success for the same unscored complexes
    allr=pd.read_csv(rm); allr['complex']=allr.entry_id.str.lower() if name=='CASF' else allr.entry_id
    miss=sorted(set(allr.complex)-set(s.complex)); mo=allr[allr.complex.isin(miss)].groupby('complex').rmsd_ang.min(); print('  unscored complexes',len(miss),'of which oracle<=2A:',int((mo<=2).sum()))
    m.to_csv(r"C:\Users\ADMIN\AppData\Local\Temp\rdk\e2e_"+name+".csv",index=False)
