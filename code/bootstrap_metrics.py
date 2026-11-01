# Estimate conditional school-bootstrap intervals and paired transfer contrasts.
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ap=argparse.ArgumentParser()
ap.add_argument('results_dir',type=Path)
ap.add_argument('--replicates',type=int,default=2000)
args=ap.parse_args(); out=args.results_dir
w=pd.read_csv(out/'within_country_predictions_all_plausible_values.csv')
e=pd.read_csv(out/'external_validation_predictions_all_plausible_values.csv')
keys=['target_code','feature_set','plausible_value','model','row_uid','school_uid']
d=w.merge(e,on=keys,suffixes=('_within','_external'),validate='one_to_one')
assert len(d)==len(w)==len(e)
assert np.allclose(d.observed_within,d.observed_external)
rows=[]; deltas=[]
for gi,((code,fs,model),g) in enumerate(d.groupby(['target_code','feature_set','model'])):
    ug=sorted(g.school_uid.unique()); ns=len(ug)
    rng=np.random.default_rng(20261022+gi)
                                                                              
    counts=rng.multinomial(ns,np.ones(ns)/ns,size=args.replicates)
    boot={c:{m:[] for m in ['RMSE','MAE','R2','mean_error']} for c in ['within','external']}
    point={c:{m:[] for m in ['RMSE','MAE','R2','mean_error']} for c in ['within','external']}
    for pv,pp in g.groupby('plausible_value'):
        y=pp.observed_within.to_numpy()
        stats=pd.DataFrame({'school':pp.school_uid.to_numpy(),'n':1,'ys':y,'yss':y*y})
        for c in boot:
            err=pp['predicted_'+c].to_numpy()-y
            stats[c+'_sse']=err*err;stats[c+'_sae']=np.abs(err);stats[c+'_err']=err
        a=stats.groupby('school').sum().reindex(ug)
        n=counts@a.n.to_numpy();ys=counts@a.ys.to_numpy();yss=counts@a.yss.to_numpy()
        for c in boot:
            sse=counts@a[c+'_sse'].to_numpy()
            vals={'RMSE':np.sqrt(sse/n),'MAE':counts@a[c+'_sae'].to_numpy()/n,
                  'R2':1-sse/(yss-ys*ys/n),'mean_error':counts@a[c+'_err'].to_numpy()/n}
            nt=a.n.sum();yt=a.ys.sum();sst=a.yss.sum()-yt*yt/nt;se=a[c+'_sse'].sum()
            pts={'RMSE':np.sqrt(se/nt),'MAE':a[c+'_sae'].sum()/nt,
                 'R2':1-se/sst,'mean_error':a[c+'_err'].sum()/nt}
            for m in vals:boot[c][m].append(vals[m]);point[c][m].append(pts[m])
    assert len(g.plausible_value.unique())==5
    for c in boot:
        for m in boot[c]:
            b=np.mean(boot[c][m],axis=0);pt=np.mean(point[c][m]);lo,hi=np.quantile(b,[.025,.975])
            rows.append(dict(target_code=code,feature_set=fs,model=model,context=c,
                metric=m,estimate=pt,lower_95=lo,upper_95=hi,schools=ns,
                students=len(g)//5,replicates=args.replicates,
                method='School percentile bootstrap; conditional on fitted models and five PV draws'))
    for m in boot['within']:
        b=np.mean(boot['external'][m],axis=0)-np.mean(boot['within'][m],axis=0)
        lo,hi=np.quantile(b,[.025,.975])
        deltas.append(dict(target_code=code,feature_set=fs,model=model,metric=m,
            external_minus_within=np.mean(point['external'][m])-np.mean(point['within'][m]),
            lower_95=lo,upper_95=hi,schools=ns,students=len(g)//5,replicates=args.replicates))
pd.DataFrame(rows).to_csv(out/'conditional_school_bootstrap_metrics.csv',index=False)
pd.DataFrame(deltas).to_csv(out/'paired_external_minus_within_bootstrap.csv',index=False)
print(f'Wrote {len(rows)} intervals and {len(deltas)} paired contrasts from {len(d):,} matched predictions.')
