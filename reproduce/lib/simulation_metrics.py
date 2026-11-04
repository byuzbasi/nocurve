"""Report a closed native simulation from saved predictions; no RNG or fitting."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

from metric_functions import losses, decomposition, mean, mcse, stratified_mcse


def write(directory,name,rows):
    with (directory/name).open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows({k:'NA' if v is None else v for k,v in r.items()} for r in rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--expected-cases',type=int,required=True)
    p.add_argument('--repetitions',type=int,required=True)
    a=p.parse_args()
    with (a.directory/'predictions.csv').open() as f:points=list(csv.DictReader(f))
    methods=['NOC','RAW','EB0','SMASH4','PS','AP','SS','KERNEL']
    data=defaultdict(list)
    for r in points:data[r['dataset']].append(r)
    assert len(data)==a.expected_cases and len(points)==a.expected_cases*256
    per_fit=[];groups=defaultdict(list)
    for dataset,rows in data.items():
        rows.sort(key=lambda r:float(r['x']));curve=rows[0]['curve'];rep=int(rows[0]['replication'])
        assert len(rows)==256
        for m in methods:
            errors=[float(r[m])-float(r['truth']) for r in rows]
            L=losses(errors)
            per_fit.append(dict(dataset=dataset,curve=curve,replication=rep,method=m,**L))
            groups[curve,m].append((rep,errors,rows,L))
    curves=sorted({c for c,m in groups});assert curves==['C1','C2','C3','C4']
    scenario=[];bias=[];means=[];lookup={}
    for (c,m),values in groups.items():
        values.sort(key=lambda x:x[0]);assert [z[0] for z in values]==list(range(1,a.repetitions+1))
        row=dict(curve=c,method=m,R=a.repetitions)
        for k in ('MSE','MAE','MAXAE'):
            v=[z[3][k] for z in values];row[k]=mean(v);row['MCSE_'+k]=mcse(v)
            lookup[c,m,k]=v
        scenario.append(row);bias.append(dict(curve=c,method=m,R=a.repetitions,
            **decomposition([z[1] for z in values])))
        for i in range(256):
            truth=float(values[0][2][i]['truth'])
            assert all(float(v[2][i]['truth'])==truth for v in values)
            means.append(dict(curve=c,method=m,R=a.repetitions,x=float(values[0][2][i]['x']),
                truth=truth,prediction=mean([float(v[2][i][m]) for v in values])))
    balanced=[];paired=[]
    for m in methods:
        row=dict(method=m,R_per_curve=a.repetitions,fixed_curves=len(curves))
        for k in ('MSE','MAE','MAXAE'):
            v=[lookup[c,m,k] for c in curves];row[k]=mean([mean(z) for z in v])
            row['MCSE_'+k]=stratified_mcse(v)
            if m!='NOC':
                diffs=[]
                for c in curves:
                    no,ref=lookup[c,'NOC',k],lookup[c,m,k]
                    d=[u-v for u,v in zip(no,ref)];diffs.append(d)
                    paired.append(dict(curve=c,reference=m,metric=k,difference=mean(d),
                        paired_MCSE=mcse(d),NOC_gain_percent=100*(1-mean(no)/mean(ref)) if mean(ref)>0 else None))
                no=mean([mean(lookup[c,'NOC',k]) for c in curves]);ref=row[k]
                paired.append(dict(curve='ALL_FIXED',reference=m,metric=k,difference=mean([mean(d) for d in diffs]),
                    paired_MCSE=stratified_mcse(diffs),NOC_gain_percent=100*(1-no/ref) if ref>0 else None))
        balanced.append(row)
    for name,rows in [('per_fit_metrics.tsv',per_fit),('scenario_metrics.tsv',scenario),
        ('balanced_metrics.tsv',balanced),('bias_variance.tsv',bias),('mean_curves.tsv',means),
        ('paired_metrics.tsv',paired)]:write(a.directory,name,rows)
    summary=dict(passed=True,new_fits=0,new_draws=0,datasets=a.expected_cases,
        complete_curves=len(per_fit),R_per_curve=a.repetitions,fixed_curves=curves,
        numerical_identity_max=max(abs(r['identity_residual'] or 0) for r in bias),
        warning='Fixed development shapes; fresh noise repetitions; not a population of functions.')
    with (a.directory/'metric_checks.json').open('x') as f:json.dump(summary,f,indent=2);f.write('\n')
    print('METRICS VERIFIED:',len(per_fit),'curves; all original-grid losses; no new fits or draws')


if __name__=='__main__':main()
