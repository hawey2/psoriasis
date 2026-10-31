import json, collections, numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
rng=np.random.default_rng(0)
GENES=json.load(open("src/mech.json"))
PATH={"Pathogenic","Likely_pathogenic","Pathogenic/Likely_pathogenic"}
BEN={"Benign","Likely_benign","Benign/Likely_benign"}
rows=[r for r in json.load(open('out/clinvar_am.json')) if 'am_score' in r]
for r in rows:
    r['y']=1 if r['clnsig'] in PATH else 0 if r['clnsig'] in BEN else None
lab=[r for r in rows if r['y'] is not None]
print(f"missense with AlphaMissense score: {len(rows)}; labelled P/B: {len(lab)}")

def boot(y,s,fn,n=2000):
    y=np.array(y);s=np.array(s);v=[]
    for _ in range(n):
        i=rng.integers(0,len(y),len(y))
        if len(set(y[i]))<2: continue
        v.append(fn(y[i],s[i]))
    return (np.percentile(v,2.5),np.percentile(v,97.5)) if v else (float('nan'),)*2

def report(name,sub):
    y=[r['y'] for r in sub]; s=[r['am_score'] for r in sub]
    npos,nneg=sum(y),len(y)-sum(y)
    if npos<3 or nneg<3: 
        print(f"{name:24s} n={len(y):5d} P={npos:4d} B={nneg:4d}  (too few for AUC)"); return None
    au=roc_auc_score(y,s); ap=average_precision_score(y,s)
    lo,hi=boot(y,s,roc_auc_score); alo,ahi=boot(y,s,average_precision_score)
    # class-call based operating point
    tp=sum(1 for r in sub if r['y']==1 and r['am_class']=='likely_pathogenic')
    fn_=sum(1 for r in sub if r['y']==1 and r['am_class']!='likely_pathogenic')
    fp=sum(1 for r in sub if r['y']==0 and r['am_class']=='likely_pathogenic')
    tn=sum(1 for r in sub if r['y']==0 and r['am_class']!='likely_pathogenic')
    sens=tp/(tp+fn_); fpr=fp/(fp+tn)
    print(f"{name:24s} n={len(y):5d} P={npos:4d} B={nneg:4d}  AUROC={au:.3f} [{lo:.3f},{hi:.3f}]  AUPRC={ap:.3f} [{alo:.3f},{ahi:.3f}]  sens@AMclass={sens:.3f} FPR={fpr:.3f}")
    return dict(name=name,n=len(y),P=npos,B=nneg,auroc=au,auroc_ci=[lo,hi],auprc=ap,auprc_ci=[alo,ahi],sens=sens,fpr=fpr)

res={}
print('\n== ALL genes pooled ==')
res['all']=report('ALL',lab)
print('\n== by curated disease mechanism ==')
for m in ['LOF','GOF','MIXED']:
    res[m]=report(f'mechanism={m}',[r for r in lab if GENES[r['gene']]==m])
print('\n== psoriasis-spectrum core genes only ==')
core={'IL36RN','MPO','AP1S3','SERPINA3','CARD14'}
res['core']=report('core (5 genes)',[r for r in lab if r['gene'] in core])
res['core_LOF']=report('core LOF (4 genes)',[r for r in lab if r['gene'] in core and GENES[r['gene']]=='LOF'])
res['core_GOF']=report('core GOF (CARD14)',[r for r in lab if r['gene']=='CARD14'])
print('\n== per gene (>=3 P and >=3 B) ==')
pg={}
for g in sorted({r['gene'] for r in lab}):
    o=report(f'{g} ({GENES[g]})',[r for r in lab if r['gene']==g])
    if o: pg[g]=o
res['per_gene']=pg

print('\n== VUS burden in psoriasis-spectrum core genes ==')
for g in sorted(core):
    sub=[r for r in rows if r['gene']==g]
    vus=[r for r in sub if r['y'] is None]
    c=collections.Counter(r['am_class'] for r in vus)
    print(f"{g:10s} missense total {len(sub):4d}  VUS/conflicting {len(vus):4d} -> AM likely_pathogenic {c['likely_pathogenic']:3d}  ambiguous {c['ambiguous']:3d}  likely_benign {c['likely_benign']:3d}")

print('\n== named literature disease alleles ==')
TARGETS={'CARD14':['G117S','E138A','E142K','E142G','R820W','M119T','R179H','H171N'],
         'IL36RN':['S113L','R48W','P76L','L27P','R10X'],
         'AP1S3':['F4C','R33W'],
         'MPO':['R569W','Y173C','M251T'],
         'SERPINA3':[]}
for g,vs in TARGETS.items():
    sub={r['prot_var']:r for r in rows if r['gene']==g}
    for v in vs:
        r=sub.get(v)
        if r: print(f"{g:8s} p.{v:6s} ClinVar={r['clnsig']:30s} AM={r['am_score']:.3f} ({r['am_class']})")
        else: print(f"{g:8s} p.{v:6s} not in ClinVar-missense∩AlphaMissense")
json.dump(res,open('out/metrics.json','w'),indent=1,default=float)
