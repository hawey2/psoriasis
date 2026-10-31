import json, re, numpy as np, collections
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
rng=np.random.default_rng(1)
GENES=json.load(open("src/mech.json"))
PATH={"Pathogenic","Likely_pathogenic","Pathogenic/Likely_pathogenic"}
BEN={"Benign","Likely_benign","Benign/Likely_benign"}
rows=[r for r in json.load(open('out/clinvar_am.json')) if 'am_score' in r]
for r in rows:
    r['y']=1 if r['clnsig'] in PATH else 0 if r['clnsig'] in BEN else None
    m=re.match(r'^([A-Z])(\d+)([A-Z])$', r.get('prot_var',''))
    r['respos']=int(m.group(2)) if m else None
lab=[r for r in rows if r['y'] is not None and r['respos'] is not None]
print(f"labelled missense with residue position: {len(lab)}")

def thr_at_fpr(y,s,target=0.05):
    neg=np.sort(np.array([si for si,yi in zip(s,y) if yi==0]))[::-1]
    k=max(1,int(np.floor(target*len(neg))))
    return neg[k-1]

def sens(y,s,t): 
    p=[si for si,yi in zip(s,y) if yi==1]; return float(np.mean([si>=t for si in p]))
def fpr(y,s,t):
    n=[si for si,yi in zip(s,y) if yi==0]; return float(np.mean([si>=t for si in n]))

print("\n=== EXPERIMENT 2: mechanism-conditional calibration of AlphaMissense (5-fold CV, FPR target 5%) ===")
res2={}
for m in ['LOF','GOF','MIXED']:
    sub=[r for r in lab if GENES[r['gene']]==m]
    y=np.array([r['y'] for r in sub]); s=np.array([r['am_score'] for r in sub])
    gy=np.array([r['y'] for r in lab]); gs=np.array([r['am_score'] for r in lab])
    skf=StratifiedKFold(5,shuffle=True,random_state=0)
    sg,sm,fg,fm=[],[],[],[]
    for tr,te in skf.split(s,y):
        # global threshold fitted on all mechanisms' training-like pool (held-out test genes' rows excluded)
        mask=np.ones(len(lab),bool)
        te_ids={id(sub[i]) for i in te}
        glob=[(r['y'],r['am_score']) for r in lab if id(r) not in te_ids]
        tg=thr_at_fpr([a for a,_ in glob],[b for _,b in glob])
        tm=thr_at_fpr(y[tr],s[tr])
        sg.append(sens(y[te],s[te],tg)); fg.append(fpr(y[te],s[te],tg))
        sm.append(sens(y[te],s[te],tm)); fm.append(fpr(y[te],s[te],tm))
    print(f"{m:6s} n={len(sub):4d}  global-threshold sens={np.mean(sg):.3f} (FPR {np.mean(fg):.3f})"
          f"   mechanism-threshold sens={np.mean(sm):.3f} (FPR {np.mean(fm):.3f})")
    res2[m]=dict(n=len(sub),sens_global=float(np.mean(sg)),fpr_global=float(np.mean(fg)),
                 sens_mech=float(np.mean(sm)),fpr_mech=float(np.mean(fm)))

print("\n=== EXPERIMENT 3: does 1D positional clustering add signal on GOF genes? (5-fold CV) ===")
gof=[r for r in lab if GENES[r['gene']]=='GOF']
y=np.array([r['y'] for r in gof]); 
print(f"GOF-gene labelled missense n={len(gof)} (P={int(y.sum())}, B={len(y)-int(y.sum())}) over genes {sorted({r['gene'] for r in gof})}")
skf=StratifiedKFold(5,shuffle=True,random_state=0)
auc_am,auc_both,auc_dist=[],[],[]
for tr,te in skf.split(np.zeros(len(gof)),y):
    trs=[gof[i] for i in tr]; tes=[gof[i] for i in te]
    byg=collections.defaultdict(list)
    for r in trs:
        if r['y']==1: byg[r['gene']].append(r['respos'])
    def dist(r):
        ps=byg.get(r['gene'],[])
        return min((abs(r['respos']-p) for p in ps), default=1000)
    Xtr=np.array([[r['am_score'], np.log1p(dist(r))] for r in trs])
    Xte=np.array([[r['am_score'], np.log1p(dist(r))] for r in tes])
    ytr=np.array([r['y'] for r in trs]); yte=np.array([r['y'] for r in tes])
    if len(set(yte))<2: continue
    lr=LogisticRegression(max_iter=1000).fit(Xtr,ytr)
    auc_both.append(roc_auc_score(yte,lr.predict_proba(Xte)[:,1]))
    auc_am.append(roc_auc_score(yte,Xte[:,0]))
    auc_dist.append(roc_auc_score(yte,-Xte[:,1]))
print(f"AlphaMissense alone      AUROC={np.mean(auc_am):.3f}")
print(f"positional clustering    AUROC={np.mean(auc_dist):.3f}")
print(f"AM + positional (logreg) AUROC={np.mean(auc_both):.3f}")
json.dump(dict(exp2=res2,exp3=dict(n=len(gof),P=int(y.sum()),
    auroc_am=float(np.mean(auc_am)),auroc_dist=float(np.mean(auc_dist)),
    auroc_both=float(np.mean(auc_both)))),open('out/metrics_exp23.json','w'),indent=1)
