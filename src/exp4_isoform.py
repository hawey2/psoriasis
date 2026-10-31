import gzip, json, collections, numpy as np
GENES=json.load(open("src/mech.json"))
PATH={"Pathogenic","Likely_pathogenic","Pathogenic/Likely_pathogenic"}
BEN={"Benign","Likely_benign","Benign/Likely_benign"}
rows=json.load(open('out/clinvar_genes.json'))
want=collections.defaultdict(list)
for r in rows:
    if 'missense_variant' in r['mc']:
        want[(r['chrom'],r['pos'],r['ref'],r['alt'])].append(r)
iso=collections.defaultdict(list)
with gzip.open('data/AlphaMissense_isoforms_hg38.tsv.gz','rt') as fh:
    for line in fh:
        if line[0]=='#': continue
        p=line.rstrip('\n').split('\t')
        if p[0]=='CHROM': continue
        k=(p[0][3:],int(p[1]),p[2],p[3])
        if k in want:
            iso[k].append((p[5],p[6],float(p[7]),p[8]))  # transcript, prot_var, score, class
print('variants with isoform-level scores:',len(iso),'| multi-isoform:',sum(1 for v in iso.values() if len(v)>1))
recs=[]
for k,v in iso.items():
    if len(v)<2: continue
    s=[x[2] for x in v]; c={x[3] for x in v}
    r0=want[k][0]
    y=1 if r0['clnsig'] in PATH else 0 if r0['clnsig'] in BEN else None
    recs.append(dict(gene=r0['gene'],mech=GENES[r0['gene']],n_iso=len(v),rng=max(s)-min(s),
                     flip=len(c)>1, y=y, clnsig=r0['clnsig'],
                     pv=v[0][1], smin=min(s), smax=max(s)))
print(f"multi-isoform variants analysed: {len(recs)}")
def summ(name,sub):
    if not sub: return
    rg=np.array([r['rng'] for r in sub]); fl=np.mean([r['flip'] for r in sub])
    print(f"{name:34s} n={len(sub):5d}  median score range={np.median(rg):.3f}  mean={rg.mean():.3f}  class-call flips across isoforms={fl*100:.1f}%")
summ('all multi-isoform',recs)
for m in ['LOF','GOF','MIXED']: summ(f'mechanism={m}',[r for r in recs if r['mech']==m])
core={'IL36RN','MPO','AP1S3','SERPINA3','CARD14'}
summ('psoriasis core 5 genes',[r for r in recs if r['gene'] in core])
for g in sorted(core): summ(f'  {g}',[r for r in recs if r['gene']==g])
summ('ClinVar pathogenic only',[r for r in recs if r['y']==1])
summ('ClinVar benign only',[r for r in recs if r['y']==0])
summ('ClinVar VUS/conflicting only',[r for r in recs if r['y'] is None])
print('\nlargest isoform disagreements in psoriasis core genes:')
for r in sorted([r for r in recs if r['gene'] in core],key=lambda x:-x['rng'])[:12]:
    print(f"  {r['gene']:9s} {r['pv']:8s} {r['n_iso']:2d} isoforms  score {r['smin']:.3f}-{r['smax']:.3f}  ClinVar={r['clnsig']}")
json.dump(dict(n_multi=len(recs),
  by_mech={m:dict(n=len([r for r in recs if r['mech']==m]),
                  median_range=float(np.median([r['rng'] for r in recs if r['mech']==m])),
                  flip_pct=float(np.mean([r['flip'] for r in recs if r['mech']==m])*100)) for m in ['LOF','GOF','MIXED']},
  core=dict(n=len([r for r in recs if r['gene'] in core]),
            median_range=float(np.median([r['rng'] for r in recs if r['gene'] in core])),
            flip_pct=float(np.mean([r['flip'] for r in recs if r['gene'] in core])*100))),
  open('out/metrics_exp4.json','w'),indent=1)
