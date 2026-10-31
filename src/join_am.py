import gzip, json
rows=json.load(open('out/clinvar_genes.json'))
want={}
for r in rows:
    if 'missense_variant' in r['mc']:
        want.setdefault((r['chrom'],r['pos'],r['ref'],r['alt']),[]).append(r)
print('missense clinvar variants wanted:',len(want))
hits=0
with gzip.open('data/AlphaMissense_hg38.tsv.gz','rt') as fh:
    for line in fh:
        if line[0]=='#': continue
        p=line.rstrip('\n').split('\t')
        if p[0]=='CHROM': continue
        k=(p[0].replace('chr',''),int(p[1]),p[2],p[3])
        if k in want:
            for r in want[k]:
                r['am_score']=float(p[8]); r['am_class']=p[9]; r['prot_var']=p[7]; r['uniprot']=p[5]; r['tx']=p[6]
            hits+=1
print('matched',hits)
json.dump([r for v in want.values() for r in v],open('out/clinvar_am.json','w'))
