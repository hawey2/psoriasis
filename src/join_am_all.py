import gzip, json, collections
rows=json.load(open('out/clinvar_genes.json'))
want=collections.defaultdict(list)
for r in rows:
    if 'missense_variant' in r['mc']:
        want[(r['chrom'],r['pos'],r['ref'],r['alt'])].append(r)
allsc=collections.defaultdict(list)
with gzip.open('data/AlphaMissense_hg38.tsv.gz','rt') as fh:
    for line in fh:
        if line[0]=='#': continue
        p=line.rstrip('\n').split('\t')
        if p[0]=='CHROM': continue
        k=(p[0][3:],int(p[1]),p[2],p[3])
        if k in want:
            allsc['|'.join(map(str,k))].append(dict(uniprot=p[5],tx=p[6],pv=p[7],score=float(p[8]),cls=p[9]))
json.dump(allsc,open('out/am_all_transcripts.json','w'))
print('keys with >=1 AM row:',len(allsc),' with >1:',sum(1 for v in allsc.values() if len(v)>1))
