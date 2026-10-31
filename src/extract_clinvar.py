import gzip, re, json, sys, collections

GENES = {
 # psoriasis-spectrum core, with the literature-reported dominant mechanism
 "IL36RN":"LOF","MPO":"LOF","AP1S3":"LOF","SERPINA3":"LOF","CARD14":"GOF",
 # broader psoriasis / IL-23-IL-17 axis and autoinflammation comparators
 "IL23R":"MIXED","ERAP1":"MIXED","TYK2":"MIXED","TNFAIP3":"LOF","NFKBIA":"MIXED",
 "TRAF3IP2":"LOF","STAT3":"MIXED","IL12B":"MIXED","NOD2":"GOF","NLRP1":"GOF",
 "NLRP3":"GOF","IL1RN":"LOF","CARD11":"MIXED","RELA":"LOF","OTULIN":"LOF",
 "IL17RA":"LOF","IL36G":"MIXED","IL36A":"MIXED","IL1RL2":"MIXED","CDSN":"LOF",
 "LCE3D":"MIXED","GJB2":"MIXED","FUT2":"LOF","POMP":"GOF","PSMB8":"LOF",
}
PATH={"Pathogenic","Likely_pathogenic","Pathogenic/Likely_pathogenic"}
BEN={"Benign","Likely_benign","Benign/Likely_benign"}

def info_dict(s):
    d={}
    for f in s.split(';'):
        if '=' in f:
            k,v=f.split('=',1); d[k]=v
    return d

rows=[]
with gzip.open('/home/claude/work/data/clinvar_20250504.vcf.gz','rt') as fh:
    for line in fh:
        if line[0]=='#': continue
        p=line.rstrip('\n').split('\t')
        chrom,pos,vid,ref,alt=p[0],p[1],p[2],p[3],p[4]
        inf=info_dict(p[7])
        gi=inf.get('GENEINFO','')
        genes=[g.split(':')[0] for g in gi.split('|') if g]
        hit=[g for g in genes if g in GENES]
        if not hit: continue
        rows.append(dict(chrom=chrom,pos=int(pos),ref=ref,alt=alt,vid=vid,
            gene=hit[0], genes=genes,
            clnsig=inf.get('CLNSIG',''), rev=inf.get('CLNREVSTAT',''),
            mc=inf.get('MC',''), clndn=inf.get('CLNDN',''),
            hgvs=inf.get('CLNHGVS',''), af_exac=inf.get('AF_EXAC',''),
            af_gnomad=inf.get('AF_ESP','')))
print('rows',len(rows))
json.dump(rows,open('/home/claude/work/out/clinvar_genes.json','w'))
c=collections.Counter((r['gene'],'P' if r['clnsig'] in PATH else 'B' if r['clnsig'] in BEN else 'U') for r in rows)
for g in GENES:
    print(g, GENES[g], 'P',c[(g,'P')],'B',c[(g,'B')],'U',c[(g,'U')])
