import sys,math; sys.path.insert(0,sys.argv[1])
from xl import read
D=sys.argv[2]
def interp(t, ts, cs):
    for i in range(1,len(ts)):
        if ts[i]>=t:
            f=(t-ts[i-1])/(ts[i]-ts[i-1]); return cs[i-1]+f*(cs[i]-cs[i-1])
def trap(ts,cs): return sum((ts[i]-ts[i-1])*(cs[i]+cs[i-1])/2 for i in range(1,len(ts)))
gm=lambda x: math.exp(sum(map(math.log,x))/len(x))
obs={}
for f in ['ER_ref_all11_96h.xlsx','ER_test6_96h.xlsx']:
    for n,rows in list(read(D+'/'+f).items())[1:]:
        if 'PK' in n: continue
        obs[n.split(' - ')[0]]=('train' if 'all11' in f else 'test',[(r[0],r[1]) for r in rows[1:] if isinstance(r[0],float) and isinstance(r[1],float) and r[0]>0])
def curve(f,col=1):
    d=list(read(D+'/'+f).items()); sim=d[0][1][1:]
    return [r[0] for r in sim],[r[col] for r in sim], [n for n,_ in d[1:]], d
bins={'0-2':(0,2),'2-6':(2,6),'6-12':(6,12.1),'12-24':(12.1,24.1)}
for f in sys.argv[3:]:
    col=2 if f.startswith('Confidence') else 1
    st,sc,names,d=curve(f,col)
    pc=max(sc); 
    print(f'== {f:26s} end={st[-1]} Cmax={pc*1e3:.3f} nM @ {st[sc.index(pc)]} h  AUC0-96={trap(st,sc)*60:.3f} µM·min  C24={interp(24,st,sc)*1e3:.3f}  sheets={len(names)}')
    for grp in ('train','test','all'):
        cr=[];ar=[];acc={k:[] for k in bins}
        for s,(g,ob) in obs.items():
            if grp!='all' and g!=grp: continue
            ot=[o[0] for o in ob]; oc=[o[1] for o in ob]
            pC=max(c for t,c in zip(st,sc) if t<=ot[-1]); pg=[interp(t,st,sc) for t in ot]
            cr.append(pC/max(oc)); ar.append(trap(ot,pg)/trap(ot,oc))
            for t,p,o in zip(ot,pg,oc):
                for k,(a,b) in bins.items():
                    if a<t<=b: acc[k].append(math.log(p/o))
        print(f'   {grp:5s} Cmax {gm(cr):.3f}  AUC {gm(ar):.3f}  |', ' '.join(f'{k}h:{math.exp(sum(v)/len(v)):.2f}' for k,v in acc.items()))
