#!/usr/bin/env python3
"""箱 2 つと等面四面体に折れる共通展開図のデータを作る（分析/conway.py の結果 + 箱の折り方 → data/data_tetra_*.json）。
usage: python3 make_tetra.py   （研究リポジトリの中で実行。箱の折りは fold_export.export、四面体は格子の折り線と二面角から作る）"""
import sys, os, json, math, itertools
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0,os.path.join(ROOT,'分析')); sys.path.insert(0,os.path.join(ROOT,'開いた箱','src'))
import conway as CW, features as FT
os.chdir(os.path.join(ROOT,'開いた箱'))
import fold_export as fe

def dihedrals(a,b,c):
    """辺² が a,b,c（各 2 本ずつ）の等面四面体の，辺² → 二面角（内側）。直方体 p×q×r の対角の 4 頂点で作る。"""
    p2,q2,r2=(a+c-b)/2,(a+b-c)/2,(b+c-a)/2; assert min(p2,q2,r2)>0
    p,q,r=map(math.sqrt,(p2,q2,r2)); V=np.array([(0,0,0),(p,q,0),(p,0,r),(0,q,r)],float)
    out={}
    for i,j in itertools.combinations(range(4),2):
        e=V[j]-V[i]; e/=np.linalg.norm(e); z,w=[V[k]-V[i] for k in range(4) if k not in (i,j)]
        z-=z.dot(e)*e; w-=w.dot(e)*e
        out[round(float((V[j]-V[i])@(V[j]-V[i])))]=math.acos(z@w/np.linalg.norm(z)/np.linalg.norm(w))
    assert sorted(out)==sorted({a,b,c}),(a,b,c)   # 辺²が等しい対は二面角も等しい（二等辺の面でもよい）
    return out

def tetra_target(cells,edges,cuts=None,tag=''):
    """cells の周から Conway の分け方を探し，edges（辺² 3 つ）の四面体になるものの格子から折り線と折り角を作る。cuts を渡せばその分け方を使う。"""
    s=CW.boundary(cells)
    if cuts is None:
        for cuts in CW.factorizations(s):
            t=CW.check(s,cuts)
            if t and t[0]==tuple(edges) and not t[1]: break
        else: raise SystemExit(f'見つからない {edges}')
    P,mB,u,v=CW.lattice(s,cuts); sx=min(x for x,y in cells); sy=min(y for x,y in cells if x==sx)
    corners={(x+i,y+j) for x,y in cells for i in (0,1) for j in (0,1)}
    assert all((x+sx,y+sy) in corners for x,y in P)   # 周の開始点の取り違えがないこと
    mB=(mB[0]+2*sx,mB[1]+2*sy); w=(v[0]-u[0],v[1]-u[1]); D=abs(u[0]*v[1]-u[1]*v[0])
    dih=dihedrals(*edges); dirs=[u,v,w]
    xs=[x for x,y in cells]; ys=[y for x,y in cells]; box=(min(xs),min(ys),max(xs)+1,max(ys)+1)
    lines=[]   # 格子点を通る 3 方向の直線（1 倍座標の 2 点）。紙の外接長方形を横切るものだけ
    R=int(2*(box[2]-box[0]+box[3]-box[1]))+2
    for d in dirs:
        for k in range(-R*10,R*10):
            # 方向 d の直線は mB + (i u + j v) を通る。d=u なら j=k，d=v なら i=k，d=w なら i=k（w の格子点は (k,0)+t w）
            base=(mB[0]+(k*v[0] if d==u else k*u[0]),mB[1]+(k*v[1] if d==u else k*u[1]))
            a=(base[0]/2,base[1]/2); b=(a[0]+d[0]/2,a[1]+d[1]/2)
            c=[(d[0]*(y-a[1])-d[1]*(x-a[0])) for x in (box[0],box[2]) for y in (box[1],box[3])]
            if min(c)<=0<=max(c) and not (min(c)==0==max(c)) or min(c)<0<max(c): lines.append((a,b))
    def angle(h0,h1,poly):
        m=((h0[0]+h1[0])-mB[0],(h0[1]+h1[1])-mB[1])      # ヒンジの中点（2 倍座標）- mB
        cx=sum(p[0] for p in poly)/len(poly); cy=sum(p[1] for p in poly)/len(poly)
        side=(h1[0]-h0[0])*(cy-h0[1])-(h1[1]-h0[1])*(cx-h0[0])   # 子が軸の左なら正
        for d in dirs:
            if abs(d[0]*(h1[1]-h0[1])-d[1]*(h1[0]-h0[0]))>1e-6: continue       # 向きが違う
            off=(d[0]*m[1]-d[1]*m[0])/D
            if abs(off-round(off))<1e-6: return (-1 if side>0 else 1)*(math.pi-dih[d[0]**2+d[1]**2])   # 紙の表を外側にして折る
        return 0.0
    def colors(P3):
        key={}; cols=[]
        for poly in P3:
            p0,p1,p2=map(np.array,poly[:3]); n=np.cross(p1-p0,p2-p0); n/=np.linalg.norm(n)
            cols.append(key.setdefault((*np.round(n,4),round(float(n@p0),4)),len(key)))
        assert len(key)==4,key      # 面が 4 つの平面に収まる
        return [fe.COLORS[c] for c in cols]
    g=math.gcd(math.gcd(edges[0],edges[1]),edges[2])
    return {'label':'等面四面体（辺 '+':'.join(f'√{e//g}' for e in edges)+'）'+tag,'lines':lines,'angle':angle,'colors':colors}

def build(conway,spec_folds,tet_pick,labels,specs,out,title):
    C=json.load(open(conway)); ent=[c for c in C if tet_pick(c)][0]; cells=[tuple(c) for c in ent['cells']]
    rows=FT.load_rows(spec_folds[0]); r=[r for r in rows if set(r['cells'])==set(cells)][0]
    row={'cells':cells,'folds':spec_folds[1](r)}
    tets=[tetra_target(cells,e) for (e,fl) in ent['tetra'] if not fl]
    fe.export(row,specs,out,title,labels=labels,extra=tets)

def tetra_folds(cells):
    """四面体ごとに，折り方の違う（回転の中心の格子が違う，紙の対称で移らない）分け方を全部 [(辺², cuts, 番号)]"""
    import tetra_all as TA
    s=CW.boundary(cells); sx=min(x for x,y in cells); base=(sx,min(y for x,y in cells if x==sx)); sym=TA.symmetries(cells)
    reps={}
    for cuts in CW.factorizations(s):
        t=CW.check(s,cuts)
        if not t or t[1]: continue
        L=TA.latt(s,cuts,base); R=reps.setdefault(t[0],[])
        imgs=[((M(*L[0])[0]+2*d[0],M(*L[0])[1]+2*d[1]),M(*L[1]),M(*L[2])) for M,d in sym]
        if not any(TA.same(r[0],I) for r in R for I in imgs): R.append((L,cuts))
    return [(e,c,j+1,len(R)) for e,R in sorted(reps.items()) for j,(L,c) in enumerate(R)]

def build_cells(cells,boxes,out,title,allfolds=False,only=None):
    """セルの一覧から直接作る（tetra_all.py の出力用）。箱の折り方は obox.kfull で全部出し，折り方ごとに 1 つの目標にする。"""
    import obox
    cells=[tuple(c) for c in cells]; folds=[]; specs=[]; labels=[]
    for b in boxes:
        K=obox.kfull(obox.surface(b),cells)
        for j,st in enumerate(K):
            folds.append([st]); specs.append(b); labels.append(b.replace('x','×')+(f'（折り方 {j+1}）' if len(K)>1 else ''))
    s=CW.boundary(cells); T=sorted({t for t in (CW.check(s,c) for c in CW.factorizations(s)) if t and not t[1]})
    tets=[tetra_target(cells,e) for e,fl in T] if not allfolds else \
         [tetra_target(cells,e,c,f'（折り方 {j}）' if n>1 else '') for e,c,j,n in tetra_folds(cells) if only is None or e==only]
    fe.export({'cells':cells,'folds':folds},specs,out,title,labels=labels,extra=tets)
    return len(folds),len(tets)

NEW=[  # tetra_all.py で見つけたもの（分析 §21）。(名前, セル, 箱, 題)
    ('tetra_a14',[(0,0),(1,0),(2,0),(3,0),(4,0),(5,0),(6,0),(4,1),(5,1),(6,1),(7,1),(8,1),(9,1),(10,1)],['1x1x3'],'面積 14：1×1×3 と等面四面体 11 種',False),
    ('tetra_a22_115',[(0,y) for y in range(11)]+[(1,y) for y in range(6,17)],['1x1x5'],'面積 22：1×1×5 と等面四面体 16 種',False),
    ('tetra_a22_115_13',[(0,y) for y in range(11)]+[(1,y) for y in range(6,17)],['1x1x5'],'面積 22：1×1×5 と，同じ等面四面体に 13 通り','MAX'),
    ('tetra_stair26',[(k,k) for k in range(13)]+[(k,k+1) for k in range(13)],[],'26 マスの階段：等面四面体 10 種',False),
    # tetra_stair26_16（同じ四面体に 16 通り）は test_tetra で閉じることを確認したが，データが 17 MB になるので載せない
    ('tetra_cube_w',[(0,0),(1,0),(1,1),(2,1),(2,2),(3,2)],['1x1x1'],'立方体の展開図（階段形）：立方体と，等面四面体 4 種に 9 通り',True),
]

if __name__=='__main__' and sys.argv[1:2]==['new']:
    A=os.path.join(ROOT,'アプリ')
    for name,cells,boxes,title,allf in NEW:
        if len(sys.argv)>2 and name not in sys.argv[2:]: continue
        only=None
        if allf=='MAX':   # 折り方が一番多い四面体だけ，その折り方を全部
            from collections import Counter
            only=Counter(e for e,c,j,n in tetra_folds(cells)).most_common(1)[0][0]
        print(name,build_cells(cells,boxes,f'{A}/data/data_{name}.json',title,allfolds=bool(allf),only=only))
elif __name__=='__main__':
    A=os.path.join(ROOT,'アプリ'); W=os.path.join(ROOT,'分析','work')
    build(f'{W}/conway_a22.json',(f'{W}/rows_a22_common.json',lambda r:[r['folds'][0],[r['folds'][1][0]],[r['folds'][1][1]]]),
          lambda c:c['k']==[1,2],['1×1×5','1×2×3（折り方 1）','1×2×3（折り方 2）'],['1x1x5','1x2x3','1x2x3'],
          f'{A}/data/data_tetra_a22.json','面積 22 の共通展開図：1×1×5，1×2×3，等面四面体')
    build(f'{W}/conway_a34.json',(os.path.join(ROOT,'開いた箱','work','rows_m34.json.gz'),lambda r:[r['folds'][2],r['folds'][3]]),
          lambda c:len(c['tetra'])==6,['1×1×8','1×2×5'],['1x1x8','1x2x5'],
          f'{A}/data/data_tetra_a34.json','面積 34 の共通展開図：1×1×8，1×2×5，等面四面体 6 種')
