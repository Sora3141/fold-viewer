#!/usr/bin/env python3
"""data/data_tetra_*.json を全ヒンジ目標角で折ったとき，本当に閉じた形になるかを確かめる（箱は fold_export の target3d と相似か，四面体は 4 面・4 頂点・面積一致・辺の比）。
usage: python3 test_tetra.py   （研究リポジトリの中で実行）"""
import sys, os, json, math, re, glob
import numpy as np
from scipy.spatial import ConvexHull
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,os.path.join(ROOT,'開いた箱','src')); os.chdir(os.path.join(ROOT,'開いた箱'))
import fold_export as fe

def area(poly): return abs(sum(poly[i][0]*poly[(i+1)%len(poly)][1]-poly[(i+1)%len(poly)][0]*poly[i][1] for i in range(len(poly))))/2
def pts(P3): return np.array([v for poly in P3 for v in poly])
def dmat(X): return np.linalg.norm(X[:,None]-X[None],axis=2)

for f in sorted(glob.glob(os.path.join(ROOT,'アプリ','data','data_tetra_*.json'))):
    d=json.load(open(f))
    for t in d['targets']:
        th={c:t['theta'][c] for c in d['order'][1:]}
        P3=fe.fold3d(d['panels'],d['parent'],d['hinge'],d['order'],th); X=pts(P3)
        if not t['label'].startswith('等面四面体'):   # 箱: fold_export が出した target3d と相似
            D1,D2=dmat(X),dmat(pts(t['target3d'])); k=(D1*D2).sum()/(D2*D2).sum()
            assert np.abs(D1-k*D2).max()<1e-6,(t['label'],k); continue
        G={}   # パネルの平面ごと（色 = 面）にまとめた平面の方程式（法線，オフセット）。重心が内側になる向きにそろえる
        for poly,c in zip(P3,t['colors']):
            p0,p1,p2=map(np.array,poly[:3]); n=np.cross(p1-p0,p2-p0); n/=np.linalg.norm(n); G.setdefault(c,[]).append((*n,-n@p0))
        assert len(G)==4
        E=[]
        for g in G.values():
            e=np.array(g); assert np.abs(e-e[0]).max()<1e-6 or np.abs(e+e[0]).max()<1e-6   # 同じ面の中は同一平面
            E.append(e[0]*(1 if e[0][:3]@X.mean(0)+e[0][3]<0 else -1))
        E=np.array(E); corners=np.array([np.linalg.solve(E[list(i),:3],-E[list(i),3]) for i in ((0,1,2),(0,1,3),(0,2,3),(1,2,3))])
        S=np.einsum('fk,pk->fp',E[:,:3],X)+E[:,3:]       # 全頂点: どの面にも内側か上（凸包の内側），かつどれかの面の上（表面上）
        assert S.max()<1e-6 and np.abs(S).min(axis=0).max()<1e-6
        # 面ごとの面積: 色（= 面）ごとのパネル面積の和が四面体の面積（紙の面積/4）に一致
        tot=sum(area(p['poly']) for p in d['panels']); A={}
        for p,c in zip(d['panels'],t['colors']): A[c]=A.get(c,0)+area(p['poly'])
        assert len(A)==4 and all(abs(a-tot/4)<1e-6 for a in A.values()),A
        # 表面積 = 紙の面積（重なり・すき間なし）
        hull=ConvexHull(corners).area; assert abs(hull-tot)<1e-6,(hull,tot)
        e=sorted(round(float(np.sum((a-b)**2)),6) for i,a in enumerate(corners) for b in corners[i+1:])
        want=[int(x[1:]) for x in re.findall(r'√\d+',t['label'])]; got=sorted({round(x/e[0]*want[0]) for x in e})
        assert got==sorted(set(want)),(got,want)
        assert e[0]==e[1] and e[2]==e[3] and e[4]==e[5]
    print('OK',os.path.basename(f),len(d["panels"]),"panels",len(d['targets']),'targets')
