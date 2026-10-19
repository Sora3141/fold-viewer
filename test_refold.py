#!/usr/bin/env python3
"""data/data_refold_*.json を全ヒンジ目標角で折ったとき，正多面体（面・頂点の数，辺 1）と等面四面体（4 面・4 頂点）に閉じるかを確かめる。
同じ名前の四面体（X* など）が 2 枚の紙で相似なことも確かめる。usage: python3 test_refold.py   （研究リポジトリの中で実行）"""
import os, json, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
from make_dodeca import fe

def shape(d, t):
    th = {c: t['theta'][c] for c in d['order'][1:]}; P3 = fe.fold3d(d['panels'], d['parent'], d['hinge'], d['order'], th); planes = []; X = []
    for poly in P3:
        p0 = np.array(poly[0]); n = sum(np.cross(np.array(poly[i]) - p0, np.array(poly[(i + 1) % len(poly)]) - p0) for i in range(len(poly))); n /= np.linalg.norm(n); q = np.array([*n, n @ p0])
        if all(np.abs(q - k).max() > 1e-6 for k in planes): planes.append(q)
        for v in map(np.array, poly):
            if all(np.linalg.norm(v - u) > 1e-7 for u in X): X.append(v)
    for q in planes: assert (np.array(X) @ q[:3] - q[3]).max() < 1e-6 or (np.array(X) @ q[:3] - q[3]).min() > -1e-6   # 凸
    V = [v for v in X if sum(abs(v @ q[:3] - q[3]) < 1e-7 for q in planes) >= 3]
    D = sorted(np.linalg.norm(V[i] - V[j]) for i in range(len(V)) for j in range(i + 1, len(V)))
    return len(planes), len(V), D

ratios = {}
for f in sorted(glob.glob(os.path.join(HERE, 'data', 'data_refold_*_X*.json'))):   # 1 枚の紙のデータ（つないだものは除く）
    d = json.load(open(f)); want = {'立方体': (6, 8), '正二十面体': (20, 12), '正八面体': (8, 6)}
    for t in d['targets']:
        F, V, D = shape(d, t); assert (F, V) == want.get(t['label'], (4, 4)), (f, t['label'], F, V)
        if V == 4: ratios.setdefault(t['label'], []).append(np.array(D[::2]) / D[0])
        else: assert abs(D[0] - 1) < 1e-6, (f, t['label'], D[0])
        print(os.path.basename(f), t['label'], F, '面', V, '頂点', 'OK')
for lab, rs in ratios.items():
    assert len(rs) == 2 and np.abs(rs[0] - rs[1]).max() < 1e-6, (lab, rs)   # 2 枚の紙の四面体が相似
    print(lab, '辺の比', np.round(rs[0], 6), '（2 枚で一致）')
