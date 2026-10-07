#!/usr/bin/env python3
"""data/data_dodeca_*.json を全ヒンジ目標角で折ったとき、12 面体の再構成（十二面体/figures/area10_ex2_creases.json の X）と合うかを確かめる。
 - 8 頂点（紙の格子点）の 3D 位置が X と剛体変換（回転と平行移動のみ）で一致  - 全パネルの頂点が多面体の表面上
 - 同じ色（面）のパネルは同一平面、面は 12  - 面ごとの面積の和 = 三角形の面積  - 表面積 = 紙の面積 = 10
 - 箱は fold_export の target3d と相似
usage: python3 test_dodeca.py   （研究リポジトリの中で実行）"""
import sys, os, json, glob
import numpy as np
from scipy.spatial import ConvexHull
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, '開いた箱', 'src')); os.chdir(os.path.join(ROOT, '開いた箱'))
import fold_export as fe

def area(poly): return abs(sum(poly[i][0] * poly[(i + 1) % len(poly)][1] - poly[(i + 1) % len(poly)][0] * poly[i][1] for i in range(len(poly)))) / 2
def pts(P3): return np.array([v for poly in P3 for v in poly])
def dmat(X): return np.linalg.norm(X[:, None] - X[None], axis=2)

def kabsch(A, B):
    """A を B に重ねる回転（det=+1）と平行移動の残差"""
    a, b = A.mean(0), B.mean(0); U, _, Vt = np.linalg.svd((A - a).T @ (B - b)); D = np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))])
    R = U @ D @ Vt; return np.abs((A - a) @ R + b - B).max(), np.linalg.det(U @ Vt)

C = json.load(open(os.path.join(ROOT, '十二面体', 'figures', 'area10_ex2_creases.json')))
# 結果 2 の例（area10_ex1.json、別のポリオミノ）と同じ多面体か: 頂点間距離の多重集合（頂点の番号・鏡映によらない）と体積
_ex1 = np.array(json.load(open(os.path.join(ROOT, '十二面体', 'figures', 'area10_ex1.json')))['X'])
for _g in C['gluings']:
    _X = np.array(_g['X']); assert np.abs(np.sort(dmat(_X)[np.triu_indices(8, 1)]) - np.sort(dmat(_ex1)[np.triu_indices(8, 1)])).max() < 1e-6
    assert abs(ConvexHull(_X).volume - ConvexHull(_ex1).volume) < 1e-6
for f in sorted(glob.glob(os.path.join(ROOT, 'アプリ', 'data', 'data_dodeca_*.json'))):
    d = json.load(open(f)); tot = sum(area(p['poly']) for p in d['panels']); assert abs(tot - 10) < 1e-9, tot
    j = 0
    for t in d['targets']:
        th = {c: t['theta'][c] for c in d['order'][1:]}
        P3 = fe.fold3d(d['panels'], d['parent'], d['hinge'], d['order'], th); Y = pts(P3)
        if not t['label'].startswith('12 面体'):
            D1, D2 = dmat(Y), dmat(pts(t['target3d'])); k = (D1 * D2).sum() / (D2 * D2).sum()
            assert np.abs(D1 - k * D2).max() < 1e-6, (t['label'], k); continue
        g = C['gluings'][j]; j += 1; X = np.array(g['X']); faces = [tuple(x) for x in g['faces']]
        # 1) 頂点: 紙の格子点 vpt[u] に来るパネルの頂点の 3D 位置（全パネルで一致する）と X
        Z = []
        for u, (px, py) in enumerate(g['vpt']):
            q = [np.array(P3[i][k]) for i, p in enumerate(d['panels']) for k, v in enumerate(p['poly']) if abs(v[0] - px) < 1e-9 and abs(v[1] - py) < 1e-9]
            assert q and max(np.abs(x - q[0]).max() for x in q) < 1e-9, (u, px, py); Z.append(q[0])   # 同じ頂点に集まる点は全部同じ位置
        res, det = kabsch(np.array(Z), X); assert det > 0 and res < 1e-6, (t['label'], res, det)
        # 2) 全頂点が多面体の表面上（凸包の内側かつどれかの面の上）
        E = []
        for a, b, c in faces:
            n = np.cross(X[b] - X[a], X[c] - X[a]); n /= np.linalg.norm(n); E.append((*n, -n @ X[a]))
        E = np.array(E)    # 外向き（重心が全部の面の内側）
        assert (np.einsum('fk,k->f', E[:, :3], X.mean(0)) + E[:, 3]).max() < 0
        # 折った形を X に重ねてから測る（剛体変換は 1) の Kabsch）
        a_, b_ = np.array(Z).mean(0), X.mean(0); U, _, Vt = np.linalg.svd((np.array(Z) - a_).T @ (X - b_)); Rm = U @ Vt
        Ym = (Y - a_) @ Rm + b_
        S = np.einsum('pk,fk->pf', Ym, E[:, :3]) + E[:, 3]   # (@ は macOS の numpy で警告が出るため)
        assert S.max() < 1e-6 and np.abs(S).min(axis=1).max() < 1e-6, (S.max(), np.abs(S).min(axis=1).max())
        # 3) 同じ色（面）のパネルは同一平面。色は 12 種
        G = {}
        for poly, c in zip(P3, t['colors']):
            p0, p1, p2 = map(np.array, poly[:3]); n = np.cross(p1 - p0, p2 - p0); n /= np.linalg.norm(n); G.setdefault(c, []).append((*n, -n @ p0))
        assert len(G) == 12
        for e in map(np.array, G.values()): assert np.abs(e - e[0]).max() < 1e-6, e
        # 4) 面ごとの面積: 色ごとのパネル面積の和が、その面（三角形）の面積と一致。面の対応は法線で取る
        A = {}
        for p, c in zip(d['panels'], t['colors']): A[c] = A.get(c, 0) + area(p['poly'])
        tri = [0.5 * np.linalg.norm(np.cross(X[b] - X[a], X[c] - X[a])) for a, b, c in faces]
        for c, e in G.items():
            i = int(np.argmax(np.einsum('fk,k->f', E[:, :3], e[0][:3] @ Rm))); assert abs(E[i, :3] @ (e[0][:3] @ Rm) - 1) < 1e-6, (c, i)
            assert abs(A[c] - tri[i]) < 1e-9, (c, i, A[c], tri[i])
        # 5) 表面積 = 紙の面積 = 10（重なり・すき間なし）
        assert abs(sum(tri) - 10) < 1e-9 and abs(ConvexHull(X).area - 10) < 1e-9
        print('  ', t['label'], 'kabsch残差 %.1e' % res, '面積 10 面ごと一致')
    print('OK', os.path.basename(f), len(d['panels']), 'panels', len(d['targets']), 'targets')
