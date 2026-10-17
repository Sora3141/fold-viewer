#!/usr/bin/env python3
"""data/data_multi_*.json（make_multi.py）を全ヒンジ目標角で折ったとき、多面体（多面数/figures/ex?_creases.json の X）と合うかを確かめる。
 - 紙の頂点（円錐点）の 3D 位置が X と剛体変換（回転と平行移動のみ）で一致  - 全パネルの頂点が多面体の表面上
 - 同じ色（面）のパネルは同一平面で、面の数は多面体の面の数  - 面ごとの面積の和 = 面の面積  - 表面積 = 紙の面積
 - 周の貼り合わせ（pairs）の両端が 3D で重なる（紙が閉じる）
usage: python3 test_multi.py   （研究リポジトリの中で実行）"""
import sys, os, json
import numpy as np
from scipy.spatial import ConvexHull
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, '多面数', 'src')); sys.path.insert(0, os.path.join(ROOT, '開いた箱', 'src')); os.chdir(os.path.join(ROOT, '開いた箱'))
import fold_export as fe, creases_multi as cm

def area(poly): return abs(sum(poly[i][0] * poly[(i + 1) % len(poly)][1] - poly[(i + 1) % len(poly)][0] * poly[i][1] for i in range(len(poly)))) / 2
def rk(p): return (round(p[0], 6) + 0.0, round(p[1], 6) + 0.0)

for name in ('poly_a10', 'ico_5solids', 'ico_9solids', 'ico_11solids', 'tri22_13solids'):
    d = json.load(open(os.path.join(ROOT, 'アプリ', 'data', f'data_multi_{name}.json')))
    C = json.load(open(os.path.join(ROOT, '多面数', 'figures', {'poly_a10': 'ex1', 'ico_5solids': 'ex2', 'ico_9solids': 'ex2_9', 'ico_11solids': 'ex_11', 'tri22_13solids': 'ex_13'}[name] + '_creases.json'))); keys = list(C['gluings'])
    tot = sum(area(p['poly']) for p in d['panels']); assert abs(tot - (len(C['cells']) if C['kind'] == 'sq' else len(C['cells']) * 3 ** .5 / 4)) < 1e-9, tot
    assert len(d['targets']) == len(keys); worst = 0
    for t, k in zip(d['targets'], keys):
        g = C['gluings'][k]; X = np.array(g['X']); faces, edges, normals = cm.polyhedron(X)
        th = {c: t['theta'][c] for c in d['order'][1:]}; P3 = fe.fold3d(d['panels'], d['parent'], d['hinge'], d['order'], th)
        img = {}   # 紙の点 → 3D（同じ点に集まるパネルの頂点は全部同じ位置）
        for p, P in zip(d['panels'], P3):
            for v, q in zip(p['poly'], P):
                q = np.array(q); r = img.setdefault(rk(v), q); assert np.abs(r - q).max() < 1e-9, (k, v)
        # 1) 円錐点: 紙の頂点の 3D 位置と X が剛体変換（回転と平行移動のみ）で一致
        Z = np.array([img[rk(v)] for v in g['vpt']]); a, b = Z.mean(0), X.mean(0)
        U, _, Vt = np.linalg.svd((Z - a).T @ (X - b)); D = np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))]); R = U @ D @ Vt
        res = np.abs((Z - a) @ R + b - X).max(); assert np.linalg.det(U @ Vt) > 0 and res < 1e-6, (k, res); worst = max(worst, res)
        # 2) 全頂点が多面体の表面上
        E = np.array([(*n, -n @ X[f[0]]) for n, f in zip(normals, faces)]); Ym = (np.array([q for P in P3 for q in P]) - a) @ R + b
        S = np.einsum('pk,fk->pf', Ym, E[:, :3]) + E[:, 3]; assert S.max() < 1e-6 and np.abs(S).min(axis=1).max() < 1e-6, (k, S.max())
        # 3) 同じ色（面）のパネルは同一平面。面の数 = 多面体の面の数
        Gp = {}
        for P, c in zip(P3, t['colors']):
            p0 = np.array(P[0]); n = sum(np.cross(np.array(P[i]) - p0, np.array(P[(i + 1) % len(P)]) - p0) for i in range(len(P))); n /= np.linalg.norm(n); Gp.setdefault(c, []).append((*n, -n @ p0))
        assert len(Gp) == len(faces), (k, len(Gp))
        for e in map(np.array, Gp.values()): assert np.abs(e - e[0]).max() < 1e-6
        # 4) 面ごとの面積: 色ごとのパネル面積の和が、その面の面積と一致（面の対応は法線で取る）
        A = {}
        for p, c in zip(d['panels'], t['colors']): A[c] = A.get(c, 0) + area(p['poly'])
        h = ConvexHull(X); fa = [sum(0.5 * np.linalg.norm(np.cross(X[s[1]] - X[s[0]], X[s[2]] - X[s[0]])) for s, q in zip(h.simplices, h.equations) if np.dot(q[:3], n) > 1 - 1e-9) for n in normals]
        for c, e in Gp.items():
            n = e[0][:3] @ R; i = int(np.argmax(normals @ n)); assert abs(normals[i] @ n - 1) < 1e-6 and abs(A[c] - fa[i]) < 1e-9, (k, c, i, A[c], fa[i])
        # 5) 表面積 = 紙の面積
        assert abs(h.area - tot) < 1e-9, (k, h.area, tot)
        # 6) 周の貼り合わせ: 半線分 i の始点・終点が、相手 j の終点・始点と 3D で重なる
        pp = cm.Paper(C['kind'], C['cells'], g['pairs'])
        for i, j in g['pairs']:
            for s_, t_ in ((i, j + 1), (i + 1, j)): assert np.abs(img[rk(pp.perim_pos(s_))] - img[rk(pp.perim_pos(t_))]).max() < 1e-9, (k, i, j)
        print('  ', t['label'], 'kabsch残差 %.1e' % res, '面 %d 面ごとの面積一致・貼り合わせ %d 組が閉じる' % (len(faces), len(g['pairs'])))
    print('OK', name, len(d['panels']), 'panels', len(d['targets']), 'targets, 最大残差 %.1e' % worst)
