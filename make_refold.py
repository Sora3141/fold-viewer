#!/usr/bin/env python3
"""再折り: 立方体と正二十面体の辺展開図が，同じ等面四面体 X* にも折れるデータを作る → data/data_refold_*.json
入力: 再折り/figures/cube_ico_X.json（再折り/src/export_cube_ico.py の出力: 周の頂点・回転中心・回転中心の格子 Λ0 の基底）。
等面四面体の辺は，回転中心の格子 Λ0 のドロネー三角形の辺の直線（簡約基底 ω1, ω2 と ω1−ω2 の 3 方向）。
パネル＝紙を水平な格子線で台形に分け、正多面体の格子線（M=4 は正方格子、M=6 は三角格子）と四面体の直線で切った凸多角形。
正多面体の折り線は格子線（辺展開図も、正八面体の円環の展開図のような一般の展開図も、紙の中の格子線はすべて正多面体の辺）。
折り角は「2 つの面の外向き法線のなす角」，ほかの境目は 0（make_multi.py と同じ向き）。
usage: python3 make_refold.py   （研究リポジトリの中で実行）"""
import sys, os, json, math
from collections import deque
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
from make_dodeca import fe                         # 開いた箱/src を path に入れ、cwd を移す
from make_multi import palette, on_seg

SRC = os.path.join(ROOT, '再折り', 'figures', 'cube_ico_X.json')

def cr(a, b): return a[0] * b[1] - a[1] * b[0]

def inside(z, poly):
    x, y = z; c = False
    for i in range(len(poly)):
        (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % len(poly)]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1: c = not c
    return c

def grid(M):
    """正多面体の格子線の向き d と、d と組んで格子の基底になる o（どちらも原点の格子点を通る）"""
    if M == 4: return [(np.array([1., 0.]), np.array([0., 1.])), (np.array([0., 1.]), np.array([1., 0.]))]
    z = np.array([.5, math.sqrt(3) / 2]); one = np.array([1., 0.])
    return [(one, z), (z, one), (z - one, one)]

def trapezoids(P, h):
    """周 P（頂点はどれも高さ h の倍数の水平線の上）を水平線で切った台形（反時計回り）"""
    lv = [round(y / h) for x, y in P]; assert all(abs(y / h - k) < 1e-9 for (x, y), k in zip(P, lv)), '頂点が水平な格子線の上にない'
    out = []
    for k in range(min(lv), max(lv)):
        y0, y1, ym = k * h, (k + 1) * h, (k + .5) * h; X = []
        for i in range(len(P)):
            (x0, a0), (x1, a1) = P[i], P[(i + 1) % len(P)]
            if (a0 - ym) * (a1 - ym) < 0: f = lambda y: x0 + (x1 - x0) * (y - a0) / (a1 - a0); X.append((f(ym), f(y0), f(y1)))
        X.sort()
        for (_, l0, l1), (_, r0, r1) in zip(X[::2], X[1::2]):
            q = [(l0, y0), (r0, y0), (r1, y1), (l1, y1)]
            out.append([v for i, v in enumerate(q) if abs(v[0] - q[i - 1][0]) > 1e-12 or abs(v[1] - q[i - 1][1]) > 1e-12])
    return out

def reduce(e1, e2):
    """ガウスの簡約。返すのはドロネー三角形の 3 辺の向き（鋭角三角形）"""
    u, v = np.array(e1), np.array(e2)
    while True:
        if u @ u > v @ v: u, v = v, u
        m = round(float(u @ v / (u @ u)))
        if m == 0: break
        v = v - m * u
    if u @ v < 0: v = -v
    return u, v, v - u

def tetra_dihedral(L):
    """辺の長さ L（3 つ，向かい合う辺は等しい）の等面四面体の，辺の長さ → 内側の二面角。直方体 p×q×r の対角の 4 頂点で作る"""
    a, b, c = (l * l for l in L)
    p, q, r = (math.sqrt(v) for v in ((a + c - b) / 2, (a + b - c) / 2, (b + c - a) / 2))
    V = np.array([(0, 0, 0), (p, q, 0), (p, 0, r), (0, q, r)], float); out = []
    for i in range(4):
        for j in range(i + 1, 4):
            e = V[j] - V[i]; ln = np.linalg.norm(e); e /= ln; z, w = [V[k] - V[i] for k in range(4) if k not in (i, j)]
            z -= z @ e * e; w -= w @ e * e; out.append((ln, math.acos(z @ w / np.linalg.norm(z) / np.linalg.norm(w))))
    return lambda l: min(out, key=lambda t: abs(t[0] - l))[1]

def make(net, out, title, solid, th_solid, nfaces, xname='X*'):
    P = [tuple(p) for p in net['poly']]; M = net['M']; G = grid(M)
    cells = trapezoids(P, 1.0 if M == 4 else math.sqrt(3) / 2); c1 = np.array(net['centers'][0]); dirs = reduce(net['e1'], net['e2'])
    others = [dirs[1], dirs[0], dirs[0]]                  # d と組んで格子の基底になるもの（ω3 = ω2 − ω1 と ω1）
    dih = tetra_dihedral([float(np.linalg.norm(d)) for d in dirs])
    pts = np.array(P, float)
    lines = []                                             # 切る直線: 正多面体の格子線と等面四面体の辺の直線（2 点）
    for (d, o), base in [(g, np.zeros(2)) for g in G] + [((d, o), c1) for d, o in zip(dirs, others)]:
        ks = [cr(d, p - base) / cr(d, o) for p in pts]
        for k in range(math.floor(min(ks)) - 1, math.ceil(max(ks)) + 2):
            a = base + k * o; lines.append((tuple(a), tuple(a + d)))
    panels = []
    for C in cells:
        polys = [C]
        for a, b in lines: polys = [q for pg in polys for q in fe.split_convex(pg, a, b)]
        for pg in polys:
            ar = sum(cr(pg[i], pg[(i + 1) % len(pg)]) for i in range(len(pg))) / 2
            if abs(ar) < 1e-10: continue
            while abs(cr((pg[1][0] - pg[0][0], pg[1][1] - pg[0][1]), (pg[2][0] - pg[0][0], pg[2][1] - pg[0][1]))) < 1e-9: pg = pg[1:] + pg[:1]   # ビューアは先頭 3 点で向きを取る
            panels.append({'poly': [[float(x), float(y)] for x, y in pg]})
    E = [[(tuple(p['poly'][k]), tuple(p['poly'][(k + 1) % len(p['poly'])])) for k in range(len(p['poly']))] for p in panels]
    adj = {}
    for i, ei in enumerate(E):
        for e in ei:
            m = ((e[0][0] + e[1][0]) / 2, (e[0][1] + e[1][1]) / 2)
            for o, eo in enumerate(E):
                if o == i: continue
                for f in eo:
                    ein = on_seg(e[0], *f) and on_seg(e[1], *f); fin = on_seg(f[0], *e) and on_seg(f[1], *e)
                    if on_seg(m, *f) and (ein or fin): adj.setdefault(i, []).append((o, e if ein else f)); break
    parent = {0: None}; hinge = {}; order = [0]; dq = deque([0])
    while dq:
        p = dq.popleft()
        for c, e in adj.get(p, []):
            if c not in parent: parent[c] = p; hinge[c] = sorted([list(e[0]), list(e[1])]); order.append(c); dq.append(c)
    assert len(parent) == len(panels), 'panel graph not connected'
    def sgn(c):
        (h0, h1), poly = hinge[c], panels[c]['poly']
        cx = sum(p[0] for p in poly) / len(poly); cy = sum(p[1] for p in poly) / len(poly)
        return -1 if (h1[0] - h0[0]) * (cy - h0[1]) - (h1[1] - h0[1]) * (cx - h0[0]) > 0 else 1   # 紙の表を外側にして折る
    def on_line(c, fam, base):
        h0, h1 = map(np.array, hinge[c]); m = (h0 + h1) / 2
        for d, o in fam:
            if abs(cr(d, h1 - h0)) > 1e-7: continue
            off = cr(d, m - base) / cr(d, o)
            if abs(off - round(off)) < 1e-7: return d
        return None
    def th_tetra(c):
        d = on_line(c, list(zip(dirs, others)), c1)
        return 0.0 if d is None else sgn(c) * (math.pi - dih(float(np.linalg.norm(d))))
    thetas = [{c: (sgn(c) * th_solid if on_line(c, G, np.zeros(2)) is not None else 0.0) for c in order[1:]},
              {c: th_tetra(c) for c in order[1:]}]
    targets = []
    for lab, th, nf in zip([solid, '等面四面体 ' + xname], thetas, [nfaces, 4]):
        P3 = fe.fold3d(panels, parent, hinge, order, th); key = []; cid = []
        for poly in P3:
            p0 = np.array(poly[0]); n = sum(np.cross(np.array(poly[i]) - p0, np.array(poly[(i + 1) % len(poly)]) - p0) for i in range(len(poly))); n /= np.linalg.norm(n)
            q = np.array([*n, n @ p0]); i = next((i for i, k in enumerate(key) if np.abs(k - q).max() < 1e-6), None)
            if i is None: key.append(q); i = len(key) - 1
            cid.append(i)
        assert len(key) == nf, (lab, len(key), nf)
        X = np.array([v for poly in P3 for v in poly]); ctr = X.mean(0)
        for q in key: assert (X @ q[:3] - q[3]).max() < 1e-6 or (X @ q[:3] - q[3]).min() > -1e-6, (lab, '凸でない')   # どの面の平面も全体の片側
        pal = palette(nf); targets.append({'label': lab, 'theta': [th.get(i, 0.0) for i in range(len(panels))], 'colors': [pal[c] for c in cid], 'target3d': P3})
    cells_out = [[sum(p[0] for p in q['poly']) / len(q['poly']) - .5, sum(p[1] for p in q['poly']) / len(q['poly']) - .5] for q in panels]   # ビューアの中心と大きさ合わせ用（面積の表示は apps.json の area_label）
    data = {'title': title, 'cells': cells_out, 'panels': panels, 'parent': [parent[i] if parent[i] is not None else -1 for i in range(len(panels))],
            'hinge': [hinge.get(i, [[0, 0], [0, 0]]) for i in range(len(panels))], 'order': order, 'targets': targets}
    json.dump(data, open(out, 'w')); print(out, len(cells), 'cells', len(panels), 'panels', [t['label'] for t in targets])

def chain(fa, fb, out, xname='X*'):
    """2 つのデータ（紙 A と紙 B。どちらも目標 1 が同じ等面四面体）をつなぐ。紙 B を面積がそろうよう拡大・縮小し、
    B の四面体を A の四面体にぴったり重ねる相似変換（鏡映もありうる）align を求め、B の四面体の面の色を A にそろえる。"""
    import itertools
    A, B = json.load(open(fa)), json.load(open(fb))
    def area(d): return sum(abs(sum(cr(p['poly'][i], p['poly'][(i + 1) % len(p['poly'])]) for i in range(len(p['poly'])))) / 2 for p in d['panels'])
    def tetra(d):
        P3 = d['targets'][1]['target3d']; X = np.array([v for p in P3 for v in p]); planes = []
        for poly in P3:
            p0 = np.array(poly[0]); n = sum(np.cross(np.array(poly[i]) - p0, np.array(poly[(i + 1) % len(poly)]) - p0) for i in range(len(poly))); n /= np.linalg.norm(n)
            if all(np.abs(np.array([*n, n @ p0]) - k).max() > 1e-6 for k in planes): planes.append(np.array([*n, n @ p0]))
        V = []
        for v in X:
            if sum(abs(v @ q[:3] - q[3]) < 1e-7 for q in planes) >= 3 and all(np.linalg.norm(v - u) > 1e-7 for u in V): V.append(v)
        return np.array(V)
    VA, VB = tetra(A), tetra(B); assert len(VA) == len(VB) == 4
    s = math.sqrt(area(A) / area(B)); best = None
    for perm in itertools.permutations(range(4)):
        P, Q = VB[list(perm)] * s, VA; pc, qc = P.mean(0), Q.mean(0)
        U, S_, Wt = np.linalg.svd((P - pc).T @ (Q - qc)); R = (U @ Wt).T          # R (P - pc) ≈ Q - qc（鏡映も許す）
        err = np.abs((R @ (P - pc).T).T - (Q - qc)).max()
        if best is None or err < best[0]: best = (err, R, pc, qc)
    err, R, pc, qc = best; assert err < 1e-6, err
    Mx = np.eye(4); Mx[:3, :3] = R * s; Mx[:3, 3] = qc - R @ pc                   # world = R (s x) + (qc − R pc)
    def planes_cols(P3, cols, M):
        out = []
        for poly, c in zip(P3, cols):
            X = np.array([(M @ np.array([*v, 1]))[:3] for v in poly]); n = np.cross(X[1] - X[0], X[2] - X[0]); n /= np.linalg.norm(n); out.append((n, n @ X[0], c))
        return out
    pa = planes_cols(A['targets'][1]['target3d'], A['targets'][1]['colors'], np.eye(4)); newc = []
    for n, o, _ in planes_cols(B['targets'][1]['target3d'], B['targets'][1]['colors'], Mx):   # B の面の色を、重ねた位置で同じ平面にある A の面の色に
        newc.append(next(c for m, p, c in pa if (np.abs(m - n).max() < 1e-6 and abs(p - o) < 1e-6) or (np.abs(m + n).max() < 1e-6 and abs(p + o) < 1e-6)))
    B['targets'][1]['colors'] = newc
    mirror = bool(np.linalg.det(R) < 0)
    A['chain'] = {'xname': xname, 'net': B, 'align': Mx.T.flatten().tolist(), 'mirror': mirror, 'scale': s}   # three.js の Matrix4.fromArray は列優先
    json.dump(A, open(out, 'w')); print(out, '倍率', round(s, 6), '鏡映' if mirror else '鏡映なし', '重ねの誤差', err)

CUBE, ICO, OCT = ('立方体', math.pi / 2, 6), ('正二十面体', math.pi - math.acos(-math.sqrt(5) / 3), 20), ('正八面体', math.pi - math.acos(-1 / 3), 8)

if __name__ == '__main__':
    F = os.path.join(ROOT, '再折り', 'figures'); D = os.path.join(HERE, 'data'); d = lambda n: os.path.join(D, f'data_refold_{n}.json')
    S = json.load(open(os.path.join(F, 'cube_ico_X.json')))
    make(S['cube'], d('cube_X'), '立方体の辺展開図 → 等面四面体 X*', *CUBE)
    make(S['ico'], d('ico_X'), '正二十面体の辺展開図 → 等面四面体 X*', *ICO)
    chain(d('cube_X'), d('ico_X'), d('cube_ico'))
    S = json.load(open(os.path.join(F, 'oct_cube_X.json')))      # 正八面体は円環を切った一般の展開図（辺展開図ではない）
    make(S['oct'], d('oct_X8'), '正八面体の展開図 → 等面四面体 X₈', *OCT, xname='X₈')
    make(S['cube'], d('cube_X8'), '立方体の辺展開図 → 等面四面体 X₈', *CUBE, xname='X₈')
    chain(d('oct_X8'), d('cube_X8'), d('oct_cube'), 'X₈')
    S = json.load(open(os.path.join(F, 'oct_ico_X.json')))
    make(S['oct'], d('oct_X82'), '正八面体の展開図 → 等面四面体 X₈₂', *OCT, xname='X₈₂')
    make(S['ico'], d('ico_X82'), '正二十面体の辺展開図 → 等面四面体 X₈₂', *ICO, xname='X₈₂')
    chain(d('oct_X82'), d('ico_X82'), d('oct_ico'), 'X₈₂')
