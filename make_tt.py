#!/usr/bin/env python3
"""再折り: 切頂四面体（辺 1）の格子の木の展開図が、等面四面体（いまは辺 √7 の正四面体）にも折れるデータを作る → data/data_refold_tt_*.json
入力: 再折り/work/ar/truncated_tetrahedron_lattice.json と、その当たり（ar_conway.py の rand=100000 → tt_lattice_big.json.hits.json）。
切頂四面体の折り角: 木で切って開いた単位三角形 28 枚（六角形は中心で 6 つに割る）を自分で並べ、展開図の周（ar_conway.word）に重ねる。
  折り目の両側の三角形が同じ面なら 0、違う面なら 2 つの面の外向き法線のなす角（六角形どうし・三角形と六角形で違う）。
等面四面体の折り目と角は make_refold.make と同じ（Conway の 4 つの中心の格子のドロネー三角形の辺）。
usage: python3 make_tt.py   （研究リポジトリの中で実行）"""
import sys, os, json, math, cmath, itertools
from collections import deque
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, '再折り', 'src'); WORK = os.path.join(ROOT, '再折り', 'work', 'ar')
sys.path.insert(0, HERE); sys.path.insert(0, SRC)
import make_refold as MR
import ar_geom as AG
GEOM = os.path.join(WORK, 'truncated_tetrahedron_lattice.json')
argv = sys.argv; sys.argv = [argv[0], GEOM]
import ar_conway as C                    # 読み込むときに sys.argv[1] の json を読む
sys.argv = argv
P = AG.G.Poly('truncated_tetrahedron')


def unit_triangles():
    """(頂点の番号 3 つ, 面の番号)。六角形の中心の番号は ar_geom と同じ（六角形の面の順に nreal から）"""
    X = {v: np.array(P.V[v], float) for v in range(P.nv)}; tri = []; c = P.nv
    for fi, f in enumerate(P.faces):
        if len(f) == 3: tri.append((tuple(f), fi)); continue
        X[c] = np.mean([X[v] for v in f], 0)
        for i in range(6): tri.append(((f[i], f[(i + 1) % 6], c), fi))
        c += 1
    return X, tri


def unfold(tree):
    """木で切って、単位三角形を平面に並べる（三角形ごとの 2 次元の頂点）"""
    X, tri = unit_triangles(); cut = {frozenset(C.CUTS[k][:2]) for k in tree}
    for k in tree: assert abs(np.linalg.norm(X[C.CUTS[k][0]] - X[C.CUTS[k][1]]) - C.CUTS[k][2]) < 1e-9   # 番号が ar_geom と合う
    ed = {}
    for t, (vs, _) in enumerate(tri):
        for a, b in itertools.combinations(vs, 2): ed.setdefault(frozenset((a, b)), []).append(t)
    assert all(len(v) == 2 for v in ed.values()) and len(ed) == len(C.CUTS)
    a, b, c = tri[0][0]; d = lambda u, w: float(np.linalg.norm(X[u] - X[w]))
    def third(pa, pb, la, lb, away):
        L = abs(pb - pa); x = (la * la - lb * lb + L * L) / (2 * L); y = math.sqrt(max(la * la - x * x, 0)); e = (pb - pa) / L
        q1, q2 = pa + e * complex(x, y), pa + e * complex(x, -y)
        return q1 if abs(q1 - away) > abs(q2 - away) else q2
    pos = {0: {a: 0j, b: complex(d(a, b), 0)}}; pos[0][c] = third(0j, pos[0][b], d(a, c), d(b, c), complex(0, -1))
    dq = deque([0])
    while dq:
        t = dq.popleft(); vs = tri[t][0]
        for u, w in itertools.combinations(vs, 2):
            if frozenset((u, w)) in cut: continue
            s = [x for x in ed[frozenset((u, w))] if x != t][0]
            if s in pos: continue
            o = [x for x in tri[s][0] if x not in (u, w)][0]; mine = [x for x in vs if x not in (u, w)][0]
            pos[s] = {u: pos[t][u], w: pos[t][w], o: third(pos[t][u], pos[t][w], d(u, o), d(w, o), pos[t][mine])}
            dq.append(s)
    assert len(pos) == len(tri)
    return X, tri, pos


def align(pos, pts):
    """自分で並べた三角形を、ar_conway.word の周の座標に重ねる合同変換（裏返しも試す）"""
    Q = [complex(round(p.real, 9), round(p.imag, 9)) for p in pts]
    near = lambda z: min(abs(z - q) for q in Q)
    allp = [p for t in pos.values() for p in t.values()]
    a, b = list(pos[0].values())[:2]
    for fl in (False, True):
        f = (lambda z: z.conjugate()) if fl else (lambda z: z)
        for p in Q:
            for q in Q:
                if abs(abs(q - p) - abs(b - a)) > 1e-9: continue
                r = (q - p) / (f(b) - f(a)); T = lambda z, r=r, p=p: p + r * (f(z) - f(a))
                if all(min(abs(T(z) - x) for z in allp) < 1e-7 for x in Q): return T
    raise SystemExit('周に重ならない')


def build(idx, out, title, tetra_name):
    h = json.load(open(os.path.join(WORK, 'tt_lattice_big.json.hits.json')))[idx]
    syms, angs, z, pts = C.word(h['tree'])
    n = len(pts); P2 = []
    for k in range(n): P2 += [pts[k], (pts[k] + pts[(k + 1) % n]) / 2]
    L = len(P2); cu = h['cuts'][0]; c = list(cu) + [cu[0] + L]
    cen = [(P2[c[k] % L] + P2[c[k + 1] % L]) / 2 for k in (1, 2, 4, 5)]          # ar_figure.py と同じ: 点対称な 4 区間の中心
    X, tri, pos = unfold(h['tree']); T = align(pos, pts)
    tri2 = [([T(pos[t][v]) for v in vs], fi) for t, (vs, fi) in enumerate(tri)]
    ctr = np.mean([X[v] for v in range(P.nv)], 0)
    def normal(fi):
        f = P.faces[fi]; nn = np.cross(X[f[1]] - X[f[0]], X[f[2]] - X[f[0]]); nn /= np.linalg.norm(nn)
        return nn if nn @ (X[f[0]] - ctr) > 0 else -nn
    def face_at(z):
        for (A, B, Cc), fi in tri2:
            s = [((q - p).conjugate() * (z - p)).imag for p, q in ((A, B), (B, Cc), (Cc, A))]
            if min(s) > -1e-12 or max(s) < 1e-12: return fi
        raise ValueError(z)
    def th_solid(hg):
        a, b = complex(*hg[0]), complex(*hg[1]); m = (a + b) / 2; nv = (b - a) / abs(b - a) * 1j * 1e-5
        f1, f2 = face_at(m + nv), face_at(m - nv)
        return 0.0 if f1 == f2 else math.acos(max(-1.0, min(1.0, float(normal(f1) @ normal(f2)))))
    # 回転の中心の格子（中心どうしの差で張る）の基底: 短いものから 2 本、3 本の差がどれも整数で書けることを確かめる
    vs = [cen[k] - cen[0] for k in (1, 2, 3)]
    cand = sorted({(round(w.real, 9), round(w.imag, 9)) for k in itertools.product(range(-3, 4), repeat=3)
                   for w in [sum(x * y for x, y in zip(k, vs))] if abs(w) > 1e-9}, key=lambda w: math.hypot(*w))
    e1 = complex(*cand[0]); e2 = next(complex(*w) for w in cand if abs((complex(*w) * e1.conjugate()).imag) > 1e-9)
    for v in vs:
        M = np.linalg.solve([[e1.real, e2.real], [e1.imag, e2.imag]], [v.real, v.imag]); assert np.abs(M - M.round()).max() < 1e-9, M
    net = {'poly': [[p.real, p.imag] for p in pts], 'M': 6, 'centers': [[q.real, q.imag] for q in cen],
           'e1': [e1.real, e1.imag], 'e2': [e2.real, e2.imag]}
    MR.make(net, out, title, '切頂四面体', th_solid, 8, xname=tetra_name)
    d = json.load(open(out)); d['targets'][1]['label'] = tetra_name; json.dump(d, open(out, 'w'))
    # 確かめ: 折った形の凸包が、切頂四面体（辺 1、体積 23√2/12）と等面四面体（体積は 4 つの中心から）に合う。表面積はどちらも 7√3
    from scipy.spatial import ConvexHull
    p = [np.array([q.real, q.imag]) for q in cen]
    a2, b2, c2 = sorted([float(np.sum((p[i] - p[j]) ** 2)) for i, j in ((0, 1), (0, 2), (0, 3))])   # 4-base の平行四辺形の 3 辺（辺² の候補）
    vols = [23 * math.sqrt(2) / 12, None]
    for tg, V0 in zip(d['targets'], vols):
        Y = np.array([v for poly in tg['target3d'] for v in poly]); H = ConvexHull(Y)
        Qs = [np.array(q) for q in tg['target3d']]
        area = sum(np.linalg.norm(np.cross(q[k - 1] - q[0], q[k] - q[0])) / 2 for q in Qs for k in range(2, len(q)))
        print(tg['label'], '凸包の体積', round(H.volume, 9), '凸包の表面積', round(H.area, 9), '紙の面積', round(area, 9), '7√3 =', round(7 * math.sqrt(3), 9),
              '' if V0 is None else f'（期待 {V0:.9f}）')
        assert abs(H.area - 7 * math.sqrt(3)) < 1e-6 and abs(area - 7 * math.sqrt(3)) < 1e-6
        if V0 is not None: assert abs(H.volume - V0) < 1e-6
    return d


if __name__ == '__main__':
    d = build(13, os.path.join(HERE, 'data', 'data_refold_tt_regular.json'), '切頂四面体の展開図 → 正四面体（辺 √7）', '正四面体（辺 √7）')
    Y = np.array([v for poly in d['targets'][1]['target3d'] for v in poly]); from scipy.spatial import ConvexHull
    H = ConvexHull(Y); E = sorted({round(float(np.linalg.norm(Y[i] - Y[j])), 6) for i, j in itertools.combinations(H.vertices, 2)})
    print('正四面体の頂点', len(H.vertices), '頂点どうしの距離', E, '√7 =', round(math.sqrt(7), 6), '体積の期待 7√7/(6√2) =', round(7 * math.sqrt(7) / (6 * math.sqrt(2)), 9))
