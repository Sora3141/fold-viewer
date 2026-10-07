#!/usr/bin/env python3
"""面積 10 のポリオミノ（1×1×2 の箱にも面 12 の凸多面体にも折れる）のデータを作る → data/data_dodeca_a10.json
入力: 十二面体/figures/area10_ex2_creases.json（十二面体/src/creases.py の出力: 多面体の座標・面と、18 本の辺が紙の上を通る線分）
折り線は格子と限らない測地線の線分。マスを線分で切り、折り角は「2 つの面の外向き法線のなす角」、辺でない境目は 0。
usage: python3 make_dodeca.py   （研究リポジトリの中で実行）"""
import sys, os, json, math
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, '開いた箱', 'src')); os.chdir(os.path.join(ROOT, '開いた箱'))
import fold_export as fe, obox

CREASES = os.path.join(ROOT, '十二面体', 'figures', 'area10_ex2_creases.json')
COLORS12 = ['#f9c9c9', '#fde2b3', '#fff5b0', '#cdebd0', '#c6dcf5', '#e2cdf3',
            '#f4a6a6', '#f5c16c', '#d9e27a', '#8fd1a0', '#8dbbe8', '#c19be0']   # fe.COLORS の 6 色 + 濃い 6 色

def _split_chord(orig):
    """fe.split_convex は直線全体で切る。ここでは線分 (a,b) の上にある弦だけで切る（別のマスを通る延長では切らない）"""
    def f(poly, a, b, eps=1e-9):
        u = (b[0] - a[0], b[1] - a[1]); n = math.hypot(*u)
        d = [(u[0] * (p[1] - a[1]) - u[1] * (p[0] - a[0])) / n for p in poly]; pts = []
        for k in range(len(poly)):
            p, q, dp, dq = poly[k], poly[(k + 1) % len(poly)], d[k], d[(k + 1) % len(poly)]
            if abs(dp) < eps: pts.append(p)
            elif dp * dq < 0: t = dp / (dp - dq); pts.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
        if len(pts) < 2: return [poly]
        s = [((x - a[0]) * u[0] + (y - a[1]) * u[1]) / n for x, y in pts]; m = (min(s) + max(s)) / 2
        return orig(poly, a, b, eps) if -eps <= m <= n + eps else [poly]
    return f

def normals(X, faces):
    """面（三角形）の外向き単位法線と、辺 → 隣り合う 2 面"""
    N = []
    for a, b, c in faces:
        n = np.cross(X[b] - X[a], X[c] - X[a]); N.append(n / np.linalg.norm(n))
    return N

def dodeca_target(g, tag):
    X = np.array(g['X']); faces = [tuple(f) for f in g['faces']]; N = normals(X, faces)
    fold = {}   # 辺 (u,v) → 折り角（法線のなす角 = π − 二面角）
    for (u, v) in map(tuple, g['edges']):
        fs = [i for i, f in enumerate(faces) if u in f and v in f]; assert len(fs) == 2
        fold[(u, v)] = math.acos(max(-1, min(1, float(N[fs[0]] @ N[fs[1]]))))
    chords = [(tuple(map(float, (Fraction(p0[0]), Fraction(p0[1])))), tuple(map(float, (Fraction(p1[0]), Fraction(p1[1])))), fold[tuple(map(int, k.split('-')))])
              for k, ch in g['chords'].items() for p0, p1 in ch]
    seen = set(); L = []
    for a, b, th in chords:
        if (a, b) not in seen and (b, a) not in seen: seen.add((a, b)); L.append((a, b, th))
    assert len(L) == len({tuple(sorted(c[:2])) for c in L})
    def angle(h0, h1, poly):
        m = ((h0[0] + h1[0]) / 2, (h0[1] + h1[1]) / 2)
        cx = sum(p[0] for p in poly) / len(poly); cy = sum(p[1] for p in poly) / len(poly)
        side = (h1[0] - h0[0]) * (cy - h0[1]) - (h1[1] - h0[1]) * (cx - h0[0])
        for a, b, th in L:
            if all(abs((b[0] - a[0]) * (q[1] - a[1]) - (b[1] - a[1]) * (q[0] - a[0])) < 1e-7 for q in (h0, h1, m)) \
               and min(a[0], b[0]) - 1e-7 <= m[0] <= max(a[0], b[0]) + 1e-7 and min(a[1], b[1]) - 1e-7 <= m[1] <= max(a[1], b[1]) + 1e-7:
                return (-1 if side > 0 else 1) * th     # 紙の表を外側にして折る（四面体と同じ向き）
        return 0.0
    def colors(P3):
        key = {}; cols = []
        for poly in P3:
            p0, p1, p2 = map(np.array, poly[:3]); n = np.cross(p1 - p0, p2 - p0); n /= np.linalg.norm(n)
            cols.append(key.setdefault((*np.round(n, 3), round(float(n @ p0), 3)), len(key)))
        assert len(key) == 12, len(key)      # 面が 12 の平面に収まる
        return [COLORS12[c] for c in cols]
    return {'label': '12 面体（面 12・頂点 8・辺 18）' + tag, 'lines': [(a, b) for a, b, _ in L], 'angle': angle, 'colors': colors}

def build(out):
    C = json.load(open(CREASES)); cells = [tuple(c) for c in C['cells']]
    K = obox.kfull(obox.surface('1x1x2'), cells); assert len(K) == 1
    fe.split_convex = _split_chord(fe.split_convex)
    n = len(C['gluings']); tg = [dodeca_target(g, f'（貼り合わせ {j + 1}）' if n > 1 else '') for j, g in enumerate(C['gluings'])]
    fe.export({'cells': cells, 'folds': [[K[0]]]}, ['1x1x2'], out, '面積 10：1×1×2 と 12 面体', labels=['1×1×2'], extra=tg)

if __name__ == '__main__':
    build(os.path.join(HERE, 'data', 'data_dodeca_a10.json'))
