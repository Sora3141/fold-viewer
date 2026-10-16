#!/usr/bin/env python3
"""1 枚の紙が面数の違う複数の凸多面体に折れる展開図のデータを作る → data/data_multi_*.json
入力: 多面数/figures/ex1_creases.json（ポリオミノ）・ex2_creases.json・ex2_9_creases.json（ポリイアモンド）。多面数/src/creases_multi.py の出力
（多面体の座標・面・辺と、各辺が紙の上を通る線分）。折り線は格子と限らない測地線。パネル＝マス（三角形）を折り線で切った凸多角形。
折り角は「2 つの面の外向き法線のなす角」、辺でない境目は 0（make_dodeca.py と同じ向き）。
usage: python3 make_multi.py   （研究リポジトリの中で実行）"""
import sys, os, json, math, colorsys
from collections import deque
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, '多面数', 'src'))
from make_dodeca import fe, _split_chord          # 開いた箱/src を path に入れ、cwd を移す
import creases_multi as cm
split = _split_chord(fe.split_convex)

def palette(n):
    """面ごとの色（隣どうしが似ないよう黄金角で色相を散らし、明るさも変える）"""
    return ['#%02x%02x%02x' % tuple(int(255 * v) for v in colorsys.hls_to_rgb(((i * 0.381966) % 1), 0.84 - 0.08 * (i % 3), 0.7)) for i in range(n)]

def on_seg(m, a, b, eps=1e-7):
    cr = (b[0] - a[0]) * (m[1] - a[1]) - (b[1] - a[1]) * (m[0] - a[0])
    return abs(cr) <= eps and min(a[0], b[0]) - eps <= m[0] <= max(a[0], b[0]) + eps and min(a[1], b[1]) - eps <= m[1] <= max(a[1], b[1]) + eps

def make(src, out, title, labels):
    C = json.load(open(src)); kind = C['kind']; G = C['gluings']; keys = list(G)
    paper = cm.Paper(kind, C['cells'], G[keys[0]]['pairs'])
    segs = [(tuple(p0), tuple(p1)) for k in keys for ch in G[k]['chords'].values() for p0, p1 in ch]
    panels = []
    for V in paper.V:
        polys = [[(z.real, z.imag) for z in V]]
        for a, b in segs: polys = [q for pg in polys for q in split(pg, a, b)]
        for pg in polys:
            while abs((pg[1][0] - pg[0][0]) * (pg[2][1] - pg[0][1]) - (pg[1][1] - pg[0][1]) * (pg[2][0] - pg[0][0])) < 1e-9: pg = pg[1:] + pg[:1]   # ビューアは先頭 3 点で向きを取るので、一直線に並ばない 3 点から始める
            panels.append({'poly': [list(p) for p in pg]})
    # 紙の上の隣り合うパネル（辺が重なる）→ 幅優先の木（根はパネル 0）
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
    targets = []
    for k, lab in zip(keys, labels):
        g = G[k]; X = np.array(g['X']); faces, edges, normals = cm.polyhedron(X); assert faces == g['faces']
        L = []
        for a, b, f1, f2 in edges:
            th = math.acos(max(-1, min(1, float(normals[f1] @ normals[f2]))))
            L += [(p0, p1, th) for p0, p1 in map(lambda c: map(tuple, c), g['chords'][f'{a}-{b}'])]
        def angle(h0, h1, poly):
            m = ((h0[0] + h1[0]) / 2, (h0[1] + h1[1]) / 2)
            cx = sum(p[0] for p in poly) / len(poly); cy = sum(p[1] for p in poly) / len(poly)
            side = (h1[0] - h0[0]) * (cy - h0[1]) - (h1[1] - h0[1]) * (cx - h0[0])
            for a, b, th in L:
                if on_seg(h0, a, b) and on_seg(h1, a, b) and on_seg(m, a, b): return (-1 if side > 0 else 1) * th   # 紙の表を外側にして折る
            return 0.0
        th = {c: angle(*hinge[c], panels[c]['poly']) for c in order[1:]}; P3 = fe.fold3d(panels, parent, hinge, order, th)
        key = {}; cid = []
        for poly in P3:
            p0 = np.array(poly[0]); n = sum(np.cross(np.array(poly[i]) - p0, np.array(poly[(i + 1) % len(poly)]) - p0) for i in range(len(poly))); n /= np.linalg.norm(n)   # 折れ線が一直線の頂点があっても向きが取れる
            cid.append(key.setdefault((*np.round(n, 3) + 0.0, round(float(n @ p0), 3) + 0.0), len(key)))
        assert len(key) == len(faces), (k, len(key), len(faces)); pal = palette(len(key))
        targets.append({'label': lab.format(F=k, V=len(X), E=len(edges)), 'theta': [th.get(i, 0.0) for i in range(len(panels))], 'colors': [pal[c] for c in cid], 'target3d': P3})
    cells = [list(c) for c in C['cells']] if kind == 'sq' else [[sum(z.real for z in V) / len(V) - .5, sum(z.imag for z in V) / len(V) - .5] for V in paper.V]
    data = {'title': title, 'cells': cells, 'panels': panels, 'parent': [parent[i] if parent[i] is not None else -1 for i in range(len(panels))],
            'hinge': [hinge.get(i, [[0, 0], [0, 0]]) for i in range(len(panels))], 'order': order, 'targets': targets}
    json.dump(data, open(out, 'w')); print(out, len(panels), 'panels', [t['label'] for t in targets])

if __name__ == '__main__':
    F = os.path.join(ROOT, '多面数', 'figures'); D = os.path.join(HERE, 'data')
    make(os.path.join(F, 'ex1_creases.json'), os.path.join(D, 'data_multi_poly_a10.json'), '面積 10：4・6・8・12 面体',
         ['{F} 面体（頂点 {V}・辺 {E}）'] * 4)
    make(os.path.join(F, 'ex2_creases.json'), os.path.join(D, 'data_multi_ico_5solids.json'), '正二十面体の辺展開図：4・6・8・12・20 面体',
         ['{F} 面体（頂点 {V}・辺 {E}）'] * 5)
    make(os.path.join(F, 'ex2_9_creases.json'), os.path.join(D, 'data_multi_ico_9solids.json'), '正二十面体の辺展開図：4・6・7・8・10・11・12・16・20 面体',
         ['{F} 面体（頂点 {V}・辺 {E}）'] * 9)
