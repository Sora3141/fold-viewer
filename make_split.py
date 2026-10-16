#!/usr/bin/env python3
"""分割展開図のデータを作る → data/data_split_*.json
入力: 分割展開図/work/n/*.json の hits[i]（cells = 直方体 A の展開図，parts = {箱: セル}）。
1 枚の展開図 A を格子線で切ると、部品がそれぞれ別の箱の展開図になる。データは 2 つの目標を持つ:
 0「A に折る」… 全体を A に折る（fold_export の通常の出力）
 1「分けて折る」… 部品ごとに自分の箱へ同時に折る。ヒンジ木が部品ごとの森なので target に parent/hinge/order を持ち、
    各部品の根パネルは roots[パネル] = [平らなときの平行移動, 折り切ったときの平行移動]（パネルは 2 つの目標で共通 = マス）。
 target.cut = パネルごとの部品番号（部品の境目を濃い線で描く）。色は部品の色相 × 面ごとの明るさ。
usage: python3 make_split.py   （研究リポジトリの中で実行）"""
import sys, os, json, colorsys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, '開いた箱', 'src')); os.chdir(os.path.join(ROOT, '開いた箱'))
import fold_export as fe, obox

def D(X): X = np.array([v for poly in X for v in poly]); return np.linalg.norm(X[:, None] - X[None], axis=2)
def dims(P3):
    X = np.array([v for poly in P3 for v in poly]); return X.max(0) - X.min(0), (X.max(0) + X.min(0)) / 2
def exp(spec, cells, tmp='/tmp/_split_part.json'):
    S = obox.surface(spec); cells = sorted(cells)
    st = next(s for s in (obox.fold_from(S, cells, s0) for s0 in range(S['NQ'] * 4)) if s is not None)
    fe.export({'cells': cells, 'folds': [[st]]}, [spec], tmp, spec); return json.load(open(tmp))
def corner(p): return (int(p['poly'][0][0]), int(p['poly'][0][1]))   # 箱はマスを切らない: パネルの最初の頂点 = マスの左下

def verify(d):
    """各目標を全ヒンジ目標角で折り、探索の面配置（target3d）と相似、箱の寸法が合うことを確かめる。2 つ目の目標は部品ごとに折る"""
    for ti, t in enumerate(d['targets']):
        t = {**{k: d[k] for k in ('parent', 'hinge', 'order')}, **t}   # 目標が持たなければ全体の木
        for part in sorted(set(t['cut'])) if ti else [None]:
            ids = [i for i in t['order'] if part is None or t['cut'][i] == part]; loc = {g: l for l, g in enumerate(ids)}
            pa = {loc[g]: (loc[t['parent'][g]] if t['parent'][g] in loc else None) for g in ids}
            P3 = fe.fold3d([d['panels'][g] for g in ids], pa, {loc[g]: t['hinge'][g] for g in ids}, [loc[g] for g in ids], {loc[g]: t['theta'][g] for g in ids})
            D1, D2 = D(P3), D([t['target3d'][g] for g in ids]); k = (D1 * D2).sum() / (D2 * D2).sum()
            assert np.abs(D1 - k * D2).max() < 1e-6, (t['label'], part)
            name = d['A'] if part is None else d['parts'][part]
            assert sorted(round(float(x), 6) for x in dims(P3)[0]) == sorted(float(x) for x in name.split('x')), (name, dims(P3)[0])
    assert all(t['cut'][i] == d['targets'][0]['cut'][i] for t in d['targets'] for i in range(len(t['cut'])))

def make(src, hit, out, title):
    rec = json.load(open(src)); h = rec['hits'][hit]; A = rec['A']
    parts = sorted(h['parts'], key=lambda s: sorted(map(int, s.split('x'))))   # 小さい箱から
    cells = [tuple(c) for c in h['cells']]; dA = exp(A, cells); panels = dA['panels']; n = len(cells)
    assert all(len(p['poly']) == 4 for p in panels)
    gid = {corner(p): i for i, p in enumerate(panels)}; assert len(gid) == n
    part_of = [None] * n
    for k, s in enumerate(parts):
        for c in h['parts'][s]: part_of[gid[tuple(c)]] = k
    parent = [-1] * n; hinge = [[[0, 0], [0, 0]] for _ in range(n)]; theta = [0.0] * n; order = []; tar3 = [None] * n; face = [0] * n; info = []
    for k, s in enumerate(parts):
        dp = exp(s, [tuple(c) for c in h['parts'][s]]); P = dp['panels']; t0 = dp['targets'][0]
        g = [gid[corner(p)] for p in P]   # 部品内の番号 → 全体の番号
        for l, gl in enumerate(g):
            assert P[l]['poly'] == panels[gl]['poly']
            parent[gl] = g[dp['parent'][l]] if dp['parent'][l] >= 0 else -1; hinge[gl] = dp['hinge'][l]
            theta[gl] = t0['theta'][l]; tar3[gl] = t0['target3d'][l]; face[gl] = fe.COLORS.index(t0['colors'][l])
        order += [g[l] for l in dp['order']]
        P3 = fe.fold3d(P, {l: (dp['parent'][l] if dp['parent'][l] >= 0 else None) for l in range(len(P))}, dp['hinge'], dp['order'], {l: t0['theta'][l] for l in range(len(P))})
        info.append((g[0], g, *dims(P3)))
    # 平らなとき: 部品の重心を全体の重心から遠ざける（マスの間に隙間 0.12 以上できる最小の倍率）
    cen = lambda ids: np.mean([[panels[i]['poly'][0][0] + .5, panels[i]['poly'][0][1] + .5] for i in ids], axis=0)
    F = cen(range(n)); fc = [cen(i[1]) for i in info]
    def gap(kk):
        pos = [(np.array(corner(panels[i])) + kk * (fc[part_of[i]] - F), part_of[i]) for i in range(n)]
        return min((max(abs(a[0] - b[0]), abs(a[1] - b[1])) for a, p in pos for b, q in pos if p != q), default=9)
    kk = next((x / 20 for x in range(2, 100) if gap(x / 20) >= 1.12), 5.0); print('離す倍率', kk, '隙間', round(gap(kk), 3))
    # 折り切ったとき: 箱を横一列（A の長い辺の向き）に重ならないよう並べる
    ax = 0 if max(c[0] for c in cells) - min(c[0] for c in cells) >= max(c[1] for c in cells) - min(c[1] for c in cells) else 1
    wid = [max(i[2][0], i[2][1]) for i in info]; mid = np.cumsum([0] + wid[:-1]) + np.arange(len(wid)) * 0.6 + np.array(wid) / 2
    mid -= (mid[-1] + wid[-1] / 2) / 2; roots = {}
    for k, (r, g, sz, ctr) in enumerate(info):
        tgt = np.array([F[0], F[1], 0.0]); tgt[ax] += mid[k]; off1 = kk * (fc[k] - F)
        roots[r] = [[round(float(off1[0]), 4), round(float(off1[1]), 4), 0.0], [round(float(x), 4) for x in tgt - ctr]]
    # 色: 部品の色相 × 面ごとの明るさ（A に折るときも同じ。切れ目が分かる）
    col = lambda k, f: '#%02x%02x%02x' % tuple(int(255 * v) for v in colorsys.hls_to_rgb((k * 0.381966 + 0.05) % 1, 0.9 - 0.045 * f, 0.62))
    colors = [col(part_of[i], face[i]) for i in range(n)]
    tA = dA['targets'][0]; tA['label'] = 'A（%s）に折る' % A.replace('x', '×'); tA['colors'] = colors; tA['cut'] = part_of
    tS = {'label': '分けて折る', 'theta': theta, 'colors': colors, 'cut': part_of, 'target3d': tar3, 'parent': parent, 'hinge': hinge, 'order': order, 'roots': roots}
    d = {'title': title, 'cells': dA['cells'], 'panels': panels, 'parent': dA['parent'], 'hinge': dA['hinge'], 'order': dA['order'], 'targets': [tA, tS], 'A': A, 'parts': parts}
    verify(d); json.dump(d, open(os.path.join(HERE, 'data', 'data_split_' + out + '.json'), 'w'))
    print('OK', out, A, parts, 'パネル', n)

WORK = os.path.join(ROOT, '分割展開図', 'work', 'n')
if __name__ == '__main__':
    for f, out in [('n4z_1x1x11.json', 'n4_1x1x11'), ('n5_1x1x17_e1.json', 'n5_1x1x17'), ('n6_1x1x21_r13.json', 'n6_1x1x21'), ('n4z_1x3x5.json', 'n4_1x3x5'), ('n5_2x2x7_e1.json', 'n5_2x2x7')]:
        make(os.path.join(WORK, f), 0, out, out)
