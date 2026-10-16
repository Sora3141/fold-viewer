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

PHI0 = float(np.arctan(np.tan(0.55) * np.sin(0.6)))   # 初期の視点（template の PITCH0・YAW0）で画面の水平に見える、展開図の面内の向き
def spin_of(cells):
    """展開図の主軸（長い方）が初期の視点で画面の横に見えるように回す角（template の group.rotation.z）"""
    c = np.array(cells, float) + .5; best = None
    for deg in range(180):   # 外接する長方形の面積が最小になる向き（長い方を x に）
        a = np.radians(deg); q = c @ np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]]); w, h = np.ptp(q, axis=0)
        if best is None or w * h < best[0] - 1e-9: best = (w * h, a if w >= h else a + np.pi / 2)
    return round(float((PHI0 - best[1] + np.pi / 2) % np.pi - np.pi / 2), 6)
MARGIN = 0.25   # 折り切った箱どうしの最小の隙間
def part_boxes(d):
    """部品ごとの折り切った箱（根を動かさないとき）: root, size, ctr（箱の中心）, flat（平らなときの部品の重心 xy）。根パネルは平面 z=0 に残る"""
    t = d['targets'][1]; out = []
    for k in range(len(d['parts'])):
        ids = [i for i in t['order'] if t['cut'][i] == k]; loc = {g: l for l, g in enumerate(ids)}
        pa = {loc[g]: (loc[t['parent'][g]] if t['parent'][g] in loc else None) for g in ids}
        sz, ctr = dims(fe.fold3d([d['panels'][g] for g in ids], pa, {loc[g]: t['hinge'][g] for g in ids}, [loc[g] for g in ids], {loc[g]: t['theta'][g] for g in ids}))
        flat = np.mean([[d['panels'][i]['poly'][0][0] + .5, d['panels'][i]['poly'][0][1] + .5] for i in ids], axis=0)
        out.append({'root': next(g for g in ids if t['parent'][g] < 0), 'size': sz, 'ctr': ctr, 'flat': flat})
    return out
def overlap(a, pa, b, pb):
    """2 つの箱（xy の中心を pa, pb に置く）が MARGIN 以内に近ければ各軸のめり込み量、離れていれば None"""
    ov = [(a['size'][i] + b['size'][i]) / 2 + MARGIN - abs((pa[i] if i < 2 else a['ctr'][2]) - (pb[i] if i < 2 else b['ctr'][2])) for i in range(3)]
    return ov if min(ov) > 0 else None

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
    # 平らなとき: 部品はその場（隙間 0。境目は濃い線）。折り切ったとき: 箱は部品の重心の真上に置き、重なるときだけ最小限ずらす
    d0 = {'panels': panels, 'parts': parts, 'targets': [None, {'parent': parent, 'hinge': hinge, 'order': order, 'theta': theta, 'cut': part_of}]}
    bx = part_boxes(d0); xy = [b['flat'].copy() for b in bx]
    for _ in range(500):   # ponytail: 重なる組を浅い軸の方向へ半分ずつ押す素朴な反復。収束しなければ assert で気づく
        moved = False
        for i in range(len(bx)):
            for j in range(i + 1, len(bx)):
                ov = overlap(bx[i], xy[i], bx[j], xy[j])
                if ov is not None:
                    a = int(np.argmin(ov[:2])); sg = 1 if xy[i][a] >= xy[j][a] else -1
                    xy[i][a] += sg * (ov[a] / 2 + 1e-3); xy[j][a] -= sg * (ov[a] / 2 + 1e-3); moved = True
        if not moved: break
    assert not moved
    roots = {b['root']: [[0.0, 0.0, 0.0], [round(float(xy[k][0] - b['ctr'][0]), 4), round(float(xy[k][1] - b['ctr'][1]), 4), 0.0]] for k, b in enumerate(bx)}
    # 色: 部品の色相 × 面ごとの明るさ（A に折るときも同じ。切れ目が分かる）
    col = lambda k, f: '#%02x%02x%02x' % tuple(int(255 * v) for v in colorsys.hls_to_rgb((k * 0.381966 + 0.05) % 1, 0.9 - 0.045 * f, 0.62))
    colors = [col(part_of[i], face[i]) for i in range(n)]
    tA = dA['targets'][0]; tA['label'] = 'A（%s）に折る' % A.replace('x', '×'); tA['colors'] = colors; tA['cut'] = part_of
    tS = {'label': '分けて折る', 'theta': theta, 'colors': colors, 'cut': part_of, 'target3d': tar3, 'parent': parent, 'hinge': hinge, 'order': order, 'roots': roots}
    d = {'title': title, 'cells': dA['cells'], 'panels': panels, 'parent': dA['parent'], 'hinge': dA['hinge'], 'order': dA['order'], 'targets': [tA, tS], 'A': A, 'parts': parts, 'spin': spin_of(cells)}
    verify(d); json.dump(d, open(os.path.join(HERE, 'data', 'data_split_' + out + '.json'), 'w'))
    print('OK', out, A, parts, 'パネル', n)

WORK = os.path.join(ROOT, '分割展開図', 'work', 'n')
if __name__ == '__main__':
    for f, out in [('n4z_1x1x11.json', 'n4_1x1x11'), ('n5_1x1x17_e1.json', 'n5_1x1x17'), ('n6_1x1x21_r13.json', 'n6_1x1x21'), ('n4z_1x3x5.json', 'n4_1x3x5'), ('n5_2x2x7_e1.json', 'n5_2x2x7'), ('hk7/n_3x4x6_s20.json', 'n7_3x4x6'), ('hk7/n_1x1x31_s1.json', 'n7_1x1x31')]:
        make(os.path.join(WORK, f), 0, out, out)
