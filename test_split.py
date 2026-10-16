#!/usr/bin/env python3
"""data/data_split_*.json を全ヒンジ目標角で折ったとき、全体は A、各部品はそれぞれの箱に閉じるか（探索の面配置と相似・箱の寸法）。make_split.verify を使う。
usage: python3 test_split.py   （研究リポジトリの中で実行）"""
import glob, json, os
import numpy as np
import make_split as ms
for f in sorted(glob.glob(os.path.join(ms.HERE, 'data', 'data_split_*.json'))):
    d = json.load(open(f)); ms.verify(d); bx = ms.part_boxes(d); r = d['targets'][1]['roots']
    xy = [b['ctr'][:2] + np.array(r[str(b['root'])][1][:2]) for b in bx]
    for i, a in enumerate(bx):
        assert np.abs(xy[i] - a['flat']).max() < 3, (f, i, xy[i], a['flat'])   # 箱の中心は部品の平らなときの重心から 3 マス以内
        for j in range(i + 1, len(bx)): assert ms.overlap(a, xy[i], bx[j], xy[j]) is None, (f, i, j)   # 折り切った箱は重ならない
    c = np.array(d['cells'], float) + .5; th = d['spin'] - ms.PHI0; q = c @ np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]]).T   # 回したあとの長い方の向きを x にそろえた座標
    assert np.ptp(q[:, 0]) >= np.ptp(q[:, 1]) - 1.5, (f, np.ptp(q, axis=0))   # 長い方が初期の視点で横（x と y の幅の差が 1.5 マス以内でも可 = ほぼ正方形）
    assert (ms.spin_of(d['cells']) == d['spin'])
    print('OK', os.path.basename(f))
