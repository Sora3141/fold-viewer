#!/usr/bin/env python3
"""data/data_split_*.json を全ヒンジ目標角で折ったとき、全体は A、各部品はそれぞれの箱に閉じるか（探索の面配置と相似・箱の寸法）。make_split.verify を使う。
usage: python3 test_split.py   （研究リポジトリの中で実行）"""
import glob, json, os
import make_split as ms
for f in sorted(glob.glob(os.path.join(ms.HERE, 'data', 'data_split_*.json'))):
    ms.verify(json.load(open(f))); print('OK', os.path.basename(f))
