# fold-viewer — 展開図フォールディング

一枚のポリオミノが複数の立体（箱・開いた箱・ポリキューブ・√k 立方体）に折れる様子を 3D で見せる小さな Web アプリ集。
GitHub Pages: https://sora3141.github.io/fold-viewer/

## 構成
- `template.html` … アプリのひな形（three.js，折り目の角度を 0→目標値で一斉に動かす）
- `data/*.json` … 展開図ごとのデータ（パネル多角形，ヒンジ木，立体ごとの折り角と色）。研究用リポジトリの `openbox/src/fold_export.py` が探索結果から生成する
- `apps.json` … 公開するアプリの一覧（題名・説明・データ）
- `build.py` … `apps.json` から `docs/<name>.html`（公開ページ）と `<name>.html`（Claude Artifact 用の本文のみ版）と `docs/index.html`（一覧）を生成
- `docs/` … GitHub Pages で配信されるファイル

## アプリを追加する
1. 探索結果の 1 行（cells と各立体への折り方）から `fold_export.py` でデータ json を作り `data/` に置く
2. `apps.json` に項目を追加
3. `python3 build.py` を実行してコミット

## データの検証
`fold_export.py` は各立体について，ヒンジ木の角度をすべて目標値にしたときの 3D 位置が探索の面配置と一致することを（剛体変換を除いて）確認している。
