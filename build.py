#!/usr/bin/env python3
"""Build every app listed in apps.json.
  <file>.html        artifact version (body only; Claude Artifacts wraps it)
  docs/<file>.html   standalone page for GitHub Pages
  docs/index.html    gallery
Add an app: put its data json in data/ (made by fold_export.py in the research repo), add an entry to apps.json, run build.py."""
import json, html
apps=json.load(open('apps.json')); tpl=open('template.html').read()
HEAD='<!doctype html>\n<html lang="ja">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width,initial-scale=1">\n'
for a in apps:
    d=json.load(open(a['data']))
    body=tpl.replace('__TITLE__',a['title']).replace('__HEADING__',a['heading']).replace('__SUB__',a['sub']).replace('/*DATA*/null',json.dumps(d,ensure_ascii=False))
    open(a['file']+'.html','w').write(body)
    # standalone: move <title>/<link>/<style> into head
    i=body.index('<div id="stage">'); head_part=body[:i]; rest=body[i:]
    open(f"docs/{a['file']}.html",'w').write(HEAD+head_part+'</head>\n<body>\n'+rest+'\n</body>\n</html>\n')
    print('built',a['file'])
cards="\n".join(f'''<a class="card" href="{a['file']}.html"><h2>{html.escape(a['title'])}</h2><p class="sub">{html.escape(a['sub'])}</p><p>{html.escape(a['note'])}</p></a>''' for a in apps)
index=f'''<!doctype html><html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>展開図フォールディング</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Zen+Kaku+Gothic+New:wght@500;700&display=swap">
<style>
:root{{--bg:#eef0f3;--ink:#1b2027;--muted:#5f6772;--line:#c9ced6;--card:#fff;--accent:#c8332b}}
@media (prefers-color-scheme:dark){{:root{{--bg:#15181d;--ink:#e8eaee;--muted:#9aa3ad;--line:#2c323a;--card:#1c2026;--accent:#e0554a}}}}
body{{margin:0;background:var(--bg);color:var(--ink);font-family:"Zen Kaku Gothic New","Hiragino Sans","Noto Sans JP",sans-serif}}
main{{max-width:880px;margin:0 auto;padding:48px 20px}}
h1{{font-size:26px;margin:0 0 6px}} .lead{{color:var(--muted);margin:0 0 32px;max-width:60ch;line-height:1.7}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:18px}}
.card{{display:block;background:var(--card);border:1px solid var(--line);border-radius:12px;padding:18px 20px;color:inherit;text-decoration:none;line-height:1.6}}
.card:hover{{border-color:var(--accent)}} .card h2{{font-size:18px;margin:0 0 4px}} .card .sub{{color:var(--muted);font-size:13px;margin:0 0 10px}} .card p{{margin:0;font-size:14px}}
footer{{color:var(--muted);font-size:12px;margin-top:40px;line-height:1.7}}
</style></head><body><main>
<h1>展開図フォールディング</h1>
<p class="lead">一枚のポリオミノが複数の立体に折れる「共通展開図」を，折り目の角度を一斉に動かして 3D で見せるページ集。JAIST 上原研究室での展開図研究（同じ箱の多重折り，箱・開いた箱・ポリキューブの共通展開図）から，見せたい例を順に追加していく。</p>
<div class="grid">
{cards}
</div>
<footer>操作：ドラッグで回転，ホイールで拡大，下のスライダーで折り進み。折り切った形が探索で得た面配置と一致することを数値的に確認したデータを使っている。途中の形は全折り目を同時に回しているだけなので，面がすれ違うことがある。</footer>
</main></body></html>'''
open('docs/index.html','w').write(index); print('built docs/index.html')
