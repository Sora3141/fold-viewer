#!/usr/bin/env python3
"""Build every app listed in apps.json.
  <file>.html        artifact version (body only; Claude Artifacts wraps it)
  docs/<file>.html   standalone page for GitHub Pages
  docs/index.html    gallery
Add an app: put its data json in data/ (made by fold_export.py in the research repo), add an entry to apps.json, run build.py."""
import json, html

import math
def _proj(pt,yaw=0.65,pitch=0.5):
    x,y,z=pt
    # rotate about vertical axis (z up in data? data uses z as the box height axis) -> view: yaw about z, then pitch
    cx,sx=math.cos(yaw),math.sin(yaw); x1=cx*x-sx*y; y1=sx*x+cx*y; z1=z
    cp,sp=math.cos(pitch),math.sin(pitch); y2=cp*y1-sp*z1; z2=sp*y1+cp*z1
    return (x1, -z2, y2)   # screen x, screen y (down), depth (larger = farther)
def svg_solid(d,ti,size=110):
    t=d['targets'][ti]; polys=[]
    for pid,poly3 in enumerate(t['target3d']):
        pts=[_proj(tuple(v)) for v in poly3]; dep=sum(p[2] for p in pts)/len(pts)
        polys.append((dep,[(p[0],p[1]) for p in pts],t['colors'][pid]))
    polys.sort(key=lambda x:-x[0])   # far first
    xs=[x for _,pg,_ in polys for x,y in pg]; ys=[y for _,pg,_ in polys for x,y in pg]
    mnx,mxx,mny,mxy=min(xs),max(xs),min(ys),max(ys); sc=(size-10)/max(mxx-mnx,mxy-mny,1e-9)
    ox=(size-(mxx-mnx)*sc)/2-mnx*sc; oy=(size-(mxy-mny)*sc)/2-mny*sc
    out=[f'<svg viewBox="0 0 {size} {size}" width="{size}" height="{size}" role="img" aria-label="{t["label"]}">']
    for _,pg,col in polys:
        out.append('<polygon points="'+' '.join(f'{x*sc+ox:.1f},{y*sc+oy:.1f}' for x,y in pg)+f'" fill="{col}" stroke="#2a2f36" stroke-width="0.8" stroke-linejoin="round"/>')
    out.append('</svg>'); return ''.join(out)
def svg_net(d,size=110):
    t=d['targets'][0]; polys=[(pan['poly'],t['colors'][i]) for i,pan in enumerate(d['panels'])]
    xs=[x for pg,_ in polys for x,y in pg]; ys=[y for pg,_ in polys for x,y in pg]
    mnx,mxx,mny,mxy=min(xs),max(xs),min(ys),max(ys); sc=(size-8)/max(mxx-mnx,mxy-mny)
    ox=(size-(mxx-mnx)*sc)/2-mnx*sc; oy=(size-(mxy-mny)*sc)/2+mxy*sc
    out=[f'<svg viewBox="0 0 {size} {size}" width="{size}" height="{size}" role="img" aria-label="展開図">']
    for pg,col in polys:
        out.append('<polygon points="'+' '.join(f'{x*sc+ox:.1f},{oy-y*sc:.1f}' for x,y in pg)+f'" fill="{col}" stroke="#2a2f36" stroke-width="0.6"/>')
    out.append('</svg>'); return ''.join(out)
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
def card(a):
    d=json.load(open(a['data'])); n=len(d['targets']); show=list(range(min(n,5)))
    thumbs='<div class="thumbs">'+f'<figure>{svg_net(d,96)}<figcaption>展開図</figcaption></figure>'+''.join(f'<figure>{svg_solid(d,i,96)}<figcaption>{html.escape(d["targets"][i]["label"])}</figcaption></figure>' for i in show)+(f'<figure class="more">+{n-5}</figure>' if n>5 else '')+'</div>'
    return f'''<a class="card" href="{a['file']}.html"><h2>{html.escape(a['title'])}</h2><p class="sub">{html.escape(a['sub'])}</p>{thumbs}<p>{html.escape(a['note'])}</p></a>'''
INTRO={'同じ箱の多重折り':'一枚の展開図が同じ箱に本質的に異なる複数の方法で折れる例。折り線の入れ替わりに注目。',
 '異なる箱の共通展開図':'一枚の展開図が形の違う複数の箱に折れる例。多重折りをもつものを選んだ。',
 'ふたの無い箱（開いた箱）':'直方体の 1 面を除いたゴミ箱型の箱。開いた箱どうし，開いた箱と閉じた箱の共通展開図と多重折り。',
 'ポリキューブ':'立方体をつなげた立体（トリ・テトラ・ペンタキューブ，十字形）と箱の共通展開図。'}
sections=[]
for c in dict.fromkeys(a['category'] for a in apps):
    grp=[a for a in apps if a['category']==c]
    sections.append(f'<section><h2 class="cat">{html.escape(c)}<span class="count">{len(grp)} 本</span></h2><p class="intro">{html.escape(INTRO.get(c,""))}</p><div class="grid">'+"\n".join(card(a) for a in grp)+'</div></section>')
cards="\n".join(sections)
index=f'''<!doctype html><html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>展開図フォールディング</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Zen+Kaku+Gothic+New:wght@500;700&display=swap">
<style>
:root{{--bg:#eef0f3;--ink:#1b2027;--muted:#5f6772;--line:#c9ced6;--card:#fff;--accent:#c8332b}}
@media (prefers-color-scheme:dark){{:root{{--bg:#15181d;--ink:#e8eaee;--muted:#9aa3ad;--line:#2c323a;--card:#1c2026;--accent:#e0554a}}}}
body{{margin:0;background:var(--bg);color:var(--ink);font-family:"Zen Kaku Gothic New","Hiragino Sans","Noto Sans JP",sans-serif}}
main{{max-width:1180px;margin:0 auto;padding:48px 20px}}
h1{{font-size:26px;margin:0 0 6px}} .lead{{color:var(--muted);margin:0 0 32px;max-width:60ch;line-height:1.7}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(360px,1fr));gap:18px}}
.cat{{font-size:20px;margin:36px 0 4px;padding-bottom:6px;border-bottom:2px solid var(--line)}} .cat .count{{font-size:13px;color:var(--muted);font-weight:500;margin-left:10px}} .intro{{color:var(--muted);margin:0 0 14px;font-size:14px}}
.thumbs{{display:flex;flex-wrap:wrap;gap:6px;margin:10px 0 12px}} .thumbs figure{{margin:0;text-align:center;width:96px}} .thumbs figcaption{{font-size:10px;color:var(--muted);line-height:1.3;margin-top:2px;word-break:keep-all}} .thumbs svg{{display:block;background:var(--bg);border-radius:6px}} .thumbs .more{{display:flex;align-items:center;justify-content:center;height:96px;color:var(--muted);font-size:18px}}
.card{{display:block;background:var(--card);border:1px solid var(--line);border-radius:12px;padding:18px 20px;color:inherit;text-decoration:none;line-height:1.6}}
.card:hover{{border-color:var(--accent)}} .card h2{{font-size:18px;margin:0 0 4px}} .card .sub{{color:var(--muted);font-size:13px;margin:0 0 10px}} .card p{{margin:0;font-size:14px}}
footer{{color:var(--muted);font-size:12px;margin-top:40px;line-height:1.7}}
</style></head><body><main>
<h1>展開図フォールディング</h1>
<p class="lead">一枚のポリオミノが複数の立体に折れる「共通展開図」を，折り目の角度を一斉に動かして 3D で見せるページ集。JAIST 上原研究室での展開図研究（同じ箱の多重折り，箱・開いた箱・ポリキューブの共通展開図）から，見せたい例を順に追加していく。</p>
{cards}
<footer>操作：ドラッグで回転，ホイールで拡大，下のスライダーで折り進み。折り切った形が探索で得た面配置と一致することを数値的に確認したデータを使っている。途中の形は全折り目を同時に回しているだけなので，面がすれ違うことがある。</footer>
</main></body></html>'''
open('docs/index.html','w').write(index); print('built docs/index.html')
