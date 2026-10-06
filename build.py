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
def _hull(P):
    P=sorted(set(P)); 
    if len(P)<3: return P
    cr=lambda o,a,b:(a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0]); lo=[]; up=[]
    for q in P:
        while len(lo)>=2 and cr(lo[-2],lo[-1],q)<=1e-12: lo.pop()
        lo.append(q)
    for q in reversed(P):
        while len(up)>=2 and cr(up[-2],up[-1],q)<=1e-12: up.pop()
        up.append(q)
    return lo[:-1]+up[:-1]
def _area(pg): return abs(sum(pg[i][0]*pg[(i+1)%len(pg)][1]-pg[(i+1)%len(pg)][0]*pg[i][1] for i in range(len(pg))))/2
def _plane(poly3):
    a,b,c=poly3[:3]; u=[b[i]-a[i] for i in range(3)]; v=[c[i]-a[i] for i in range(3)]
    n=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]; L=math.sqrt(sum(x*x for x in n)) or 1; n=[x/L for x in n]
    return tuple(round(x,3) for x in n)+(round(sum(n[i]*a[i] for i in range(3)),3),)
def svg_solid(d,ti,size=110):
    """立体の絵。面ごとに塗り，パネルの継ぎ目（折り線や他の立体の折り目）は描かず，面の縁だけを描く。"""
    t=d['targets'][ti]; faces={}
    for pid,poly3 in enumerate(t['target3d']):
        faces.setdefault((_plane(poly3),t['colors'][pid]),[]).append([_proj(tuple(v)) for v in poly3])
    polys=[]
    for (pl,col),pans in faces.items():
        dep=sum(p[2] for pg in pans for p in pg)/sum(len(pg) for pg in pans)
        flat=[[(p[0],p[1]) for p in pg] for pg in pans]; H=_hull([q for pg in flat for q in pg])
        convex=len(H)>=3 and abs(_area(H)-sum(_area(pg) for pg in flat))<1e-6*max(1,_area(H))
        polys.append((dep,flat,H if convex else None,col,(pl,col)))
    polys.sort(key=lambda x:-x[0])   # far first
    # 凸な立体（箱・四面体）は裏を向いた面を描かない。細長い四面体では面の重心の深さで並べると裏の面が手前に来るため
    P=[p for pans in faces.values() for pg in pans for p in pg]; cen=[sum(p[k] for p in P)/len(P) for k in range(3)]
    def facing(pans):
        a,b,c=pans[0][:3]; u=[b[k]-a[k] for k in range(3)]; v=[c[k]-a[k] for k in range(3)]
        n=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]; off=sum(n[k]*a[k] for k in range(3))
        side=[sum(n[k]*q[k] for k in range(3))-off for q in P]; L=math.sqrt(sum(x*x for x in n)) or 1
        if sum(n[k]*(a[k]-cen[k]) for k in range(3))<0: n=[-x for x in n]; side=[-x for x in side]
        return max(side)<=1e-6*L*max(1,max(abs(x) for x in P[0])), n[2]<0   # (立体がこの面の内側にある, 視点（深さの小さい側）を向く)
    F=[facing(pans) for pans in faces.values()]
    if all(cv for cv,_ in F):
        vis={k for k,(cv,fr) in zip(faces,F) if fr}; polys=[pp for pp in polys if pp[4] in vis]
    xs=[x for _,fl,_,_,_ in polys for pg in fl for x,y in pg]; ys=[y for _,fl,_,_,_ in polys for pg in fl for x,y in pg]
    mnx,mxx,mny,mxy=min(xs),max(xs),min(ys),max(ys); sc=(size-10)/max(mxx-mnx,mxy-mny,1e-9)
    ox=(size-(mxx-mnx)*sc)/2-mnx*sc; oy=(size-(mxy-mny)*sc)/2-mny*sc
    pts=lambda pg:' '.join(f'{x*sc+ox:.1f},{y*sc+oy:.1f}' for x,y in pg)
    out=[f'<svg viewBox="0 0 {size} {size}" width="{size}" height="{size}" role="img" aria-label="{t["label"]}">']
    for _,fl,H,col,_ in polys:
        if H: out.append(f'<polygon points="{pts(H)}" fill="{col}" stroke="#2a2f36" stroke-width="0.8" stroke-linejoin="round"/>')
        else:   # 凸でない面（ポリキューブ）は継ぎ目を塗りの色で隠す
            out+= [f'<polygon points="{pts(pg)}" fill="{col}" stroke="{col}" stroke-width="0.6" stroke-linejoin="round"/>' for pg in fl]
    out.append('</svg>'); return ''.join(out)
def svg_net(d,size=110):
    """展開図の絵。立体 0 の面の色で塗り，折り線は描かず，紙の縁だけを描く。"""
    t=d['targets'][0]; polys=[(pan['poly'],t['colors'][i]) for i,pan in enumerate(d['panels'])]
    xs=[x for pg,_ in polys for x,y in pg]; ys=[y for pg,_ in polys for x,y in pg]
    mnx,mxx,mny,mxy=min(xs),max(xs),min(ys),max(ys); sc=(size-8)/max(mxx-mnx,mxy-mny)
    ox=(size-(mxx-mnx)*sc)/2-mnx*sc; oy=(size-(mxy-mny)*sc)/2+mxy*sc
    out=[f'<svg viewBox="0 0 {size} {size}" width="{size}" height="{size}" role="img" aria-label="展開図">']
    for pg,col in polys:
        out.append('<polygon points="'+' '.join(f'{x*sc+ox:.1f},{oy-y*sc:.1f}' for x,y in pg)+f'" fill="{col}" stroke="{col}" stroke-width="0.6"/>')
    E=[(pg[k],pg[(k+1)%len(pg)]) for pg,_ in polys for k in range(len(pg))]; G={}
    for i,(a,b) in enumerate(E):   # 辺を通るマスに登録（中点の近くだけ調べる）
        for gx in range(math.floor(min(a[0],b[0])-1e-9),math.floor(max(a[0],b[0])+1e-9)+1):
            for gy in range(math.floor(min(a[1],b[1])-1e-9),math.floor(max(a[1],b[1])+1e-9)+1): G.setdefault((gx,gy),[]).append(i)
    def on(m,a,b):
        cr=(b[0]-a[0])*(m[1]-a[1])-(b[1]-a[1])*(m[0]-a[0])
        return abs(cr)<1e-7 and min(a[0],b[0])-1e-9<=m[0]<=max(a[0],b[0])+1e-9 and min(a[1],b[1])-1e-9<=m[1]<=max(a[1],b[1])+1e-9
    seg=[]
    for i,(a,b) in enumerate(E):   # 紙の縁 = 中点がほかのパネルの辺に乗らない辺
        m=((a[0]+b[0])/2,(a[1]+b[1])/2)
        if not any(j!=i and on(m,*E[j]) for j in G.get((math.floor(m[0]),math.floor(m[1])),())):
            seg.append(f'M{a[0]*sc+ox:.1f} {oy-a[1]*sc:.1f}L{b[0]*sc+ox:.1f} {oy-b[1]*sc:.1f}')
    out.append(f'<path d="{"".join(seg)}" stroke="#2a2f36" stroke-width="0.8" fill="none" stroke-linecap="round"/>')
    out.append('</svg>'); return ''.join(out)
apps=json.load(open('apps.json')); tpl=open('template.html').read()
HEAD='<!doctype html>\n<html lang="ja">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">\n'
for a in apps:
    d=json.load(open(a['data']))
    sub=a['sub']+(('　／　出典：'+a['cite']) if a.get('cite') else '')
    body=tpl.replace('__TITLE__',a['title']).replace('__HEADING__',a['heading']).replace('__SUB__',sub).replace('__NOTE__',html.escape(a['note'])).replace('/*DATA*/null',json.dumps(d,ensure_ascii=False))
    open(a['file']+'.html','w').write(body.replace('__HOME__',''))   # artifact: no gallery to go back to
    body=body.replace('__HOME__','<a class="home" href="./" aria-label="一覧へ戻る" title="一覧へ戻る">←</a>')
    # standalone: move <title>/<link>/<style> into head
    i=body.index('<div id="stage">'); head_part=body[:i]; rest=body[i:]
    open(f"docs/{a['file']}.html",'w').write(HEAD+head_part+'</head>\n<body>\n'+rest+'\n</body>\n</html>\n')
    print('built',a['file'])
def card(a):
    d=json.load(open(a['data'])); n=len(d['targets']); show=list(range(min(n,5)))
    thumbs='<div class="thumbs">'+f'<figure>{svg_net(d,96)}<figcaption>展開図</figcaption></figure>'+''.join(f'<figure>{svg_solid(d,i,96)}<figcaption>{html.escape(d["targets"][i]["label"])}</figcaption></figure>' for i in show)+(f'<figure class="more">+{n-5}</figure>' if n>5 else '')+'</div>'
    area=len(d['cells'])
    return f'''<a class="card" href="{a['file']}.html"><h2><span class="area">面積 {area}</span>{html.escape(a['title'])}</h2><p class="sub">{html.escape(a['sub'])}</p>{thumbs}<p>{html.escape(a['note'])}</p>{('<p class="cite">出典：'+html.escape(a['cite'])+'</p>') if a.get('cite') else ''}</a>'''
INTRO={'先行研究の例':'論文に掲載されている共通展開図・多重折りの例を，図から読み取って折り方を検証したもの。出典を各カードに示す。',
 '同じ箱の多重折り':'一枚の展開図が同じ箱に本質的に異なる複数の方法で折れる例。折り線の入れ替わりに注目。',
 '異なる箱の共通展開図':'一枚の展開図が形の違う複数の箱に折れる例。多重折りをもつものを選んだ。',
 'ふたの無い箱（開いた箱）':'直方体の 1 面を除いたゴミ箱型の箱。開いた箱どうし，開いた箱と閉じた箱の共通展開図と多重折り。',
 '箱と四面体':'一枚の展開図が 2 種類の箱と，4 面が合同な鋭角三角形の四面体（等面四面体）に折れる例。四面体の折り線は斜めの格子の直線で，二面角から折り角を出している。',
 'ポリキューブ':'立方体をつなげた立体（トリ・テトラ・ペンタキューブ，十字形）と箱の共通展開図。'}
sections=[]
for c in dict.fromkeys(a['category'] for a in apps):
    grp=[a for a in apps if a['category']==c]
    sections.append(f'<section id="c{len(sections)}"><h2 class="cat">{html.escape(c)}<span class="count">{len(grp)} 本</span></h2><p class="intro">{html.escape(INTRO.get(c,""))}</p><div class="grid">'+"\n".join(card(a) for a in grp)+'</div></section>')
cards="\n".join(sections)
cats=list(dict.fromkeys(a['category'] for a in apps))
nav='<nav class="toc" aria-label="分類">'+''.join(f'<a href="#c{i}">{html.escape(c)}<span>{sum(a["category"]==c for a in apps)}</span></a>' for i,c in enumerate(cats))+'</nav>'
index=f'''<!doctype html><html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>展開図フォールディング</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Zen+Kaku+Gothic+New:wght@500;700&display=swap">
<style>
:root{{--bg:#eef0f3;--ink:#1b2027;--muted:#5f6772;--line:#c9ced6;--card:#fff;--accent:#c8332b}}
@media (prefers-color-scheme:dark){{:root{{--bg:#15181d;--ink:#e8eaee;--muted:#9aa3ad;--line:#2c323a;--card:#1c2026;--accent:#e0554a}}}}
body{{margin:0;background:var(--bg);color:var(--ink);font-family:"Zen Kaku Gothic New","Hiragino Sans","Noto Sans JP",sans-serif}}
main{{max-width:1180px;margin:0 auto;padding:48px 16px}}
h1{{font-size:26px;margin:0 0 6px}} h1 .total{{font-size:13px;color:var(--muted);font-weight:500;margin-left:12px}} .lead{{color:var(--muted);margin:0 0 32px;max-width:60ch;line-height:1.7}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(min(360px,100%),1fr));gap:18px}}
.toc{{position:sticky;top:0;z-index:1;display:flex;gap:6px;overflow-x:auto;padding:10px 0;margin:0 0 8px;background:var(--bg);scrollbar-width:none}} .toc::-webkit-scrollbar{{display:none}}
.toc a{{flex:0 0 auto;font-size:13px;color:var(--ink);text-decoration:none;border:1px solid var(--line);background:var(--card);border-radius:999px;padding:5px 12px}} .toc a:hover{{border-color:var(--accent)}} .toc span{{color:var(--muted);font-size:11px;margin-left:6px}}
section{{scroll-margin-top:56px}}
.cat{{font-size:20px;margin:36px 0 4px;padding-bottom:6px;border-bottom:2px solid var(--line)}} .cat .count{{font-size:13px;color:var(--muted);font-weight:500;margin-left:10px}} .intro{{color:var(--muted);margin:0 0 14px;font-size:14px}}
.thumbs{{display:flex;flex-wrap:wrap;gap:6px;margin:10px 0 12px}} .thumbs figure{{margin:0;text-align:center;width:96px}} .thumbs figcaption{{font-size:10px;color:var(--muted);line-height:1.3;margin-top:2px;word-break:keep-all}} .thumbs svg{{display:block;background:var(--bg);border-radius:6px}} .thumbs .more{{display:flex;align-items:center;justify-content:center;height:96px;color:var(--muted);font-size:18px}}
.card{{display:block;background:var(--card);border:1px solid var(--line);border-radius:12px;padding:18px 20px;color:inherit;text-decoration:none;line-height:1.6}}
.card{{transition:border-color .15s,transform .15s,box-shadow .15s}} .card:hover{{border-color:var(--accent);transform:translateY(-2px);box-shadow:0 6px 20px #0001}} .card:focus-visible{{outline:2px solid var(--accent);outline-offset:2px}} @media (prefers-reduced-motion:reduce){{.card{{transition:none}} .card:hover{{transform:none}}}} .card h2{{font-size:18px;margin:0 0 4px;display:flex;align-items:center;gap:8px;flex-wrap:wrap}} .card .area{{font-size:11px;font-weight:600;color:var(--accent);border:1px solid var(--accent);border-radius:999px;padding:1px 8px;letter-spacing:.04em}} .card .sub{{color:var(--muted);font-size:13px;margin:0 0 10px}} .card .cite{{color:var(--muted);font-size:11px;margin-top:8px;border-top:1px dashed var(--line);padding-top:6px}} .card p{{margin:0;font-size:14px}}
footer{{color:var(--muted);font-size:12px;margin-top:40px;line-height:1.7}}
</style></head><body><main>
<h1>展開図フォールディング<span class="total">{len(apps)} 本</span></h1>
<p class="lead">一枚のポリオミノが複数の立体に折れる「共通展開図」を，折り目の角度を一斉に動かして 3D で見せるページ集。JAIST 上原研究室での展開図研究（同じ箱の多重折り，箱・開いた箱・ポリキューブの共通展開図）から，見せたい例を順に追加していく。</p>
{nav}
{cards}
<footer>操作：ドラッグで回転，ホイール・ピンチで拡大，ダブルクリックで視点を戻す，下のスライダーで折り進み，Space で再生，←→ で立体を切り替え。折り切った形が探索で得た面配置と一致することを数値的に確認したデータを使っている。途中の形は全折り目を同時に回しているだけなので，面がすれ違うことがある。</footer>
</main></body></html>'''
open('docs/index.html','w').write(index); print('built docs/index.html')
