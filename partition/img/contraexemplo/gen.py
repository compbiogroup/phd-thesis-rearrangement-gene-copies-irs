#!/usr/bin/env python3
"""Gera as figuras dos contraexemplos do Algoritmo 1 (CompatibilityWithText.md).

Uso: python3 gen.py   (escreve os .svg nesta pasta)

Cada cromossomo é (circular?, genes, IRs). As adjacências são numeradas
a partir de 1, como (c,i) no texto. Marcações por figura:
  cls[(g,c,i)]  -> classificação da adjacência em T^min(G,H): 'bc' | 'perf'
  cuts[(g,c,i)] -> tipo de corte: 'p1' | 'p2' | 'alt' | 'opt'
  notes         -> [(g, c, i_ini, i_fim, texto)] colchete sobre genes i_ini..i_fim
  graph         -> componente perfeita: {'G': [(c, i, adj)], 'H': [...], 'edges': [(iG, iH)]};
                   vértices cortados na figura são riscados e o emparelhamento é recalculado
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

GW, GH = 46, 30        # gene
AW = 62                # espaço da adjacência
ROWH = 104             # altura de uma linha de cromossomo
LEFT = 70

COL = {
    'bc': '#2563eb', 'perf': '#16a34a', 'none': '#9ca3af',
    'p1': '#dc2626', 'p2': '#ea580c', 'alt': '#9333ea', 'opt': '#111827',
}
CLS_TXT = {'bc': 'adjacência em algum B_C', 'perf': 'componente perfeita (fora de T^min(G,H))'}
CUT_TXT = {'p1': 'corte da fase 1 (em B_C)', 'p2': 'corte da fase 2 (fora de B_C)',
           'alt': 'corte da saída alternativa', 'opt': 'breakpoint do ótimo'}


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def gene(x, y, g):
    """Desenha um gene em (x,y) = canto superior esquerdo."""
    out = []
    if g.lstrip('+-−').startswith('#'):
        lab = g.lstrip('+-−')
        out.append(f'<rect x="{x+8}" y="{y}" width="{GW-16}" height="{GH}" rx="3" '
                   f'fill="#e5e7eb" stroke="#4b5563"/>')
        out.append(f'<text x="{x+GW/2}" y="{y+GH/2+5}" text-anchor="middle" '
                   f'font-size="13">{esc(lab)}</text>')
        return out
    neg = g[0] in '-−'
    lab = g.lstrip('+-−')
    t = 10
    if neg:
        pts = [(x, y+GH/2), (x+t, y), (x+GW, y), (x+GW, y+GH), (x+t, y+GH)]
        fill = '#fde68a'
    else:
        pts = [(x, y), (x+GW-t, y), (x+GW, y+GH/2), (x+GW-t, y+GH), (x, y+GH)]
        fill = '#bfdbfe'
    p = ' '.join(f'{a:.1f},{b:.1f}' for a, b in pts)
    out.append(f'<polygon points="{p}" fill="{fill}" stroke="#374151"/>')
    sign = '−' if neg else '+'
    out.append(f'<text x="{x+GW/2}" y="{y+GH/2+5}" text-anchor="middle" '
               f'font-size="14" font-weight="bold">{sign}{esc(lab)}</text>')
    return out


def cut_mark(x, y, kind):
    c = COL[kind]
    dash = ' stroke-dasharray="5,3"' if kind == 'opt' else ''
    return [f'<line x1="{x}" y1="{y+2}" x2="{x}" y2="{y+GH-2}" stroke="{c}" '
            f'stroke-width="4"{dash}/>']


def chrom(x, y, gi, ci, ch, cls, cuts, used):
    circ, genes, irs = ch
    out = []
    n = len(genes)
    nadj = n if circ else n - 1
    gx = [x + k * (GW + AW) for k in range(n)]
    for i in range(1, nadj + 1):
        key = (gi, ci, i)
        k = cls.get(key, 'none')
        if k != 'none':
            used['cls'].add(k)
        c = COL[k]
        w = 3 if k != 'none' else 1.5
        x1 = gx[i-1] + GW
        x2 = x1 + AW
        ym = y + GH / 2
        out.append(f'<line x1="{x1}" y1="{ym}" x2="{x2}" y2="{ym}" stroke="{c}" stroke-width="{w}"/>')
        if circ and i == n:   # fecho do circular: volta por baixo
            out.append(f'<path d="M {x2} {ym} C {x2+20} {ym}, {x2+20} {y+GH+22}, {x2} {y+GH+22} '
                       f'L {gx[0]-12} {y+GH+22} C {gx[0]-30} {y+GH+22}, {gx[0]-30} {ym}, {gx[0]} {ym}" '
                       f'fill="none" stroke="{c}" stroke-width="{w}"/>')
        mx = x1 + AW / 2
        if key in cuts:
            used['cut'].add(cuts[key])
            out += cut_mark(mx, y, cuts[key])
        out.append(f'<text x="{mx}" y="{y-5}" text-anchor="middle" font-size="12" '
                   f'fill="#111827">{esc(irs[i-1])}</text>')
        out.append(f'<text x="{mx}" y="{y+GH+12}" text-anchor="middle" font-size="9" '
                   f'fill="#6b7280">({ci},{i})</text>')
    for k in range(n):
        out += gene(gx[k], y, genes[k])
    width = n * GW + nadj * AW + (25 if circ else 0)
    return out, gx, width


def max_matching(edges, nG):
    """Emparelhamento máximo (caminhos aumentantes); edges = [(iG, iH)]."""
    adj = {g: [h for a, h in edges if a == g] for g in range(nG)}
    mh = {}

    def aug(g, seen):
        for h in adj[g]:
            if h not in seen:
                seen.add(h)
                if h not in mh or aug(mh[h], seen):
                    mh[h] = g
                    return True
        return False
    for g in range(nG):
        aug(g, set())
    return {(g, h) for h, g in mh.items()}


def perfect_graph(x0, y0, h, graph, cuts):
    """Painel com o grafo de compatibilidade da componente perfeita de T^min(G,H),
    sem as adjacências já cortadas na figura."""
    c = COL['perf']
    gv, hv = graph['G'], graph['H']
    n = max(len(gv), len(hv))
    sp = 112
    w = max(n * sp, 230) + 20
    alive = {'G': [('G', ci, i) not in cuts for ci, i, _ in gv],
             'H': [('H', ci, i) not in cuts for ci, i, _ in hv]}
    live = [(a, b) for a, b in graph['edges'] if alive['G'][a] and alive['H'][b]]
    mat = max_matching(live, len(gv))
    mG = {a for a, _ in mat}
    mH = {b for _, b in mat}
    perfect = len(mat) == sum(alive['G']) == sum(alive['H'])
    out = [f'<rect x="{x0}" y="{y0-34}" width="{w}" height="{h+34}" rx="8" fill="#f0fdf4" '
           f'stroke="{c}" stroke-dasharray="4,3"/>',
           f'<text x="{x0+w/2}" y="{y0-14}" text-anchor="middle" font-size="12" font-weight="bold" '
           f'fill="{c}">componente perfeita (passo 1)</text>']
    yg, yh = y0 + 40, y0 + h - 68

    def xs(k, m):
        return x0 + w / 2 + (k - (m - 1) / 2) * sp
    pg = [xs(k, len(gv)) for k in range(len(gv))]
    ph = [xs(k, len(hv)) for k in range(len(hv))]
    for a, b in graph['edges']:
        if (a, b) in live:
            sw, op = (3.5, 1) if (a, b) in mat else (1.5, 0.6)
            out.append(f'<line x1="{pg[a]}" y1="{yg}" x2="{ph[b]}" y2="{yh}" stroke="{c}" '
                       f'stroke-width="{sw}" stroke-opacity="{op}"/>')
        else:
            out.append(f'<line x1="{pg[a]}" y1="{yg}" x2="{ph[b]}" y2="{yh}" stroke="#d1d5db" '
                       f'stroke-width="1.2" stroke-dasharray="3,3"/>')
    for side, vs, ps, yy, dy1, dy2, ms in (('G', gv, pg, yg, -27, -13, mG), ('H', hv, ph, yh, 24, 38, mH)):
        for k, ((ci, i, adj), px) in enumerate(zip(vs, ps)):
            nm = f'{side}({ci},{i})'
            if not alive[side][k]:
                cc = COL[cuts[(side, ci, i)]]
                out.append(f'<circle cx="{px}" cy="{yy}" r="7" fill="white" stroke="#9ca3af" '
                           f'stroke-width="1.5" stroke-dasharray="2,2"/>')
                out.append(f'<path d="M {px-8} {yy-8} L {px+8} {yy+8} M {px-8} {yy+8} L {px+8} {yy-8}" '
                           f'stroke="{cc}" stroke-width="3"/>')
                tc = '#9ca3af'
            else:
                sc = c if (k in ms or perfect) else COL['p1']
                fill = 'white' if sc == c else '#fee2e2'
                out.append(f'<circle cx="{px}" cy="{yy}" r="7" fill="{fill}" stroke="{sc}" stroke-width="2.5"/>')
                tc = '#111827'
            out.append(f'<text x="{px}" y="{yy+dy1}" text-anchor="middle" font-size="11" '
                       f'font-weight="bold" fill="{tc}">{esc(nm)}</text>')
            out.append(f'<text x="{px}" y="{yy+dy2}" text-anchor="middle" '
                       f'font-size="11" fill="{tc if tc != "#111827" else "#374151"}">{esc(adj)}</text>')
    if not any(alive['G']) and not any(alive['H']):
        msg, mc = 'todas as adjacências foram cortadas', '#6b7280'
    elif perfect:
        msg, mc = 'arestas grossas: um emparelhamento perfeito', '#6b7280'
    else:
        msg, mc = 'sem emparelhamento perfeito (vértices em vermelho sem par)', COL['p1']
    out.append(f'<text x="{x0+w/2}" y="{y0+h-8}" text-anchor="middle" font-size="10" fill="{mc}">'
               f'{msg}</text>')
    return out, w


def figure(name, title, G, H, cls=None, cuts=None, notes=None, caption=None, graph=None):
    cls, cuts, notes = cls or {}, cuts or {}, notes or []
    used = {'cls': set(), 'cut': set()}
    body = []
    y = 84
    pos = {}
    maxw = 0
    for gi, gen in (('G', G), ('H', H)):
        body.append(f'<text x="18" y="{y+GH/2+6}" font-size="18" font-weight="bold">{gi}</text>')
        for ci, ch in enumerate(gen, 1):
            o, gx, w = chrom(LEFT, y, gi, ci, ch, cls, cuts, used)
            body += o
            pos[(gi, ci)] = (gx, y)
            maxw = max(maxw, LEFT + w + 40)
            y += ROWH
        y += 14
    for gi, ci, a, b, txt, *kind in notes:
        gx, yy = pos[(gi, ci)]
        yb = yy - 24
        c = '#be185d'
        # bloco que passa pelo fecho de um circular: um colchete aberto em cada ponta
        spans = [(a, b, 'MLLL')] if a <= b else [(a, len(gx), 'MLL-'), (1, b, '-LLL')]
        for sa, sb, _ in spans:
            x1, x2 = gx[sa-1] - 3, gx[sb-1] + GW + 3
            if a <= b:
                d = f'M {x1} {yb+6} L {x1} {yb} L {x2} {yb} L {x2} {yb+6}'
            elif sa == a:
                d = f'M {x1} {yb+6} L {x1} {yb} L {x2+20} {yb}'
            else:
                d = f'M {x1-20} {yb} L {x2} {yb} L {x2} {yb+6}'
            body.append(f'<path d="{d}" fill="none" stroke="{c}" stroke-width="1.5"/>')
            body.append(f'<text x="{(x1+x2)/2}" y="{yb-4}" text-anchor="middle" font-size="12" '
                        f'fill="{c}" font-weight="bold">{esc(txt)}</text>')
    if graph:
        o, gw = perfect_graph(maxw + 10, 84, y - 14 - 84, graph, cuts)
        body += o
        maxw += gw + 10
    # legenda
    ly = y + 4
    items = [('line', k, CLS_TXT[k]) for k in ('bc', 'perf') if k in used['cls']]
    items += [('cut', k, CUT_TXT[k]) for k in ('p1', 'p2', 'alt', 'opt') if k in used['cut']]
    lx = LEFT
    for typ, k, txt in items:
        wtxt = 7 * len(txt) + 40
        if lx + wtxt > max(maxw, 640):
            lx = LEFT
            ly += 22
        if typ == 'line':
            body.append(f'<line x1="{lx}" y1="{ly}" x2="{lx+24}" y2="{ly}" stroke="{COL[k]}" stroke-width="3"/>')
        else:
            dash = ' stroke-dasharray="4,2"' if k == 'opt' else ''
            body.append(f'<line x1="{lx+12}" y1="{ly-9}" x2="{lx+12}" y2="{ly+9}" stroke="{COL[k]}" '
                        f'stroke-width="3.5"{dash}/>')
        body.append(f'<text x="{lx+30}" y="{ly+4}" font-size="12">{esc(txt)}</text>')
        lx += wtxt
    H_ = ly + 18
    if caption:
        for line in caption.split('\n'):
            H_ += 4
            body.append(f'<text x="18" y="{H_+10}" font-size="12" fill="#374151">{esc(line)}</text>')
            H_ += 14
        H_ += 8
    W = max(maxw, 660)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H_}" '
           f'viewBox="0 0 {W} {H_}" font-family="DejaVu Sans, Arial, sans-serif">',
           f'<rect width="{W}" height="{H_}" fill="white"/>',
           f'<text x="18" y="26" font-size="16" font-weight="bold">{esc(title)}</text>']
    svg += body + ['</svg>']
    with open(os.path.join(HERE, name + '.svg'), 'w') as f:
        f.write('\n'.join(svg))


def tag(g, lst, kind):
    return {(g, c, i): kind for c, i in lst}


def merge(*ds):
    r = {}
    for d in ds:
        r.update(d)
    return r


# ---------------------------------------------------------------- circular
CG = [(True, ['+1', '+1'], ['2', '1']),
      (True, ['−1', '+1', '−1'], ['0', '0', '2'])]
CH = [(True, ['−1', '+1', '−1'], ['[1:2]', '[2:2]', '[1:1]']),
      (True, ['−1', '−1'], ['[2:2]', '[1:2]'])]
c_cls = merge(tag('G', [(2, 1), (2, 2)], 'bc'), tag('H', [(1, 1), (1, 2)], 'bc'),
              tag('G', [(1, 1), (1, 2), (2, 3)], 'perf'), tag('H', [(1, 3), (2, 1), (2, 2)], 'perf'))
c_p1 = merge(tag('G', [(2, 1), (2, 2)], 'p1'), tag('H', [(1, 1), (1, 2)], 'p1'))
c_p2 = merge(tag('G', [(1, 2)], 'p2'), tag('H', [(2, 2)], 'p2'))

c_graph = {'G': [(1, 1, '+1 2 +1'), (1, 2, '+1 1 +1'), (2, 3, '−1 2 −1')],
           'H': [(1, 3, '−1 [1:1] −1'), (2, 1, '−1 [2:2] −1'), (2, 2, '−1 [1:2] −1')],
           'edges': [(0, 1), (0, 2), (1, 0), (1, 2), (2, 1), (2, 2)]}
figure('circ_1_tmin', 'Contraexemplo circular — T^min(G,H) e conjuntos B_C', CG, CH, c_cls, graph=c_graph,
       caption='As 4 junções +1/−1 formam componentes de um vértice (caso 1); os pares +1 +1 e −1 −1\n'
               'formam uma componente perfeita e não estão em nenhum B_C.')
figure('circ_2_fase1', 'Contraexemplo circular — fase 1', CG, CH, c_cls, c_p1, graph=c_graph,
       notes=[('G', 2, 3, 1, 'sem par'), ('H', 1, 3, 1, 'sem par')],
       caption='Os cortes são os breakpoints do ótimo nesses cromossomos e T^min fica vazio, mas\n'
               'G[−1 ─2─ −1] e H[−1 ─[1:1]─ −1] não têm correspondente como blocos.')
figure('circ_3_fase2', 'Contraexemplo circular — fase 2 (antes da escolha de D)', CG, CH, c_cls, graph=c_graph,
       cuts=merge(c_p1, c_p2),
       notes=[('G', 2, 3, 1, 'B'), ('H', 2, 1, 2, 'D'), ('G', 1, 1, 2, 'parceiro de D')],
       caption='D é cortado em H(2,2) e o parceiro em G(1,2): os dois cortes caem na componente perfeita,\n'
               'fora de qualquer B_C (contradiz partition.tex:270). Depois a saída alternativa corta\n'
               'H(1,3) e G(1,1), custo final 8.')
figure('circ_4_otimo', 'Contraexemplo circular — ótimo (custo 6)', CG, CH, c_cls, graph=c_graph,
       cuts=merge(tag('G', [(2, 1), (2, 2), (1, 1)], 'opt'), tag('H', [(1, 1), (1, 2), (2, 1)], 'opt')),
       caption='O ótimo corta os dois circulares na outra IR: G(1,1) e H(2,1).\n'
               '[+1 ─1─ +1] ↔ [−1 ─[1:1]─ −1] (invertido), [−1 ─2─ −1] ↔ [−1 ─[1:2]─ −1], [+1] ↔ [+1].')

# ---------------------------------------------------------------- linear
LG = [(False, ['#', '+1', '−1', '−1', '+1', '#'], ['2', '2', '2', '2', '1']),
      (False, ['#', '−1', '+1', '−1', '+1', '#'], ['2', '2', '2', '0', '1'])]
LH = [(False, ['#', '+1', '−1', '−1', '+1', '#'], ['[0:1]', '[2:2]', '[0:0]', '[1:2]', '[1:2]']),
      (False, ['#', '−1', '+1', '−1', '+1', '#'], ['[0:1]', '[2:2]', '[1:1]', '[1:1]', '[2:3]'])]
l_bcG = [(1, 1), (1, 3), (2, 4), (1, 2), (2, 3)]
l_bcH = [(1, 2), (1, 1), (1, 3), (2, 3), (2, 4)]
l_cls = merge(tag('G', l_bcG, 'bc'), tag('H', l_bcH, 'bc'),
              tag('G', [(1, 4), (2, 2)], 'perf'), tag('H', [(1, 4), (2, 2)], 'perf'))
l_p1 = merge(tag('G', l_bcG, 'p1'), tag('H', l_bcH, 'p1'))
l_p2a = merge(tag('G', [(1, 4)], 'p2'), tag('H', [(2, 2)], 'p2'))
l_p2b = merge(tag('G', [(2, 2)], 'p2'), tag('H', [(1, 4)], 'p2'))

l_graph = {'G': [(1, 4, '−1 2 +1'), (2, 2, '−1 2 +1')],
           'H': [(1, 4, '−1 [1:2] +1'), (2, 2, '−1 [2:2] +1')],
           'edges': [(0, 0), (0, 1), (1, 0), (1, 1)]}
figure('lin_1_tmin', 'Contraexemplo linear — T^min(G,H) e conjuntos B_C', LG, LH, l_cls, graph=l_graph,
       caption='8 componentes do caso 1 (azul). As adjacências −1 +1 em verde formam uma componente\n'
               'perfeita (2 × 2) e não estão em nenhum B_C.')
figure('lin_2_fase1', 'Contraexemplo linear — fase 1 (10 cortes em B_C)', LG, LH, l_cls, l_p1, graph=l_graph,
       notes=[('G', 2, 5, 6, 'sem par'), ('H', 2, 5, 6, 'sem par')],
       caption='T^min fica vazio; sobram sem par G[+1 ─1─ #] e H[+1 ─[2:3]─ #] (1 ∉ [2:3]).')
figure('lin_3_fase2a', 'Contraexemplo linear — fase 2, primeira quebra', LG, LH, l_cls, graph=l_graph,
       cuts=merge(l_p1, l_p2a),
       notes=[('G', 2, 5, 6, 'B'), ('H', 2, 1, 3, 'D'), ('G', 1, 4, 6, 'parceiro de D')],
       caption='Com a escolha de D: corte em H(2,2) e no parceiro em G(1,4), ambos na componente perfeita.')
figure('lin_4_fase2b', 'Contraexemplo linear — fase 2, segunda quebra (custo 14)', LG, LH, l_cls, graph=l_graph,
       cuts=merge(l_p1, l_p2a, l_p2b),
       notes=[('G', 2, 5, 6, 'B'), ('H', 1, 4, 6, 'D'), ('G', 2, 1, 3, 'parceiro de D')],
       caption='Corte em H(1,4) e G(2,2). Os quatro cortes da fase 2 estão fora de qualquer B_C\n'
               '(contradiz partition.tex:270), embora a saída alternativa não seja usada.')
figure('lin_5_otimo', 'Contraexemplo linear — ótimo (custo 10)', LG, LH, l_cls, graph=l_graph,
       cuts=merge(tag('G', [(1, 1), (1, 3), (2, 2), (2, 3), (2, 4)], 'opt'),
             tag('H', [(1, 1), (1, 3), (2, 2), (2, 3), (2, 4)], 'opt')),
       caption='O ótimo mantém G(1,2) e H(1,2) (+1 ─2─ −1 ↔ +1 ─[2:2]─ −1) e corta só G(2,3) da\n'
               'componente {G(1,2), G(2,3), H(1,2)}; o Approx corta os três já na fase 1.')

# ---------------------------------------------------------------- telômeros com sinal
TG = [(False, ['#L', '+1', '+1', '+2', '−2', '−1', '#R'], ['2', '1', '0', '1', '1', '2'])]
TH = [(False, ['#L', '+1', '+2', '−2', '−1', '−1', '#R'],
       ['[0:0]', '[0:0]', '[0:0]', '[0:1]', '[0:1]', '[0:1]'])]
t_bcG, t_bcH = [(1, 1), (1, 4), (1, 6)], [(1, 1), (1, 3), (1, 6)]
t_cls = merge(tag('G', t_bcG, 'bc'), tag('H', t_bcH, 'bc'),
              tag('G', [(1, 2)], 'perf'), tag('H', [(1, 5)], 'perf'))
figure('tel_1_fase2', 'Exemplo com telômeros com sinal — fases 1 e 2 (custo 8 = ótimo)', TG, TH, t_cls,
       merge(tag('G', t_bcG, 'p1'), tag('H', t_bcH, 'p1'),
             tag('G', [(1, 2)], 'p2'), tag('H', [(1, 5)], 'p2')),
       notes=[('G', 1, 5, 6, 'B'), ('H', 1, 4, 6, 'D'), ('G', 1, 2, 4, 'parceiro de D')],
       graph={'G': [(1, 2, '+1 1 +1')], 'H': [(1, 5, '−1 [0:1] −1')], 'edges': [(0, 0)]},
       caption='Fase 1 corta os B_C (telômeros e +2 −2). Na fase 2, D é cortado em H(1,5) e o parceiro\n'
               'em G(1,2): adjacências de uma componente perfeita, fora de qualquer B_C.')
