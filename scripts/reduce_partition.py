"""Turn the partitions produced by MCFISP-exe (partition/) into reduced genomes.

usage: python3 reduce_partition.py RED_FILE COST_FILE OUT_1 OUT_2 ...

OUT_i are the output files of MCFISP-exe for the instances of one input file,
in the order of the instances (MCFISP-exe -o PREFIX writes PREFIX_0001, ...).

RED_FILE gets one reduced instance per input instance (4 lines, in the format of
the thesis data set: genomes without caps, n + 1 intergenic regions), which can
be given to dec/dist of cyc with the option -e. Each block becomes a single
gene (see adj_graph.tex, "Fixando Adjacências com uma Partição de Genomas"):
  - the blocks are paired by a maximum matching between compatible blocks, and
    each pair gets a new label, with the sign of the block of the target
    inverted when the two blocks are compatible in reverse;
  - a block without a pair (exclusive block, unbalanced genomes) gets a label
    that appears only in its genome, so cyc treats it as an indel;
  - the intergenic regions of the breakpoints are kept and the ones inside a
    block are dropped, except the one between a cap and the first (last) gene of
    the block that contains the cap, which becomes the first (last) region.
The labels of the pairs make the reduced genomes permutations (plus indels), so
random mappings are not needed (dist -k 1).

COST_FILE gets one line per instance:
  breakpoints_origin breakpoints_target cost blocks_origin blocks_target exclusive_origin exclusive_target time
where cost = breakpoints_origin + breakpoints_target is the cost of the
partition (thesis definition). An instance without a partition (error or time
limit of MCFISP-exe) gets the line "NA" and a reduced instance with 4 lines
"NA", so the order of the instances is kept (dec and dist of cyc write "NA" as
the result of such an instance).
"""
import sys


def parse_side(genes_line, irs_line, flex):
    """List of chromosomes of one genome. Each chromosome is (circular, blocks,
    breakpoint_irs), each block is (genes, internal_irs)."""
    def ir(tok):
        if flex:
            lo, hi = tok.split(':')
            return (int(lo), int(hi))
        return int(tok)
    chroms = []
    gene_chrs = [c for c in genes_line.split(';') if c.strip()]
    ir_chrs = irs_line.split(';')
    for k, gc in enumerate(gene_chrs):
        toks = gc.split()
        circular = toks[0] == 'C'
        blocks_g = ' '.join(toks[1:]).split('|')
        parts_i = ir_chrs[k].split('|')
        blocks = []
        bp_irs = []
        for b, g in enumerate(blocks_g):
            genes = [int(x) for x in g.split()]
            internal = [ir(x) for x in parts_i[2 * b].split()]
            blocks.append((genes, internal))
            if 2 * b + 1 < len(parts_i):
                t = parts_i[2 * b + 1].split()
                if t:
                    bp_irs.append(ir(t[0]))
        chroms.append((circular, blocks, bp_irs))
    return chroms


def fits(v, iv):
    return iv[0] <= v <= iv[1]


def orientation(bg, bh):
    """'d' if block bg (origin) is compatible with bh (target) directly, 'r' if
    in reverse, None otherwise (the direct one is preferred)."""
    (gg, gi), (hg, hi) = bg, bh
    if len(gg) != len(hg):
        return None
    if gg == hg and all(fits(v, iv) for v, iv in zip(gi, hi)):
        return 'd'
    if [-x for x in reversed(gg)] == hg and all(fits(v, iv) for v, iv in zip(reversed(gi), hi)):
        return 'r'
    return None


def max_matching(n_g, adj):
    match_h = {}

    def augment(u, seen):
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                if v not in match_h or augment(match_h[v], seen):
                    match_h[v] = u
                    return True
        return False
    sys.setrecursionlimit(10000)
    for u in range(n_g):
        augment(u, set())
    return {u: v for v, u in match_h.items()}


def reduce_side(chroms, labels, flex):
    """Reduced chromosomes (genes, irs) of one genome. labels[(c, b)] is the
    signed label of block b of chromosome c. Only linear chromosomes."""
    def fmt(r):
        return f'{r[0]}:{r[1]}' if flex else str(r)
    CAP = None
    out_g, out_i = [], []
    for c, (circular, blocks, bp_irs) in enumerate(chroms):
        if circular:
            raise NotImplementedError('circular chromosomes are not supported')
        # items of the chromosome (caps and reduced genes) and the region between
        # each pair of consecutive items
        items, irs = [], []
        for b, (bg, bi) in enumerate(blocks):
            if all(x == 0 for x in bg):  # a cap alone
                sub, sub_irs = [CAP], []
            else:
                sub, sub_irs = [labels[(c, b)]], []
                if bg[0] == 0:  # the cap is glued to the block: keep its region
                    sub, sub_irs = [CAP] + sub, [bi[0]] + sub_irs
                if bg[-1] == 0:
                    sub, sub_irs = sub + [CAP], sub_irs + [bi[-1]]
            if items:
                irs.append(bp_irs[b - 1])
            items += sub
            irs += sub_irs
        assert items[0] is CAP and items[-1] is CAP and len(irs) == len(items) - 1
        out_g.append('L ' + ' '.join(map(str, items[1:-1])))
        out_i.append(' '.join(fmt(r) for r in irs))
    return ' ; '.join(out_g) + ' ;', ' ; '.join(out_i) + ' ;'


def reduce_instance(lines):
    g = parse_side(lines[0], lines[1], False)
    h = parse_side(lines[2], lines[3], True)
    bg = [(c, b) for c, ch in enumerate(g) for b in range(len(ch[1]))]
    bh = [(c, b) for c, ch in enumerate(h) for b in range(len(ch[1]))]
    block = lambda side, cb: side[cb[0]][1][cb[1]]
    # cap blocks ([0]) are not elements of the reduced genome and are not paired
    real_g = [i for i, cb in enumerate(bg) if any(x != 0 for x in block(g, cb)[0])]
    real_h = [j for j, cb in enumerate(bh) if any(x != 0 for x in block(h, cb)[0])]
    ori = {}
    adj = {i: [] for i in range(len(bg))}
    for i in real_g:
        for j in real_h:
            o = orientation(block(g, bg[i]), block(h, bh[j]))
            if o:
                ori[(i, j)] = o
                adj[i].append(j)
    match = max_matching(len(bg), adj)
    labels_g, labels_h = {}, {}
    label = 0
    for i in real_g:
        if i in match:
            label += 1
            labels_g[bg[i]] = label
            labels_h[bh[match[i]]] = label if ori[(i, match[i])] == 'd' else -label
    matched_h = set(match.values())
    excl_g = [i for i in real_g if i not in match]
    excl_h = [j for j in real_h if j not in matched_h]
    for i in excl_g:
        label += 1
        labels_g[bg[i]] = label
    for j in excl_h:
        label += 1
        labels_h[bh[j]] = label
    rg = reduce_side(g, labels_g, False)
    rh = reduce_side(h, labels_h, True)
    bps_g = sum(len(ch[2]) for ch in g)
    bps_h = sum(len(ch[2]) for ch in h)
    stats = (bps_g, bps_h, bps_g + bps_h, len(real_g), len(real_h), len(excl_g), len(excl_h))
    return [rg[0], rg[1], rh[0], rh[1]], stats


def main():
    red_file, cost_file, outs = sys.argv[1], sys.argv[2], sys.argv[3:]
    with open(red_file, 'w') as red, open(cost_file, 'w') as cost:
        for out in outs:
            with open(out) as f:
                lines = [l.rstrip('\n') for l in f]
            time = next((l.split(':', 1)[1].strip() for l in lines if l.startswith('# Time')), 'NA')
            if len(lines) < 4 or lines[0].startswith('#'):
                red.write('NA\nNA\nNA\nNA\n')
                cost.write('NA\n')
                continue
            reduced, stats = reduce_instance(lines[:4])
            red.write('\n'.join(reduced) + '\n')
            cost.write(' '.join(map(str, stats)) + f' {time}\n')


if __name__ == '__main__':
    main()
