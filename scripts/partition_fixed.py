"""Write the partitions produced by MCFISP-exe (partition/) in positions of the
original genomes, without reducing them (for the adjacency graph with fixed
adjacencies, and optionally fixed characters, of the thesis).

usage: python3 partition_fixed.py FIX_FILE COST_FILE PART_FILE INPUT_FILE

PART_FILE has the outputs of MCFISP-exe for the instances of INPUT_FILE (the
input of MCFISP-exe), joined in the order of the instances (the PART files of
run_part_tese.sh). Each output ends with a line "# Time: ...". INPUT_FILE is
needed to place the blocks of circular chromosomes (MCFISP-exe writes first the
block that crosses the origin).

FIX_FILE gets 3 lines per instance, with the positions of the genomes of the
input file. Genes are numbered 1..n in each chromosome (caps excluded) and
intergenic region i is the one after gene i:
  - linear chromosome: regions 0..n, region 0 is the one before gene 1 (after
    the cap) and region n the one after gene n (before the cap);
  - circular chromosome: regions 1..n, region n is the one between gene n and
    gene 1 (the last region of the input).
Multiple chromosomes are separated by ';' as in the input. The lines are:
  1. the breakpoints of the origin: the indices of the intergenic regions whose
     adjacency is a breakpoint of the partition of the origin (all the other
     adjacencies are inside a block, i.e., fixed);
  2. the breakpoints of the target, in the same way;
  3. the correspondence of genes given by a bijection between the common blocks:
     for each gene j of the origin, the position of its pair in the target,
     negative when the two blocks are paired in reverse, and 0 when the gene is in
     an exclusive block (no pair). Genes are counted over all chromosomes.
A line without breakpoints is "-" (so no line is empty). An instance without a
partition (error or time limit of MCFISP-exe) gets the three lines "NA".

COST_FILE gets one line per instance, as in reduce_partition.py:
  breakpoints_origin breakpoints_target cost blocks_origin blocks_target exclusive_origin exclusive_target time
"""
import sys

from reduce_partition import fits, max_matching, orientation, parse_side


def read_input_side(genes_line, irs_line, flex):
    """Chromosomes of one genome of the input file: (circular, genes, irs), with
    the intergenic regions as text in the format of the output of MCFISP-exe (a
    rigid region v of the target is written v:v)."""
    chrs = [c.split() for c in genes_line.split(';') if c.strip()]
    irs = [c.split() for c in irs_line.split(';')]
    out = []
    for k, c in enumerate(chrs):
        circular = c[0] == 'C'
        genes = [int(x) for x in (c[1:] if c[0] in 'LC' else c)]
        regions = [r if not flex or ':' in r else f'{r}:{r}' for r in irs[k]]
        out.append((circular, genes, regions))
    return out


def ir_token(r):
    return f'{r[0]}:{r[1]}' if isinstance(r, tuple) else str(r)


def side_positions(chroms, inp):
    """For each chromosome: the breakpoints (indices of intergenic regions of the
    input genome: region i is after gene i, caps excluded) and, for each block,
    the positions (1-based, counted over all chromosomes) of its genes, None for
    the caps. inp has the chromosomes of the input file, used to find where the
    first block of a circular chromosome starts (it may cross the origin)."""
    bps_chrs, block_pos = [], []
    offset = 0  # genes of the previous chromosomes
    for (circular, blocks, bp_irs), (in_circ, in_genes, in_irs) in zip(chroms, inp):
        assert circular == in_circ
        bps = []
        if not circular:
            k = 0  # index of the next gene, with the caps (0 is the first cap)
            n_genes = sum(len(bg) for bg, _ in blocks) - 2
            for b, (bg, _) in enumerate(blocks):
                pos = []
                for x in bg:
                    pos.append(None if k == 0 or k == n_genes + 1 else offset + k)
                    k += 1
                block_pos.append(pos)
                if b < len(blocks) - 1:
                    # cut between the genes k - 1 and k (with caps) = region k - 1
                    bps.append(k - 1)
        else:
            # genes and regions in the order of the blocks (the region after each gene)
            genes = [x for bg, _ in blocks for x in bg]
            regions = []
            for b, (bg, bi) in enumerate(blocks):
                regions += [ir_token(r) for r in bi]
                if b < len(bp_irs):
                    regions.append(ir_token(bp_irs[b]))
            n_genes = len(genes)
            rot = next(
                (r for r in range(n_genes)
                 if all(in_genes[(r + t) % n_genes] == genes[t] and in_irs[(r + t) % n_genes] == regions[t]
                        for t in range(n_genes))),
                None)
            if rot is None:
                raise ValueError('partition does not match the circular chromosome of the input')
            t = 0
            for b, (bg, _) in enumerate(blocks):
                pos = [offset + (rot + t + d) % n_genes + 1 for d in range(len(bg))]
                block_pos.append(pos)
                t += len(bg)
                if b < len(bp_irs):  # cut after the last gene of the block
                    bps.append(pos[-1] - offset)
        bps_chrs.append(bps)
        offset += n_genes
    return bps_chrs, block_pos


def orientation_circular(bg, bh):
    """Pairing of two whole circular chromosomes: ('d', r) if the target is the
    origin rotated by r (target[t] = origin[t + r]), ('r', r) if it is the reverse
    of the origin rotated by r, None if they are not compatible."""
    (gg, gi), (hg, hi) = bg, bh
    n = len(gg)
    if n != len(hg):
        return None
    rg = [-x for x in reversed(gg)]
    ri = list(reversed(gi[:-1])) + [gi[-1]]  # region after each gene of the reverse
    for r in range(n):
        if all(gg[(t + r) % n] == hg[t] and fits(gi[(t + r) % n], hi[t]) for t in range(n)):
            return ('d', r)
    for r in range(n):
        if all(rg[(t + r) % n] == hg[t] and fits(ri[(t + r) % n], hi[t]) for t in range(n)):
            return ('r', r)
    return None


def fix_instance(lines, inp):
    g = parse_side(lines[0], lines[1], False)
    h = parse_side(lines[2], lines[3], True)
    blocks_g = [blk for _, bl, _ in g for blk in bl]
    blocks_h = [blk for _, bl, _ in h for blk in bl]
    # a block is a whole circular chromosome if its chromosome is circular and has no breakpoint
    whole_g = [circ and not bps for circ, bl, bps in g for _ in bl]
    whole_h = [circ and not bps for circ, bl, bps in h for _ in bl]
    real_g = [i for i, (bg, _) in enumerate(blocks_g) if any(x != 0 for x in bg)]
    real_h = [j for j, (bh, _) in enumerate(blocks_h) if any(x != 0 for x in bh)]
    ori, adj = {}, {i: [] for i in range(len(blocks_g))}
    for i in real_g:
        for j in real_h:
            if whole_g[i] != whole_h[j]:
                continue
            if whole_g[i]:
                o = orientation_circular(blocks_g[i], blocks_h[j])
            else:
                o = orientation(blocks_g[i], blocks_h[j])
                o = (o, 0) if o else None
            if o:
                ori[(i, j)] = o
                adj[i].append(j)
    match = max_matching(len(blocks_g), adj)
    bps_g, pos_g = side_positions(g, inp[0])
    bps_h, pos_h = side_positions(h, inp[1])
    n_g = sum(p is not None for blk in pos_g for p in blk)
    corr = [0] * n_g
    for i, j in match.items():
        pg, ph = pos_g[i], pos_h[j]
        d, r = ori[(i, j)]
        n = len(pg)
        for u, a in enumerate(pg):
            if a is None:
                continue
            if whole_g[i]:
                t = (u - r) % n if d == 'd' else (n - 1 - u - r) % n
            else:
                t = u if d == 'd' else n - 1 - u
            b = ph[t]
            if b is not None:
                corr[a - 1] = b if d == 'd' else -b
    fmt = lambda chrs: ' ; '.join(' '.join(map(str, c)) if c else '-' for c in chrs)
    excl_g = [i for i in real_g if i not in match]
    excl_h = [j for j in real_h if j not in set(match.values())]
    nb_g = sum(len(c) for c in bps_g)
    nb_h = sum(len(c) for c in bps_h)
    stats = (nb_g, nb_h, nb_g + nb_h, len(real_g), len(real_h), len(excl_g), len(excl_h))
    return [fmt(bps_g), fmt(bps_h), ' '.join(map(str, corr))], stats


def split_outputs(path):
    """Outputs of the instances in a PART file (each ends with '# Time')."""
    outs, cur = [], []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            cur.append(line)
            if line.startswith('# Time'):
                outs.append(cur)
                cur = []
    return outs


def main():
    fix_file, cost_file, part_file, input_file = sys.argv[1:5]
    inp = [l.rstrip('\n') for l in open(input_file) if l.strip()]
    with open(fix_file, 'w') as fix, open(cost_file, 'w') as cost:
        for k, lines in enumerate(split_outputs(part_file)):
            time = lines[-1].split(':', 1)[1].strip()
            if len(lines) < 5 or lines[0].startswith('#'):
                fix.write('NA\nNA\nNA\n')
                cost.write('NA\n')
                continue
            sides = (read_input_side(inp[4 * k], inp[4 * k + 1], False), read_input_side(inp[4 * k + 2], inp[4 * k + 3], True))
            fixed, stats = fix_instance(lines[:4], sides)
            fix.write('\n'.join(fixed) + '\n')
            cost.write(' '.join(map(str, stats)) + f' {time}\n')


if __name__ == '__main__':
    main()
