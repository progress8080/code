import numpy as np

from shipsupport.model import Block, BlockType, Section, WeightItem
from shipsupport.solvers import BlockSprings
from shipsupport.weights import build_load_model, build_nodes


def block_type(k=1e9, **kw):
    return BlockType("test", contact_width=1.0, contact_length=1.0, stiffness_override=k, **kw)


def make_blocks(xs, k=1e9, gaps=None, **kw):
    t = block_type(k, **kw)
    gaps = [0.0] * len(xs) if gaps is None else gaps
    return [Block(float(x), t, float(g), f"B{i + 1}") for i, (x, g) in enumerate(zip(xs, gaps))]


def uniform_setup(length, q, block_x, k=1e9, EI=1e12, kGA=None, gaps=None, max_len=0.25,
                  extra_items=(), **block_kw):
    """Uniform load q (N/m) over [0, length] on blocks at block_x."""
    items = [WeightItem("uniform", q * length, length / 2, 0.0, length), *extra_items]
    blocks = make_blocks(block_x, k, gaps, **block_kw)
    keys = [0.0, length, *block_x] + [w.xa for w in items] + [w.xb for w in items]
    nodes = build_nodes(keys, 0.0, length, max_len)
    load = build_load_model(items, nodes)
    sections = [Section(0.0, length, EI, kGA)]
    return load, sections, blocks, BlockSprings.from_blocks(blocks)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-30)


def interp_v(result, x):
    return float(np.interp(x, result.nodes, result.v))
