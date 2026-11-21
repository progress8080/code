import numpy as np
import pytest

from shipsupport.model import InputError, WeightCurve, WeightItem
from shipsupport.weights import build_load_model, build_nodes, node_index, trapezoid_ends


def curve_load(curve, x_aft, x_fwd, max_len=0.5):
    items = curve.items()
    nodes = build_nodes([w.xa for w in items] + [w.xb for w in items], x_aft, x_fwd, max_len)
    return build_load_model(items, nodes)


def test_trapezoid_preserves_weight_and_lcg():
    item = WeightItem("a", 1000.0, 4.5, 2.0, 8.0)
    wa, wb = trapezoid_ends(item)
    L = item.length
    assert L * (wa + wb) / 2 == pytest.approx(1000.0)
    assert item.xa + L * (wa + 2 * wb) / (3 * (wa + wb)) == pytest.approx(4.5)


def test_trapezoid_limits_are_triangles():
    wa, wb = trapezoid_ends(WeightItem("a", 600.0, 2.0, 0.0, 6.0))   # x̄ = L/3
    assert wb == pytest.approx(0.0, abs=1e-12) and wa == pytest.approx(200.0)
    wa, wb = trapezoid_ends(WeightItem("a", 600.0, 4.0, 0.0, 6.0))   # x̄ = 2L/3
    assert wa == pytest.approx(0.0, abs=1e-12) and wb == pytest.approx(200.0)


def test_lcg_outside_middle_third_rejected():
    with pytest.raises(InputError):
        WeightItem("bad", 100.0, 1.0, 0.0, 6.0)


def test_random_items_integrate_exactly():
    rng = np.random.default_rng(1)
    items = []
    for i in range(30):
        xa = rng.uniform(0, 90)
        L = rng.uniform(1, 30)
        lcg = xa + L * rng.uniform(1 / 3, 2 / 3)
        items.append(WeightItem(f"w{i}", rng.uniform(1e4, 1e6), lcg, xa, xa + L))
    items.append(WeightItem("point", 5e5, 33.3, 33.3, 33.3))
    keys = [w.xa for w in items] + [w.xb for w in items]
    nodes = build_nodes(keys, 0.0, 125.0, 0.7)
    load = build_load_model(items, nodes)

    W = sum(w.weight for w in items)
    M = sum(w.weight * w.lcg for w in items)
    assert load.total == pytest.approx(W, rel=1e-12)
    assert load.moment == pytest.approx(M, rel=1e-12)
    assert np.all(np.diff(nodes) <= 0.7 + 1e-12)


def test_nodes_contain_key_points_and_lookup():
    nodes = build_nodes([3.3, 7.25, 7.25], 0.0, 10.0, 1.0)
    assert node_index(nodes, [0.0, 3.3, 7.25, 10.0]).size == 4
    with pytest.raises(InputError):
        node_index(nodes, 3.31)
    with pytest.raises(InputError):
        build_nodes([11.0], 0.0, 10.0, 1.0)


# --------------------------------------------------------------------------- weight curve

def test_curve_segments_are_exact_trapezoids():
    """Every segment [x_k, x_k+1] becomes one item with weight (q_k + q_k+1) L / 2 and
    LCG x_k + L (q_k + 2 q_k+1) / (3 (q_k + q_k+1)); the load model integrates it exactly."""
    x = [0.0, 10.0, 25.0, 40.0]
    q = [5e3, 20e3, 0.0, 12e3]
    curve = WeightCurve("c", x, q)
    items = curve.items()
    assert [w.name for w in items] == ["c 1/3", "c 2/3", "c 3/3"]
    assert items[0].weight == pytest.approx(12.5e4) and items[0].lcg == pytest.approx(
        10.0 * (5e3 + 40e3) / (3 * 25e3))
    assert items[1].lcg == pytest.approx(10.0 + 15.0 / 3)          # triangle, zero at the right
    assert items[2].lcg == pytest.approx(25.0 + 2 * 15.0 / 3)      # triangle, zero at the left
    # closed form of ∫q dx and ∫x q dx for the piecewise-linear curve
    W = sum(0.5 * (q[k] + q[k + 1]) * (x[k + 1] - x[k]) for k in range(3))
    M = sum((x[k + 1] - x[k]) * (q[k] * (2 * x[k] + x[k + 1]) + q[k + 1] * (x[k] + 2 * x[k + 1])) / 6
            for k in range(3))
    load = curve_load(curve, -5.0, 45.0)
    assert load.total == pytest.approx(W, rel=1e-12)
    assert load.moment == pytest.approx(M, rel=1e-12)
    # the intensity of the load model is the curve itself (and zero outside it)
    probe = np.array([-2.0, 3.0, 17.5, 30.0, 39.0, 43.0])
    expected = np.interp(probe, x, q, left=0.0, right=0.0)
    e = np.searchsorted(load.nodes, probe) - 1                      # element of each probe
    s = (probe - load.nodes[e]) / load.lengths[e]
    assert np.allclose(load.q1[e] + s * (load.q2[e] - load.q1[e]), expected, rtol=1e-12,
                       atol=1e-9)


def test_curve_with_a_jump():
    """The same x twice is a jump: 10 t/m on [0, 20], 30 t/m on [20, 30]."""
    curve = WeightCurve("jump", [0.0, 20.0, 20.0, 30.0], [10.0, 10.0, 30.0, 30.0])
    items = curve.items()
    assert len(items) == 2                                       # the zero-length segment is skipped
    load = curve_load(curve, 0.0, 30.0)
    assert load.total == pytest.approx(10 * 20 + 30 * 10, rel=1e-12)
    assert load.lcg == pytest.approx((200 * 10 + 300 * 25) / 500, rel=1e-12)
    i = int(np.searchsorted(load.nodes, 20.0))
    assert load.q2[i - 1] == pytest.approx(10.0) and load.q1[i] == pytest.approx(30.0)


@pytest.mark.parametrize("x, q, message", [
    ([0.0, 10.0, 5.0], [1.0, 1.0, 1.0], "point 3 \\(x = 5 m\\).*must not be smaller"),
    ([0.0, 10.0, 20.0], [1.0, -2.0, 1.0], "point 2 \\(x = 10 m\\).*negative"),
    ([0.0, float("nan")], [1.0, 1.0], "point 2"),
    ([0.0], [1.0], "at least two points"),
    ([0.0, 10.0], [0.0, 0.0], "carries no weight"),
])
def test_invalid_curve_points_are_named(x, q, message):
    with pytest.raises(InputError, match=message):
        WeightCurve("bad", x, q)
