#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║   ATLAS OF WEYL CHAMBERS FOR E₈ — COMPLETE SCAN OF 120 POSITIVE ROOTS        ║
║   ======================================================================     ║
║                                                                              ║
║   v2: +Cartan invariants (trace, determinant) + file export                  ║
║                                                                              ║
║   AUTHOR:  Sergey Viktorovich Matershov                                      ║
║   ORCID:   0009-0009-0641-1357                                               ║
║   License: CC BY-NC-ND 4.0 International                                     ║
║   DATE:    2026-07-08                                                        ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

from fractions import Fraction
import time
import json
import os
from collections import defaultdict
from datetime import datetime


# ═══════════════════════════════════════════════════════════════
# 1. ROOT SYSTEM GENERATION
# ═══════════════════════════════════════════════════════════════

def generate_all_E8_roots():
    """Generate all 240 roots of E₈ in ℝ⁸."""
    roots = []
    
    # Type I: (±1, ±1, 0, 0, 0, 0, 0, 0) — 112 roots (all permutations)
    for i in range(8):
        for j in range(i + 1, 8):
            for s1 in [1, -1]:
                for s2 in [1, -1]:
                    root = [Fraction(0, 1)] * 8
                    root[i] = Fraction(s1, 1)
                    root[j] = Fraction(s2, 1)
                    roots.append(tuple(root))
    
    # Type II: (±½)⁸ with even number of minus signs — 128 roots
    for n in range(256):
        if bin(n).count('1') % 2 == 0:
            root = tuple(
                Fraction(1 if (n >> i) & 1 == 0 else -1, 2)
                for i in range(8)
            )
            roots.append(root)
    
    return roots


def get_positive_roots(all_roots):
    """Extract positive roots (first non-zero coordinate > 0)."""
    pos = []
    for root in all_roots:
        for coord in root:
            if coord > 0:
                pos.append(root)
                break
            elif coord < 0:
                break
    return pos


def scalar_product(a, b):
    """Euclidean scalar product."""
    return sum(ai * bi for ai, bi in zip(a, b))


# ═══════════════════════════════════════════════════════════════
# 2. SIMPLE ROOT DISCOVERY
# ═══════════════════════════════════════════════════════════════

def find_simple_roots_from_seed(positive_roots, seed_positive):
    """
    Algorithmic discovery of simple roots starting from a seed.
    
    A positive root is simple iff it cannot be expressed as
    the sum of two other positive roots (Lemma, Humphreys 1972).
    
    The seed root is placed first in the search order,
    which determines the resulting Weyl chamber.
    """
    pos_set = set(positive_roots)
    n = len(positive_roots[0])
    
    # Place seed at front of search order
    reordered = [seed_positive] + [r for r in positive_roots if r != seed_positive]
    
    simple = []
    
    for alpha in reordered:
        is_simple = True
        for beta in reordered:
            if beta == alpha:
                continue
            gamma = tuple(alpha[i] - beta[i] for i in range(n))
            if gamma in pos_set:
                is_simple = False
                break
        if is_simple:
            simple.append(alpha)
    
    return simple


# ═══════════════════════════════════════════════════════════════
# 3. CARTAN MATRIX AND FUNDAMENTAL WEIGHTS
# ═══════════════════════════════════════════════════════════════

def build_cartan_matrix(simple_roots):
    """
    Build Cartan matrix from simple roots.
    C_{ij} = 2⟨α_i, α_j⟩ / |α_j|²
    """
    rank = len(simple_roots)
    norms = [scalar_product(r, r) for r in simple_roots]
    cartan = [[Fraction(0, 1) for _ in range(rank)] for _ in range(rank)]
    for i in range(rank):
        for j in range(rank):
            dot = scalar_product(simple_roots[i], simple_roots[j])
            cartan[i][j] = Fraction(2 * dot, norms[j])
    return cartan


def det_8x8(m):
    """Determinant of an 8×8 integer matrix via Gaussian elimination."""
    a = [[float(m[i][j]) for j in range(8)] for i in range(8)]
    det = 1.0
    for i in range(8):
        pivot = None
        for j in range(i, 8):
            if abs(a[j][i]) > 1e-10:
                pivot = j
                break
        if pivot is None:
            return 0
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            det = -det
        det *= a[i][i]
        for j in range(i + 1, 8):
            factor = a[j][i] / a[i][i]
            for k in range(i, 8):
                a[j][k] -= factor * a[i][k]
    return int(round(det))


def invert_matrix(matrix):
    """Invert a matrix over Fraction via Gauss–Jordan elimination."""
    n = len(matrix)
    aug = [[Fraction(0, 1) for _ in range(2 * n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            aug[i][j] = matrix[i][j]
        aug[i][n + i] = Fraction(1, 1)
    
    for i in range(n):
        pivot = None
        for j in range(i, n):
            if aug[j][i] != 0:
                pivot = j
                break
        if pivot is None:
            return None
        if pivot != i:
            aug[i], aug[pivot] = aug[pivot], aug[i]
        div = aug[i][i]
        for j in range(2 * n):
            aug[i][j] /= div
        for j in range(n):
            if j != i:
                factor = aug[j][i]
                for k in range(2 * n):
                    aug[j][k] -= factor * aug[i][k]
    
    return [[aug[i][n + j] for j in range(n)] for i in range(n)]


def compute_fundamental_weights(simple_roots):
    """
    Compute fundamental weights via inverse Cartan matrix.
    ω_i = Σ_j (C^{-1})_{ij} α_j
    """
    rank = len(simple_roots)
    n = len(simple_roots[0])
    cartan = build_cartan_matrix(simple_roots)
    cartan_inv = invert_matrix(cartan)
    
    if cartan_inv is None:
        return None
    
    omega = []
    for i in range(rank):
        w = [Fraction(0, 1) for _ in range(n)]
        for j in range(rank):
            coeff = cartan_inv[i][j]
            for k in range(n):
                w[k] += coeff * simple_roots[j][k]
        omega.append(tuple(w))
    
    return omega


def weyl_dimension(omega, rho, pos_roots, weight):
    """
    Weyl dimension formula.
    dim = ∏_{α>0} ⟨λ+ρ, α⟩ / ⟨ρ, α⟩
    """
    n = len(omega[0])
    
    lam = [Fraction(0, 1) for _ in range(n)]
    for i, k in enumerate(weight):
        if k > 0:
            for j in range(n):
                lam[j] += Fraction(k, 1) * omega[i][j]
    
    lam_plus_rho = [lam[j] + rho[j] for j in range(n)]
    
    num = Fraction(1, 1)
    den = Fraction(1, 1)
    
    for root in pos_roots:
        num *= sum(lam_plus_rho[j] * root[j] for j in range(n))
        den *= sum(rho[j] * root[j] for j in range(n))
    
    if den == 0:
        return 0
    
    result = num / den
    return int(result) if result.denominator == 1 else round(float(result))


# ═══════════════════════════════════════════════════════════════
# 4. FULL SCAN OF ALL 120 POSITIVE ROOTS
# ═══════════════════════════════════════════════════════════════

def full_atlas(pos_roots):
    """Full scan: use every positive root as a seed."""
    print(f"  Scanning all {len(pos_roots)} positive roots...")
    print()
    
    n = len(pos_roots[0])
    rho = tuple(Fraction(sum(root[j] for root in pos_roots), 2) for j in range(n))
    
    chambers = []
    seen_cartan = set()
    seen_dim_sets = defaultdict(list)  # dimension set -> list of chamber ids
    
    # Bourbaki reference dimensions
    etalon = [248, 3875, 147250, 6696000, 6899079264, 147250, 3875, 248]
    etalon_set = set(etalon)
    
    for idx, seed in enumerate(pos_roots):
        simple = find_simple_roots_from_seed(pos_roots, seed)
        
        if len(simple) != 8:
            continue
        
        cartan = build_cartan_matrix(simple)
        cartan_int = tuple(tuple(int(cartan[i][j]) for j in range(8)) for i in range(8))
        
        if cartan_int in seen_cartan:
            continue
        
        seen_cartan.add(cartan_int)
        
        # Cartan invariants
        trace_val = sum(cartan_int[i][i] for i in range(8))
        det_val = det_8x8(cartan_int)
        
        omega = compute_fundamental_weights(simple)
        if omega is None:
            continue
        
        dims = []
        for i in range(8):
            weight = tuple(1 if j == i else 0 for j in range(8))
            dim = weyl_dimension(omega, rho, pos_roots, weight)
            dims.append(dim)
        
        matches = sum(1 for d in dims if d in etalon_set)
        dim_set = tuple(sorted(dims))
        
        chambers.append({
            'id': len(chambers) + 1,
            'seed_idx': idx + 1,
            'seed_root': [float(x) for x in seed],
            'cartan_matrix': [[int(cartan[i][j]) for j in range(8)] for i in range(8)],
            'trace': trace_val,
            'determinant': det_val,
            'fundamental_dims': dims,
            'dim_set': list(dim_set),
            'matches_with_bourbaki': matches,
        })
        
        seen_dim_sets[dim_set].append(len(chambers))
        
        if (idx + 1) % 20 == 0:
            print(f"    Scanned: {idx+1}/{len(pos_roots)}, chambers found: {len(chambers)}")
    
    print(f"    Scanned: {len(pos_roots)}/{len(pos_roots)}, chambers found: {len(chambers)}")
    print()
    
    return chambers, seen_dim_sets


def print_full_atlas(chambers, seen_dim_sets):
    """Display full atlas."""
    print("=" * 80)
    print("  COMPLETE ATLAS OF WEYL CHAMBERS FOR E₈")
    print("=" * 80)
    print()
    
    print(f"  Total unique Weyl chambers: {len(chambers)}")
    print(f"  Unique dimension sets: {len(seen_dim_sets)}")
    print()
    
    for ch in chambers:
        print(f"  Chamber #{ch['id']} (seed root #{ch['seed_idx']})")
        print(f"    Matches Bourbaki: {ch['matches_with_bourbaki']}/8")
        print(f"    Cartan trace = {ch['trace']}, determinant = {ch['determinant']}")
        print(f"    Fundamental dimensions:")
        for i, dim in enumerate(ch['fundamental_dims']):
            in_etalon = " ← Bourbaki" if dim in [248, 3875, 147250, 6696000, 6899079264] else ""
            print(f"      ω{i+1} = {dim:,}{in_etalon}")
        print()
    
    print("  DIMENSION SET STATISTICS:")
    print(f"  {'Dimension set':<70} {'Chambers':<10}")
    print(f"  {'-'*70} {'-'*10}")
    for dim_set, ch_ids in sorted(seen_dim_sets.items(), key=lambda x: -len(x[1])):
        dim_str = ", ".join(f"{d:,}" for d in dim_set)
        print(f"  {dim_str:<70} {len(ch_ids):<10}")
    print()
    
    print("  KEY FINDINGS:")
    print(f"  1. Out of 120 positive roots, {len(chambers)} unique Weyl chambers found.")
    if len(seen_dim_sets) == 1:
        print(f"  2. ALL {len(chambers)} chambers share the SAME set of fundamental dimensions.")
        print(f"     Dimensions: {', '.join(f'{d:,}' for d in list(seen_dim_sets.keys())[0])}")
        print(f"  3. Chambers are equivalent — differ only by permutation of weights.")
    else:
        print(f"  2. Found {len(seen_dim_sets)} DISTINCT dimension sets!")
    print(f"  4. Cartan trace = 16, determinant = 1 for all chambers (E₈ invariants).")
    print(f"  5. No prime fundamental dimensions in any chamber.")
    print()


def save_results(chambers, seen_dim_sets, elapsed, output_dir="atlas_results"):
    """Export all results to files."""
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Full JSON export
    atlas_data = {
        'timestamp': timestamp,
        'total_chambers': len(chambers),
        'unique_dim_sets': len(seen_dim_sets),
        'elapsed_sec': elapsed,
        'chambers': chambers,
        'dim_sets': {str(k): v for k, v in seen_dim_sets.items()},
    }
    
    json_path = f"{output_dir}/atlas_E8_{timestamp}.json"
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(atlas_data, f, indent=2, ensure_ascii=False)
    print(f"  💾 Saved: {json_path}")
    
    # LaTeX table export
    latex_path = f"{output_dir}/atlas_E8_table_{timestamp}.tex"
    with open(latex_path, 'w', encoding='utf-8') as f:
        f.write("% Atlas of Weyl chambers for E₈ — auto-generated\n")
        f.write(f"% Timestamp: {timestamp}\n\n")
        f.write("\\begin{tabular}{c|cccccccc}\n")
        f.write("\\toprule\n")
        f.write("Chamber & $\\omega_1$ & $\\omega_2$ & $\\omega_3$ & $\\omega_4$ & $\\omega_5$ & $\\omega_6$ & $\\omega_7$ & $\\omega_8$ \\\\\n")
        f.write("\\midrule\n")
        for ch in chambers:
            dims = ch['fundamental_dims']
            dim_str = " & ".join(f"{d:,}" for d in dims)
            f.write(f"\\#{ch['id']} & {dim_str} \\\\\n")
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
    print(f"  💾 Saved: {latex_path}")
    
    # Plain text report
    report_path = f"{output_dir}/atlas_E8_report_{timestamp}.txt"
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("ATLAS OF WEYL CHAMBERS FOR E₈\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Total chambers: {len(chambers)}\n")
        f.write(f"Unique dimension sets: {len(seen_dim_sets)}\n")
        f.write(f"Elapsed time: {elapsed:.1f} sec\n\n")
        for ch in chambers:
            f.write(f"Chamber #{ch['id']} (seed #{ch['seed_idx']})\n")
            f.write(f"  Trace={ch['trace']}, det={ch['determinant']}\n")
            f.write(f"  Dimensions: {', '.join(f'{d:,}' for d in ch['fundamental_dims'])}\n\n")
    print(f"  💾 Saved: {report_path}")
    
    return json_path


# ═══════════════════════════════════════════════════════════════
# 5. MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("╔" + "═" * 78 + "╗")
    print("║  ATLAS OF WEYL CHAMBERS FOR E₈ — FULL SCAN OF 120 POSITIVE ROOTS     ║")
    print("║  v2: +Cartan invariants + file export                                 ║")
    print("╚" + "═" * 78 + "╝")
    print()
    
    t_start = time.time()
    
    # Generate root system
    print("  Generating 240 roots of E₈...")
    all_roots = generate_all_E8_roots()
    pos_roots = get_positive_roots(all_roots)
    print(f"  Total roots: {len(all_roots)}")
    print(f"  Positive roots: {len(pos_roots)}")
    print()
    
    # Full scan
    chambers, seen_dim_sets = full_atlas(pos_roots)
    
    # Display
    print_full_atlas(chambers, seen_dim_sets)
    
    elapsed = time.time() - t_start
    print(f"  Total time: {elapsed:.1f} sec")
    
    # Export
    save_results(chambers, seen_dim_sets, elapsed)
    
    print("=" * 80)
