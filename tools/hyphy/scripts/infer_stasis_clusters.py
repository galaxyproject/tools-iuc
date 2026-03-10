#!/usr/bin/env python3
"""
B-STILL Stasis Cluster Inference Tool
====================================
Identifies regional footprints of extreme purifying selection (stasis) in B-STILL 
JSON results using a FWER-controlled Hypergeometric Scan Statistic.

Usage:
    python3 infer_stasis_clusters.py input.json --ebf 10 --permutations 10000 --output results.json
"""

import sys
import json
import argparse
import time
import numpy as np
from scipy.stats import hypergeom

def get_sf_optimized(n, d, L, K, cache):
    """Retrieves or computes Hypergeometric Survival Function value."""
    key = (n, d)
    if key not in cache:
        cache[key] = hypergeom.sf(n - 1, L, K, d)
    return cache[key]

def scan_intervals(indices, L, K, max_size, sf_cache, threshold=None):
    """
    Scans all possible intervals [i, j] anchored by stasis events.
    Returns the minimum p-value if threshold is None, else returns all significant segments.
    """
    best_p = 1.0
    segments = []
    num_events = len(indices)
    
    for n in range(3, min(max_size + 1, num_events + 1)):
        for i in range(num_events - n + 1):
            d = indices[i + n - 1] - indices[i] + 1
            p = get_sf_optimized(n, d, L, K, sf_cache)
            
            if threshold is None:
                if p < best_p: best_p = p
            else:
                if p <= threshold:
                    segments.append({
                        "start": int(indices[i] + 1), 
                        "end": int(indices[i + n - 1] + 1),
                        "p_value": p,
                        "k": n,
                        "d": int(d)
                    })
    
    return best_p if threshold is None else segments

def merge_segments(segments, merge_dist=15):
    """Merges overlapping or nearby significant segments."""
    if not segments: return []
    segments.sort(key=lambda x: x['start'])
    
    merged = []
    curr = segments[0]
    for next_s in segments[1:]:
        if next_s['start'] <= curr['end'] + merge_dist:
            curr['end'] = max(curr['end'], next_s['end'])
            curr['p_value'] = min(curr['p_value'], next_s['p_value'])
            curr['d'] = curr['end'] - curr['start'] + 1
        else:
            merged.append(curr)
            curr = next_s
    merged.append(curr)
    return merged

def main():
    parser = argparse.ArgumentParser(description="Infer stasis clusters from B-STILL JSON.")
    parser.add_argument("input", help="Path to B-STILL JSON result file")
    parser.add_argument("--ebf", type=float, default=10.0, help="EBF threshold for defining stasis sites (default: 10.0)")
    parser.add_argument("--permutations", type=int, default=10000, help="Number of permutations for FWER control (default: 10000)")
    parser.add_argument("--alpha", type=float, default=0.05, help="Family-wise error rate threshold (default: 0.05)")
    parser.add_argument("--max-cluster", type=int, default=30, help="Maximum number of stasis sites per interval scan (default: 30)")
    parser.add_argument("--merge", type=int, default=15, help="Distance in codons to merge adjacent clusters (default: 15)")
    parser.add_argument("--output", help="Path to save results in JSON format")
    
    args = parser.parse_args()

    try:
        with open(args.input, "r") as f:
            data = json.load(f)
    except Exception as e:
        print("Error loading JSON: {0}".format(e))
        sys.exit(1)

    sites = data.get("MLE", {}).get("content", {}).get("0", [])
    ebfs = [s[12] if (len(s) > 12 and isinstance(s[12], (int, float))) else 0 for s in sites]
    L = len(ebfs)
    
    if L < 10:
        print("Alignment too short for cluster analysis.")
        sys.exit(0)

    stasis_indices = np.array([i for i, val in enumerate(ebfs) if val >= args.ebf])
    K = len(stasis_indices)
    
    print("--- B-STILL Cluster Inference ---")
    print("Input: {0}".format(args.input))
    print("Gene Length (L): {0} codons".format(L))
    print("Stasis Sites (K): {0} (EBF >= {1})".format(K, args.ebf))
    
    if K < 3:
        print("Insufficient stasis sites to form clusters (minimum 3 required).")
        sys.exit(0)

    print("Running {0} permutations for FWER control...".format(args.permutations))
    null_min_ps = []
    all_positions = np.arange(L)
    sf_cache = {}
    
    start_time = time.time()
    for i in range(args.permutations):
        if i > 0 and i % 1000 == 0:
            elapsed = time.time() - start_time
            print("  Processed {0} permutations... ({1:.1f} per sec)".format(i, i/elapsed))
            
        shuffled = sorted(np.random.choice(all_positions, K, replace=False))
        min_p = scan_intervals(shuffled, L, K, args.max_cluster, sf_cache)
        null_min_ps.append(min_p)
    
    crit_p = np.percentile(null_min_ps, args.alpha * 100)
    print("Gene-specific Critical P-value (FWER {0}): {1:.2e}".format(args.alpha, crit_p))

    print("Scanning observed sequence for significant clusters...")
    raw_segments = scan_intervals(stasis_indices, L, K, args.max_cluster, sf_cache, threshold=crit_p)
    
    final_clusters = merge_segments(raw_segments, merge_dist=args.merge)
    
    for c in final_clusters:
        c['k'] = sum(1 for idx in stasis_indices if c['start'] <= idx+1 <= c['end'])

    print("\nFound {0} significant stasis clusters:".format(len(final_clusters)))
    if final_clusters:
        print("\nLegend:")
        print("  k : Number of high-confidence stasis sites within the cluster")
        print("  d : Total span of the cluster in codons")
        print("\n{:<8} | {:<8} | {:<5} | {:<5} | {:<10}".format("Start", "End", "k", "d", "P-value"))
        print("-" * 45)
        for c in final_clusters:
            print("{:<8} | {:<8} | {:<5} | {:<5} | {:.2e}".format(c['start'], c['end'], c['k'], c['d'], c['p_value']))

    if args.output:
        output_data = {
            "input_file": args.input,
            "parameters": vars(args),
            "summary": {
                "gene_length": L,
                "total_stasis_sites": K,
                "critical_p_value": float(crit_p),
                "num_clusters": len(final_clusters)
            },
            "clusters": final_clusters
        }
        with open(args.output, "w") as f:
            json.dump(output_data, f, indent=4)
        print("\nDetailed results saved to {0}".format(args.output))

if __name__ == "__main__":
    main()
