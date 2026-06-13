#!/usr/bin/env python3
"""
Fetch Allen Mouse Connectivity tracing data via the official REST API.

This script avoids the AllenSDK dependency by using:
  - structure graph download
  - source-search service for injection experiments
  - ProjectionStructureUnionize queries for area-level target weights

Outputs:
  1. allen_tracing_matrix__<metric>__<aggregate>.csv
  2. allen_tracing_long.csv
  3. allen_tracing_experiment_counts.csv
  4. allen_tracing_metadata.json

Recommended metric:
  normalized_projection_volume

References:
  - https://brain-map.org/support/documentation/api-allen-brain-connectivity-atlas
  - https://brain-map.org/support/documentation/connected-services-and-pipes
  - https://api.brain-map.org/doc/ProjectionStructureUnionize.html
"""

from __future__ import annotations

import argparse
import csv
import http.client
import json
import math
import statistics
import sys
import time
import urllib.parse
import urllib.request
import urllib.error
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


API_BASE = "https://api.brain-map.org/api/v2"
PRODUCT_ID = 5  # Mouse Connectivity Projection
GRAPH_ID = 1    # Mouse adult brain ontology


def fetch_json(url: str, pause_s: float = 0.0, timeout_s: float = 60.0, max_retries: int = 5) -> dict:
    if pause_s > 0:
        time.sleep(pause_s)
    last_err = None
    headers = {
        "User-Agent": "snpdc-ibl-allen-fetch/1.0",
        "Accept": "application/json",
        "Connection": "close",
    }
    for attempt in range(max_retries):
        req = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=timeout_s) as resp:
                data = resp.read().decode("utf-8")
            return json.loads(data)
        except (urllib.error.URLError,
                urllib.error.HTTPError,
                http.client.RemoteDisconnected,
                http.client.IncompleteRead,
                TimeoutError) as exc:
            last_err = exc
            if isinstance(exc, urllib.error.HTTPError) and 400 <= exc.code < 500 and exc.code != 429:
                raise
            if attempt == max_retries - 1:
                break
            sleep_s = min(30.0, (2 ** attempt) * 1.5)
            print(f"[allen] request failed, retrying in {sleep_s:.1f}s: {type(exc).__name__}: {exc}", file=sys.stderr)
            time.sleep(sleep_s)
    raise RuntimeError(f"Allen API request failed after {max_retries} attempts: {url}\nLast error: {last_err}")


def fetch_paged_service(criteria: str, pause_s: float = 0.0, page_size: int = 2000) -> List[dict]:
    rows: List[dict] = []
    start = 0
    while True:
        crit = f"{criteria}[start_row$eq{start}][num_rows$eq{page_size}]"
        url = f"{API_BASE}/data/query.json?criteria={urllib.parse.quote(crit, safe='[],$:=')}"
        payload = fetch_json(url, pause_s=pause_s)
        msg = payload.get("msg", [])
        rows.extend(msg)
        total = payload.get("total_rows", len(rows))
        if len(rows) >= total or not msg:
            break
        start += page_size
    return rows


def fetch_paged_model_query(model_path: str, criteria: str, only_fields: Iterable[str], pause_s: float = 0.0,
                            page_size: int = 2000) -> List[dict]:
    rows: List[dict] = []
    start = 0
    only = ",".join(only_fields)
    while True:
        qs = {
            "criteria": criteria,
            "only": only,
            "num_rows": str(page_size),
            "start_row": str(start),
        }
        url = f"{API_BASE}/data/{model_path}/query.json?{urllib.parse.urlencode(qs)}"
        payload = fetch_json(url, pause_s=pause_s)
        msg = payload.get("msg", [])
        rows.extend(msg)
        total = payload.get("total_rows", len(rows))
        if len(rows) >= total or not msg:
            break
        start += page_size
    return rows


def flatten_structure_graph(node: dict, out: List[dict]) -> None:
    rec = {k: node.get(k) for k in ("id", "acronym", "name", "graph_id", "graph_order", "structure_id_path", "parent_structure_id")}
    out.append(rec)
    for child in node.get("children", []) or []:
        flatten_structure_graph(child, out)


def download_structure_graph(pause_s: float = 0.0) -> Tuple[List[dict], Dict[str, int]]:
    url = f"{API_BASE}/structure_graph_download/{GRAPH_ID}.json"
    payload = fetch_json(url, pause_s=pause_s)
    roots = payload["msg"]
    rows: List[dict] = []
    for root in roots:
        flatten_structure_graph(root, rows)
    acronym_to_id = {r["acronym"]: int(r["id"]) for r in rows if r.get("acronym")}
    return rows, acronym_to_id


def read_area_list(path: Path) -> List[str]:
    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text())
        if isinstance(data, dict) and "areas" in data:
            vals = data["areas"]
        else:
            vals = data
        return [str(x).strip() for x in vals if str(x).strip()]
    areas: List[str] = []
    with path.open("r", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        has_comma = "," in sample
        if has_comma:
            reader = csv.reader(f)
            header = next(reader, None)
            if header and len(header) == 1 and header[0].strip().lower() in {"area", "areas", "acronym", "area_name"}:
                pass
            else:
                if header:
                    for item in header:
                        item = item.strip()
                        if item:
                            areas.append(item)
            for row in reader:
                for item in row:
                    item = item.strip()
                    if item:
                        areas.append(item)
        else:
            for line in f:
                item = line.strip()
                if item:
                    areas.append(item)
    return list(dict.fromkeys(areas))


def find_source_experiments(src: str, pause_s: float = 0.0, primary_only: bool = True) -> List[dict]:
    primary = "true" if primary_only else "false"
    crit = (
        f"service::mouse_connectivity_injection_structure"
        f"[injection_structures$eq{src}]"
        f"[primary_structure_only$eq{primary}]"
        f"[product_ids$eq{PRODUCT_ID}]"
        f"[transgenic_lines$eq0]"
    )
    return fetch_paged_service(crit, pause_s=pause_s)


def chunked(seq: List[int], size: int) -> Iterable[List[int]]:
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def fetch_unionizes_for_source(exp_ids: List[int], target_ids: List[int], metrics: List[str], hemisphere_id: int,
                               pause_s: float = 0.0, exp_chunk: int = 20, target_chunk: int = 50) -> List[dict]:
    rows: List[dict] = []
    only_fields = ["section_data_set_id", "structure_id", "hemisphere_id", "is_injection"] + metrics
    for exp_batch in chunked(exp_ids, exp_chunk):
        exp_part = ",".join(str(x) for x in exp_batch)
        for tgt_batch in chunked(target_ids, target_chunk):
            tgt_part = ",".join(str(x) for x in tgt_batch)
            criteria = (
                f"[section_data_set_id$in{exp_part}]"
                f",[structure_id$in{tgt_part}]"
                f",[hemisphere_id$eq{hemisphere_id}]"
                f",[is_injection$eqfalse]"
            )
            batch = fetch_paged_model_query(
                "ProjectionStructureUnionize",
                criteria,
                only_fields=only_fields,
                pause_s=pause_s,
                page_size=2000,
            )
            rows.extend(batch)
    return rows


def safe_float(x) -> float:
    try:
        v = float(x)
    except Exception:
        return math.nan
    return v


def aggregate_source_target(unionizes: List[dict], metrics: List[str]) -> Dict[str, Dict[Tuple[int, int], float]]:
    by_pair: Dict[str, Dict[Tuple[int, int], List[float]]] = {m: {} for m in metrics}
    for row in unionizes:
        exp_id = int(row["section_data_set_id"])
        tgt_id = int(row["structure_id"])
        for metric in metrics:
            val = safe_float(row.get(metric))
            if not math.isfinite(val):
                continue
            by_pair[metric].setdefault((exp_id, tgt_id), []).append(val)

    out: Dict[str, Dict[Tuple[int, int], float]] = {}
    for metric, pairs in by_pair.items():
        out[metric] = {}
        for key, vals in pairs.items():
            out[metric][key] = statistics.mean(vals)
    return out


def aggregate_across_experiments(vals: List[float], aggregate: str) -> float:
    if not vals:
        return math.nan
    aggregate = aggregate.lower()
    if aggregate == "median":
        return statistics.median(vals)
    if aggregate == "mean":
        return statistics.mean(vals)
    if aggregate == "max":
        return max(vals)
    if aggregate == "support_frac":
        return sum(v > 0 for v in vals) / len(vals)
    raise ValueError(f"Unknown aggregate: {aggregate}")


def slug(name: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in name).strip("_")


@dataclass
class SourceSummary:
    acronym: str
    structure_id: int
    experiment_ids: List[int]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--areas-file", required=True, help="Text/CSV/JSON file listing acronyms to compare.")
    ap.add_argument("--out-dir", required=True, help="Output directory for Allen tracing files.")
    ap.add_argument("--metrics", default="normalized_projection_volume,projection_density,projection_energy",
                    help="Comma-separated Allen unionize metrics.")
    ap.add_argument("--aggregates", default="median,mean,max,support_frac",
                    help="Comma-separated aggregation methods across experiments.")
    ap.add_argument("--hemisphere-id", type=int, default=3, help="Allen hemisphere id. 3 means both hemispheres.")
    ap.add_argument("--pause-s", type=float, default=0.1, help="Pause between API requests.")
    ap.add_argument("--primary-only", action="store_true", default=False,
                    help="Restrict source experiments to primary injection structure only.")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    areas = read_area_list(Path(args.areas_file))
    if not areas:
        raise SystemExit("No areas found in areas file.")

    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]
    valid_metrics = {"normalized_projection_volume", "projection_density", "projection_energy"}
    bad_metrics = [m for m in metrics if m not in valid_metrics]
    if bad_metrics:
        raise SystemExit(f"Unsupported metrics: {bad_metrics}")

    aggregates = [a.strip() for a in args.aggregates.split(",") if a.strip()]
    valid_aggs = {"median", "mean", "max", "support_frac"}
    bad_aggs = [a for a in aggregates if a not in valid_aggs]
    if bad_aggs:
        raise SystemExit(f"Unsupported aggregates: {bad_aggs}")

    print(f"[allen] requested areas: {len(areas)}")
    structures, acronym_to_id = download_structure_graph(pause_s=args.pause_s)
    print(f"[allen] structure graph loaded: {len(structures)} structures")

    valid_areas = [a for a in areas if a in acronym_to_id]
    missing_areas = [a for a in areas if a not in acronym_to_id]
    target_ids = [acronym_to_id[a] for a in valid_areas]

    source_summaries: List[SourceSummary] = []
    experiment_counts: List[dict] = []

    matrices: Dict[Tuple[str, str], Dict[Tuple[str, str], float]] = {
        (metric, agg): {} for metric in metrics for agg in aggregates
    }
    long_rows: List[dict] = []

    for src in valid_areas:
        src_id = acronym_to_id[src]
        experiments = find_source_experiments(src, pause_s=args.pause_s, primary_only=args.primary_only)
        exp_ids = sorted({int(r["id"]) for r in experiments if "id" in r})
        source_summaries.append(SourceSummary(src, src_id, exp_ids))
        experiment_counts.append({
            "source_acronym": src,
            "source_structure_id": src_id,
            "n_experiments": len(exp_ids),
        })
        print(f"[allen] source {src}: {len(exp_ids)} experiments")
        if not exp_ids:
            continue

        unionizes = fetch_unionizes_for_source(
            exp_ids, target_ids, metrics, args.hemisphere_id,
            pause_s=args.pause_s,
        )
        exp_target_by_metric = aggregate_source_target(unionizes, metrics)

        for tgt in valid_areas:
            tgt_id = acronym_to_id[tgt]
            for metric in metrics:
                exp_target = exp_target_by_metric[metric]
                vals = [val for (exp_id, tid), val in exp_target.items() if tid == tgt_id]
                n_support = len(vals)
                for agg in aggregates:
                    weight = aggregate_across_experiments(vals, agg)
                    matrices[(metric, agg)][(tgt, src)] = weight
                    long_rows.append({
                        "target_acronym": tgt,
                        "target_structure_id": tgt_id,
                        "source_acronym": src,
                        "source_structure_id": src_id,
                        "metric": metric,
                        "aggregate": agg,
                        "weight": weight,
                        "n_supporting_experiments": n_support,
                    })

    written_matrices = []
    for metric in metrics:
        for agg in aggregates:
            matrix_csv = out_dir / f"allen_tracing_matrix__{slug(metric)}__{slug(agg)}.csv"
            with matrix_csv.open("w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["target_acronym"] + valid_areas)
                for tgt in valid_areas:
                    row = [tgt]
                    for src in valid_areas:
                        v = matrices[(metric, agg)].get((tgt, src), math.nan)
                        row.append("" if not math.isfinite(v) else f"{v:.17g}")
                    writer.writerow(row)
            written_matrices.append(matrix_csv.name)

    # Backward-compatible default path for the main recommended choice.
    default_metric = metrics[0]
    default_agg = aggregates[0]
    legacy_csv = out_dir / "allen_tracing_matrix.csv"
    legacy_src = out_dir / f"allen_tracing_matrix__{slug(default_metric)}__{slug(default_agg)}.csv"
    legacy_csv.write_text(legacy_src.read_text())

    long_csv = out_dir / "allen_tracing_long.csv"
    with long_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(long_rows[0].keys()) if long_rows else [
            "target_acronym", "target_structure_id", "source_acronym",
            "source_structure_id", "metric", "aggregate", "weight", "n_supporting_experiments"
        ])
        writer.writeheader()
        writer.writerows(long_rows)

    counts_csv = out_dir / "allen_tracing_experiment_counts.csv"
    with counts_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["source_acronym", "source_structure_id", "n_experiments"])
        writer.writeheader()
        writer.writerows(experiment_counts)

    meta = {
        "metrics": metrics,
        "aggregates": aggregates,
        "hemisphere_id": args.hemisphere_id,
        "primary_only": args.primary_only,
        "product_id": PRODUCT_ID,
        "graph_id": GRAPH_ID,
        "requested_areas": areas,
        "valid_areas": valid_areas,
        "missing_areas": missing_areas,
        "n_requested_areas": len(areas),
        "n_valid_areas": len(valid_areas),
        "n_missing_areas": len(missing_areas),
        "api_base": API_BASE,
        "written_matrices": written_matrices,
        "legacy_matrix": legacy_csv.name,
    }
    (out_dir / "allen_tracing_metadata.json").write_text(json.dumps(meta, indent=2))

    print(f"[allen] valid areas: {len(valid_areas)}")
    print(f"[allen] missing areas: {len(missing_areas)}")
    print(f"[allen] wrote {len(written_matrices)} matrix files to: {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
