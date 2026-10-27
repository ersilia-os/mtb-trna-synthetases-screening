#!/usr/bin/env python3
"""
Builds the compound list for the Enamine order: cuts script 101's selection lists at a per-list
top-N, prices the result from an Enamine quote, and flags which compounds to buy.

1. From script 101's `merged_selections.csv`, keeps the entries ranked <= `--top-n` (single and
   dual lists) or <= `--top-n-multi` (multi lists) -- on the ranks as they are, no re-ranking
   after the price rule. Script 101's 4 hand-picked near-misses (ranks 11-12) fall outside any
   `--top-n` <= 10 and drop out automatically.
2. Joins the quote's catalog ID, collection, price, lead time and purity on compound ID. Every
   selected compound must be in the quote, otherwise the script stops.
3. `to_buy` = quote price < `--price-cutoff`. Non-buyable compounds are dropped from the output;
   `--keep-non-buyable` keeps them in the table, flagged False.
4. `origins`, `n_selections`, `best_rank`, `selected_targets` are recomputed from the kept
   entries only (script 101's own versions describe the uncut top-10/top-50 selection); the
   remaining columns (physchem/ADMET, on-target scores, off-target selectivity) are script
   101's `audit_input.csv`, unchanged.

The quote's prices are confidential (Enamine quote terms), and `output/` is eosvc-synced and
public -- so the quote path and the output path are arguments, and the script refuses an output
path inside the repository.

Usage:
    python 102_order_selection.py --quote <quote.xlsx> --output <outside-repo>/selection.csv \
        --price-cutoff <eur>
"""
import argparse
import os
import sys

import pandas as pd

ROOT = os.path.realpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

MERGED_CSV = os.path.join(ROOT, "output", "101_merge_selections", "merged_selections.csv")
AUDIT_CSV = os.path.join(ROOT, "output", "101_merge_selections", "audit_input.csv")

TOP_N = 7
TOP_N_MULTI = 35
PRICE_COLUMN = "Price per 3 mg, EUR"

QUOTE_RENAME = {"Customer Code": "compound_id", "Catalog ID": "enamine_catalog_id",
                "Collection": "enamine_collection", "Lead time": "enamine_lead_time",
                "Purity, %": "enamine_purity_pct"}
SCRIPT_101_PROVENANCE = ["origins", "n_selections", "best_rank", "hit_types_catalytic",
                         "hit_types_noncatalytic", "selected_targets", "starred_columns", "prioritized"]


def cut_selections(merged, top_n, top_n_multi):
    """Entries within the per-list top-N (`top_n_multi` for multi lists)."""
    limit = merged["hit_type"].map(lambda t: top_n_multi if t == "multi" else top_n)
    return merged[merged["rank"] <= limit]


def read_quote(path, price_column):
    """One row per compound in the quote, with prices rounded to the cent."""
    quote = pd.read_excel(path).rename(columns=QUOTE_RENAME)
    quote = quote[quote["compound_id"].notna()]
    if price_column not in quote.columns:
        sys.exit(f"Price column '{price_column}' not in the quote, found: {list(quote.columns)}")
    quote["enamine_price_eur"] = quote[price_column].astype(float).round(2)
    if not quote["compound_id"].is_unique:
        sys.exit("Duplicate customer codes in the quote")
    return quote[["compound_id", "enamine_catalog_id", "enamine_collection", "enamine_price_eur",
                  "enamine_lead_time", "enamine_purity_pct"]]


def provenance(kept):
    """Per-compound summary of the kept entries."""
    grouped = kept.groupby("compound_id")
    return pd.DataFrame({
        "origins": grouped["origin"].agg(lambda s: "|".join(sorted(set(s)))),
        "n_selections": grouped.size(),
        "best_rank": grouped["rank"].min(),
        "selected_targets": grouped["targets"].agg(lambda s: "|".join(sorted({t for v in s for t in v.split("|")}))),
    }).reset_index()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--quote", required=True, help="Enamine quote xlsx (confidential)")
    parser.add_argument("--output", required=True, help="output CSV, must be outside the repository")
    parser.add_argument("--top-n", type=int, default=TOP_N)
    parser.add_argument("--top-n-multi", type=int, default=TOP_N_MULTI)
    parser.add_argument("--price-cutoff", type=float, required=True, help="to_buy = price < cutoff (EUR)")
    parser.add_argument("--price-column", default=PRICE_COLUMN)
    parser.add_argument("--keep-non-buyable", action="store_true", help="keep to_buy = False rows in the output")
    args = parser.parse_args()

    out_path = os.path.realpath(os.path.expanduser(args.output))
    if os.path.commonpath([ROOT, out_path]) == ROOT:
        sys.exit(f"Refusing to write quote prices inside the repository: {out_path}")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    kept = cut_selections(pd.read_csv(MERGED_CSV), args.top_n, args.top_n_multi)
    quote = read_quote(os.path.expanduser(args.quote), args.price_column)

    missing = sorted(set(kept["compound_id"]) - set(quote["compound_id"]))
    if missing:
        sys.exit(f"{len(missing)} selected compounds are not in the quote, e.g. {missing[:3]}")
    unpriced = quote[quote["compound_id"].isin(kept["compound_id"]) & quote["enamine_price_eur"].isna()]
    if len(unpriced):
        sys.exit(f"{len(unpriced)} selected compounds have no price in the quote, e.g. {unpriced['compound_id'].tolist()[:3]}")

    audit = pd.read_csv(AUDIT_CSV).drop(columns=SCRIPT_101_PROVENANCE)
    table = (quote[quote["compound_id"].isin(kept["compound_id"])]
             .merge(provenance(kept), on="compound_id", validate="1:1")
             .merge(audit, on="compound_id", validate="1:1"))
    table.insert(table.columns.get_loc("enamine_price_eur") + 1, "to_buy",
                 table["enamine_price_eur"] < args.price_cutoff)
    first = ["compound_id", "smiles"]
    table = table[first + [c for c in table.columns if c not in first]]

    buy = table[table["to_buy"]]
    print(f"top-{args.top_n} per list, top-{args.top_n_multi} per multi list, to_buy = price < {args.price_cutoff:g} EUR")
    print(f"Selected compounds: {len(table):,} | to_buy: {len(buy):,} | not to_buy: {len(table) - len(buy):,}")
    print(f"to_buy compound cost: {buy['enamine_price_eur'].sum():,.2f} EUR")
    print(f"to_buy collections: {buy['enamine_collection'].value_counts().to_dict()}")
    (table if args.keep_non_buyable else buy).to_csv(out_path, index=False)
    print(f"Saved {len(table) if args.keep_non_buyable else len(buy):,} rows -> {out_path}")


if __name__ == "__main__":
    main()
