"""Print the fraction of entries per pT bin for samples of a coffea output file.

By default only the mu-enriched QCD b and bb samples (the proxies) are evaluated,
in the pt250msd50mreg50to200 category (the full selected region, no pass/fail split).

Histograms are read from output["variables"][variable][sample][dataset] and have
axes (cat, [variation], variable). Datasets of the same sample are summed, the
requested category (and variation, for MC) is selected and the entries are
integrated over the requested bin ranges. The last bin can extend to "inf", in
which case the overflow is included.

Example:
    python pt_bin_fractions.py output_all.coffea
    python pt_bin_fractions.py output_all.coffea -c msd-100to150_Pt-250toInf_globalParT3_XbbVsQCD-HHbbbb_HP-pass
    python pt_bin_fractions.py output_all.coffea -s all
"""
import argparse

import numpy as np
from coffea.util import load


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", nargs="?", default="output_all.coffea", help="Coffea output file")
    parser.add_argument("-v", "--variable", default="FatJetGood_pt", help="Histogram in output['variables']")
    parser.add_argument("-c", "--category", nargs="+", default=["pt250msd50mreg50to200"], help="Categories to evaluate")
    parser.add_argument("-b", "--bins", nargs="+", type=float, default=[250, 300, 400, 500, np.inf],
                        help="Bin edges, e.g. 250 300 400 500 inf")
    parser.add_argument("--variation", default="nominal", help="Variation to use for MC histograms")
    parser.add_argument("-s", "--samples", nargs="+",
                        default=["QCD_MuEnriched__QCD_MuEnriched_b", "QCD_MuEnriched__QCD_MuEnriched_bb"],
                        help="Exact sample names to evaluate ('all' for every sample). A combined row is added.")
    parser.add_argument("--list-categories", action="store_true", help="List available categories and exit")
    return parser.parse_args()


def edge_index(edges, value):
    """Index of `value` in the flow-including values array (0 = underflow)."""
    if value == -np.inf:
        return 0
    if value == np.inf:
        return len(edges) + 1
    idx = np.flatnonzero(np.isclose(edges, value))
    if len(idx) == 0:
        raise ValueError(f"Bin edge {value} is not a histogram edge. Available edges: {edges}")
    return idx[0] + 1


def project(hists, category, variable, variation):
    """Sum the datasets of one sample and return the 1D (values, variances) including flow."""
    values, variances, edges = None, None, None
    for h in hists.values():
        sel = {"cat": category}
        if "variation" in h.axes.name:
            sel["variation"] = variation
        h1d = h[sel]
        v = h1d.values(flow=True)
        w = h1d.variances(flow=True)
        if values is None:
            values, variances, edges = v.copy(), w.copy(), h1d.axes[variable].edges
        else:
            values += v
            variances += w
    return values, variances, edges


def main():
    args = parse_args()
    output = load(args.input)
    histos = output["variables"][args.variable]

    if args.list_categories:
        h = next(iter(next(iter(histos.values())).values()))
        print("\n".join(h.axes["cat"]))
        return

    if args.samples == ["all"]:
        samples = list(histos.keys())
    else:
        missing = [s for s in args.samples if s not in histos]
        if missing:
            raise KeyError(f"Samples {missing} not found. Available: {list(histos.keys())}")
        samples = args.samples

    ranges = list(zip(args.bins[:-1], args.bins[1:]))
    labels = [f"[{lo:g}, {hi:g}]" for lo, hi in ranges]

    for category in args.category:
        grouped = {}
        for sample in samples:
            values, variances, edges = project(histos[sample], category, args.variable, args.variation)
            grouped[sample] = [values, variances]

        mc_keys = [k for k in grouped if not k.startswith("DATA")]
        if len(mc_keys) > 1:
            name = "Sum (" + " + ".join(k.split("_")[-1] for k in mc_keys) + ")" if args.samples != ["all"] else "Total MC"
            grouped[name] = [sum(grouped[k][0] for k in mc_keys), sum(grouped[k][1] for k in mc_keys)]

        name_width = max(len(k) for k in grouped) + 2
        print(f"\nVariable: {args.variable}   Category: {category}   Variation (MC): {args.variation}")
        header = f"{'Sample':<{name_width}}{'Total':>14}" + "".join(f"{l:>18}" for l in labels)
        print(header)
        print("-" * len(header))
        for key, (values, variances) in grouped.items():
            total = values.sum()
            row = f"{key:<{name_width}}{total:>14.1f}"
            for lo, hi in ranges:
                content = values[edge_index(edges, lo):edge_index(edges, hi)].sum()
                frac = 100 * content / total if total > 0 else np.nan
                row += f"{frac:>17.2f}%"
            print(row)
        print("Total = all entries incl. under/overflow; percentages are relative to it.")


if __name__ == "__main__":
    main()
