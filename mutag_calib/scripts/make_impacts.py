"""Make combine impact plots for one fit folder.

Run inside a folder containing workspace.root (and ideally fit_status.json):

    python3 make_impacts.py [--parallel 8] [--mcstat-total]

Fit options (POIs, frozen parameters, parameter ranges, minimizer strategy) are
taken from the FitDiagnostics command stored in fit_status.json, so the impacts
are consistent with the fit.

Steps:
  1. initial fit                          (combineTool.py -M Impacts --doInitialFit)
  2. one fit per nuisance                 (--doFits), including every per-bin
                                          autoMCStats nuisance (prop_bin*)
  3. collect into impacts.json; nuisance fits that failed are retried once with
     more robust minimizer settings (and dropped with a warning if still failing)
  4. only with --mcstat-total: the prop_bin* entries are replaced by one entry
     "autoMCStats_total" (fit with the whole autoMCStats group frozen at the
     best-fit snapshot, subtracted in quadrature). This only changes the plot,
     not the fit.
  5. plot impacts_<POI>.pdf for every POI
"""
import argparse
import json
import math
import os
import shlex
import subprocess
import sys

import ROOT

ROOT.gROOT.SetBatch(True)

MASS = "120"
# DEFAULT_FIT_OPTS = {
#     "--redefineSignalPOIs": "SF_b,SF_c,SF_light",
#     "--setParameters": "r=1,SF_light=1",
#     "--freezeParameters": "r,SF_light",
#     "--setParameterRanges": None,
#     "--cminDefaultMinimizerStrategy": "0",
# }
DEFAULT_FIT_OPTS = {
    "--redefineSignalPOIs": "SF_b,SF_c_fail,SF_c_pass,SF_light",
    "--setParameters": "r=1,SF_light=1,SF_c_pass=1",
    "--freezeParameters": "r,SF_light,SF_c_pass",
    "--setParameterRanges": None,
    "--cminDefaultMinimizerStrategy": "0",
}
ROBUST_OPTS = [
    "--cminDefaultMinimizerStrategy", "1",
    "--X-rtd", "MINIMIZER_analytic",
    "--cminFallbackAlgo", "Minuit2,Migrad,0:0.1",
]
MCSTAT_ENTRY = "autoMCStats_total"


def run(cmd, log=None):
    print(">> " + " ".join(shlex.quote(c) for c in cmd), flush=True)
    if log is None:
        return subprocess.run(cmd).returncode
    with open(log, "a") as f:
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT).returncode


def read_fit_options():
    opts = dict(DEFAULT_FIT_OPTS)
    if not os.path.exists("fit_status.json"):
        print("fit_status.json not found, using default fit options")
        return opts
    tokens = shlex.split(json.load(open("fit_status.json"))["command"])
    for i, tok in enumerate(tokens[:-1]):
        if tok in opts:
            opts[tok] = tokens[i + 1]
    return opts


def build_options(fit_opts):
    frozen = [p for p in (fit_opts["--freezeParameters"] or "").split(",") if p]
    pois = [p for p in fit_opts["--redefineSignalPOIs"].split(",") if p and p not in frozen]
    common = ["--redefineSignalPOIs", ",".join(pois)]
    for key in ("--setParameters", "--freezeParameters", "--setParameterRanges"):
        if fit_opts[key]:
            common += [key, fit_opts[key]]
    strategy = ["--cminDefaultMinimizerStrategy", fit_opts["--cminDefaultMinimizerStrategy"] or "0"]
    return pois, frozen, common, strategy


def impacts_cmd(common, *extra):
    return ["combineTool.py", "-M", "Impacts", "-d", "workspace.root", "-m", MASS] + common + list(extra)


def n_entries(fname):
    if not os.path.exists(fname):
        return 0
    f = ROOT.TFile.Open(fname)
    if not f or f.IsZombie():
        return 0
    tree = f.Get("limit")
    n = tree.GetEntries() if tree else 0
    f.Close()
    return n


def collect(exclude):
    """Collect into impacts.json; return the nuisances whose fits are missing or failed."""
    if os.path.exists("impacts.json"):
        os.remove("impacts.json")
    out = subprocess.run(impacts_cmd(COMMON, *STRATEGY, "--exclude", exclude, "-o", "impacts.json"),
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith("Missing inputs: "):
            return line[len("Missing inputs: "):].split(",")
    return []


def singles_intervals(fname, pois):
    """Return {poi: (best, lo, hi)} from a MultiDimFit --algo singles output."""
    f = ROOT.TFile.Open(fname)
    tree = f.Get("limit")
    vals = {p: [] for p in pois}
    for entry in tree:
        for p in pois:
            vals[p].append(getattr(entry, p))
    f.Close()
    # entry 0 is the best fit, then (lo, hi) pairs, one per POI in order
    res = {}
    for i, p in enumerate(pois):
        best = vals[p][0]
        lo, hi = sorted(vals[p][1 + 2 * i:3 + 2 * i])
        res[p] = (best, lo, hi)
    return res


def mcstat_total_impact(pois):
    """Impact of all autoMCStats nuisances together, via freezing at best fit."""
    print("\n=== Total autoMCStats impact ===")
    # Total intervals come from the step 1 initial fit
    nominal_file = "higgsCombine_initialFit_Test.MultiDimFit.mH%s.root" % MASS
    snapshot_file = "higgsCombine_mcstatSnapshot.MultiDimFit.mH%s.root" % MASS
    frozen_file = "higgsCombine_mcstatFrozen.MultiDimFit.mH%s.root" % MASS
    # The frozen interval scan is sensitive to the exact starting point, so try
    # combinations of minimizer settings for the snapshot and the frozen fit
    ok = False
    for snap_min in (STRATEGY, ROBUST_OPTS):
        run(["combine", "-M", "MultiDimFit", "-d", "workspace.root", "-m", MASS, "-n", "_mcstatSnapshot",
             "--algo", "none", "--saveWorkspace"] + COMMON + snap_min, log="impacts_mcstat.log")
        for frozen_min in (STRATEGY, ROBUST_OPTS):
            with open("impacts_mcstat.log", "a") as f:
                start = f.tell()
            run(["combine", "-M", "MultiDimFit", "-d", snapshot_file, "-m", MASS, "-n", "_mcstatFrozen",
                 "--algo", "singles", "--robustFit", "1", "--snapshotName", "MultiDimFit",
                 "--freezeNuisanceGroups", "autoMCStats"] + COMMON + frozen_min, log="impacts_mcstat.log")
            with open("impacts_mcstat.log") as f:
                f.seek(start)
                if "No valid" not in f.read():
                    ok = True
                    break
        if ok:
            break
    if not ok:
        print("WARNING: an interval scan hit the parameter range in the autoMCStats fits "
              "(see impacts_mcstat.log); no total entry added")
        return None
    if n_entries(nominal_file) < 1 + 2 * len(pois) or n_entries(frozen_file) < 1 + 2 * len(pois):
        print("WARNING: autoMCStats fits failed, see impacts_mcstat.log; no total entry added")
        return None
    nom = singles_intervals(nominal_file, pois)
    frz = singles_intervals(frozen_file, pois)
    impacts = {}
    for p in pois:
        best, lo, hi = nom[p]
        _, flo, fhi = frz[p]
        down = math.sqrt(max((best - lo) ** 2 - (best - flo) ** 2, 0.0))
        up = math.sqrt(max((hi - best) ** 2 - (fhi - best) ** 2, 0.0))
        impacts[p] = (best, down, up)
        print("  %-8s = %.4f  total -%.4f/+%.4f  | MC-stat -%.4f/+%.4f" % (p, best, best - lo, hi - best, down, up))
    return impacts


def add_mcstat_entry(impacts):
    data = json.load(open("impacts.json"))
    data["params"] = [p for p in data["params"] if p["name"] != MCSTAT_ENTRY]
    # "Unconstrained" so plotImpacts draws no pull for this pseudo-nuisance
    entry = {"name": MCSTAT_ENTRY, "type": "Unconstrained", "groups": ["autoMCStats"],
             "prefit": [-1.0, 0.0, 1.0], "fit": [-1.0, 0.0, 1.0]}
    for p, (best, down, up) in impacts.items():
        entry[p] = [best - down, best, best + up]
        entry["impact_" + p] = max(down, up)
    data["params"].append(entry)
    json.dump(data, open("impacts.json", "w"), indent=2)


def main():
    global COMMON, STRATEGY
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--parallel", type=int, default=8, help="parallel nuisance fits (default 8)")
    parser.add_argument("--mcstat-total", action="store_true",
                        help="show one combined autoMCStats entry instead of one per bin (prop_bin*)")
    args = parser.parse_args()

    if not os.path.exists("workspace.root"):
        sys.exit("No workspace.root in %s" % os.getcwd())

    pois, frozen, COMMON, STRATEGY = build_options(read_fit_options())
    exclude = ",".join(frozen + (["rgx{prop_bin.*}"] if args.mcstat_total else []))
    print("POIs: %s | frozen: %s | exclude: %s" % (pois, frozen, exclude))
    for log in ("impacts_initialFit.log", "impacts_doFits.log", "impacts_mcstat.log"):
        if os.path.exists(log):
            os.remove(log)

    print("\n=== Step 1: initial fit ===")
    run(impacts_cmd(COMMON, *STRATEGY, "--doInitialFit", "--robustFit", "1"), log="impacts_initialFit.log")
    os.system("grep -E 'ERROR|Warning|^ +SF_' impacts_initialFit.log")

    print("\n=== Step 2: nuisance fits ===")
    run(impacts_cmd(COMMON, *STRATEGY, "--doFits", "--robustFit", "1",
                    "--parallel", str(args.parallel), "--exclude", exclude), log="impacts_doFits.log")

    print("\n=== Step 3+4: collect, retry failed fits ===")
    failed = collect(exclude)
    if failed:
        print("Retrying %d failed fits: %s" % (len(failed), ", ".join(failed)))
        run(impacts_cmd(COMMON, *ROBUST_OPTS, "--doFits", "--robustFit", "1",
                        "--parallel", str(args.parallel), "--named", ",".join(failed)), log="impacts_doFits.log")
        still = collect(exclude)
        if still:
            print("WARNING: still failing, excluded from impacts: %s" % ", ".join(still))
            exclude = ",".join([exclude] + still)
            collect(exclude)
    else:
        print("All nuisance fits OK")
    if not os.path.exists("impacts.json"):
        sys.exit("impacts.json was not produced, check the output above")

    if args.mcstat_total:
        impacts = mcstat_total_impact(pois)
        if impacts:
            add_mcstat_entry(impacts)

    print("\n=== Step 5: plots ===")
    for poi in pois:
        run(["plotImpacts.py", "-i", "impacts.json", "-o", "impacts_%s" % poi, "--POI", poi])
    print("\nDone: " + ", ".join("impacts_%s.pdf" % p for p in pois))


if __name__ == "__main__":
    main()
