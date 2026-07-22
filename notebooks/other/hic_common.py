"""Shared settings and helpers for the compartment / insulation notebooks.

The number of processes and the number of threads per process are NOT set
here: every notebook chooses them itself and passes them to run_parallel and
to the workers.
"""

import os

# Math libraries are pinned to one thread, because parallelism comes from the
# process pool. This has to happen before numpy is imported.
for _variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                  "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import re
import gc
import glob
import multiprocessing
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

import cooler
import cooltools
import bioframe


# ---------------------------------------------------------------- folders ---

RAW_ROOT = Path("../raw")                       # mcool files
DATA_ROOT = Path("../data")
COMPARTMENT_CACHE = DATA_ROOT / "compartments"  # E1 and saddle per cooler
INSULATION_CACHE = DATA_ROOT / "insulation"     # insulation table per cooler
DATASETS_DIR = DATA_ROOT / "datasets"           # summary tables
TABLES_DIR = Path("../tables")
FIGURES_DIR = Path("../figures")

GC_COV_FILE = Path("other/files/hg38_gc_cov_100kb.tsv")

SUFFIX = "hg38.mapq_30.1000"                    # end of every mcool file name


# --------------------------------------------------------------- settings ---

RES_COMPARTMENTS = 100_000     # resolution for E1 and saddle
RES_INSULATION = 10_000        # resolution for insulation

INSUL_WINDOWS = [50_000, 100_000]   # only the two windows used in the paper
TOLERANCES = [1, 2, 3, 5]           # boundary matching tolerance, in bins

MAIN_WINDOW = 100_000          # window of figure 5 and table S5
MAIN_TOLERANCE = 2             # tolerance of figure 5 and table S5

DEPTHS = [200_000_000, 100_000_000, 50_000_000, 25_000_000, 10_000_000, 1_000_000]
FULL_DEPTH_READS = 350_000_000  # nominal depth of a full library, used as x value

N_GROUPS = 38                  # number of E1 bins in the saddle
Q_LO, Q_HI = 0.025, 0.975      # E1 quantile range of the saddle
STRENGTH_FRAC = 0.20           # corner fraction used for compartment strength


# --------------------------------------------------------------- datasets ---

LIBRARIES = [
    ("U54-ESC4DN-DSG-DpnII-20190530-R1-T1",     "SRR13601502"),
    ("U54-ESC4DN-DSG-DpnII-20190530-R2-T1",     "SRR13601511"),
    ("U54-ESC4DN-FA-DpnII-2017524-R1-T1",       "SRR13601520"),
    ("U54-ESC4DN-FA-DpnII-2017524-R1-T2",       "SRR13601529"),
    ("U54-HFFc6-DSG-DdeI-DpnII-20190711-R1-T1", "SRR13601574"),
    ("U54-HFFc6-DSG-DdeI-DpnII-20191219-R3-T1", "SRR13601583"),
    ("U54-HFFc6-DSG-DpnII-20180319-R1-T1",      "SRR13601592"),
    ("U54-HFFc6-DSG-DpnII-20190102-R2-T1",      "SRR13601599"),
    ("kbm7",                                    "SRR1658708"),
]

# saddle plots in figure 4
SADDLE_LIBRARY = ("U54-ESC4DN-DSG-DpnII-20190530-R1-T1", "SRR13601502")

BIO_GROUPS = {
        "ESC4DN-DSG-DpnII": [
        ("U54-ESC4DN-DSG-DpnII-20190530-R1-T1", "SRR13601502"),
        ("U54-ESC4DN-DSG-DpnII-20190530-R2-T1", "SRR13601511"),
    ],
    "HFFc6-DSG-DdeI-DpnII": [
        ("U54-HFFc6-DSG-DdeI-DpnII-20190711-R1-T1", "SRR13601574"),
        ("U54-HFFc6-DSG-DdeI-DpnII-20191219-R3-T1", "SRR13601583"),
    ],
    "HFFc6-DSG-DpnII": [
        ("U54-HFFc6-DSG-DpnII-20180319-R1-T1", "SRR13601592"),
        ("U54-HFFc6-DSG-DpnII-20190102-R2-T1", "SRR13601599"),
    ],
    "IMR90": [                                       
        ("IMR90", "SRR1658672"), ("IMR90", "SRR1658673"),
        ("IMR90", "SRR1658674"), ("IMR90", "SRR1658675"),
        ("IMR90", "SRR1658676"),  
        ("IMR90", "SRR1658677"), ("IMR90", "SRR1658678"),
        ("IMR90", "SRR1658679"),
    ],
    "HMEC": [                                        
        ("HMEC", "SRR1658680"),  
        ("HMEC", "SRR1658681"), ("HMEC", "SRR1658682"),
        ("HMEC", "SRR1658683"), ("HMEC", "SRR1658684"),
        ("HMEC", "SRR1658685"),
    ],
    "k562": [                                       
        ("k562", "SRR1658693"), ("k562", "SRR1658694"),
        ("k562", "SRR1658695"),    
        ("k562", "SRR1658697"),    
        ("k562", "SRR1658699"),    
        ("k562", "SRR1658701"),    
    ],
    "kbm7": [                                       
        ("kbm7", "SRR1658703"),
        ("kbm7", "SRR1658705"),    
        ("kbm7", "SRR1658706"), ("kbm7", "SRR1658707"),
        ("kbm7", "SRR1658708"),    
    ],
    "HUVEC": [                                     
        ("HUVEC", "SRR1658709"),   
        ("HUVEC", "SRR1658712"),
        ("HUVEC", "SRR1658714"),  
        ("HUVEC", "SRR1658715"),
    ]}

TECH_GROUPS = {
    "U54-ESC4DN-DSG-DpnII-20190530-R1-T1": [
        "SRR13601502", "SRR13601503", "SRR13601504", "SRR13601505", "SRR13601506",
        "SRR13601507", "SRR13601508", "SRR13601509", "SRR13601510",
    ],
    "U54-ESC4DN-DSG-DpnII-20190530-R2-T1": [
        "SRR13601511", "SRR13601512", "SRR13601513", "SRR13601514", "SRR13601515",
        "SRR13601516", "SRR13601517", "SRR13601518", "SRR13601519",
    ],
    "U54-HFFc6-DSG-DdeI-DpnII-20190711-R1-T1": [
        "SRR13601574", "SRR13601575", "SRR13601576", "SRR13601577", "SRR13601578",
        "SRR13601579", "SRR13601580", "SRR13601581", "SRR13601582",
    ],
    "U54-HFFc6-DSG-DdeI-DpnII-20191219-R3-T1": [
        "SRR13601583", "SRR13601584", "SRR13601585", "SRR13601586", "SRR13601587",
        "SRR13601588", "SRR13601589", "SRR13601590", "SRR14276401",
    ],
    "U54-HFFc6-DSG-DpnII-20180319-R1-T1": [
        "SRR13601592", "SRR13601593", "SRR13601594", "SRR13601595", "SRR13601596",
        "SRR13601597", "SRR13601598",
    ],
    "U54-HFFc6-DSG-DpnII-20190102-R2-T1": [
        "SRR13601599", "SRR13601600", "SRR13601601", "SRR13601602", "SRR13601603",
        "SRR13601604", "SRR13601605", "SRR13601606",
    ],
}


def all_libraries():
    '''Every (sample_name, SRR) that has to be computed, without duplicates.'''
    pairs = list(LIBRARIES)
    for group in BIO_GROUPS.values():
        pairs += group
    for sample_name, srrs in TECH_GROUPS.items():
        pairs += [(sample_name, srr) for srr in srrs]

    unique = []
    for pair in pairs:
        if pair not in unique:
            unique.append(pair)
    return unique


# ------------------------------------------------------------------ files ---

def full_mcool(sample_name, srr):
    '''Path to the full-depth mcool of one library.'''
    return str(RAW_ROOT / sample_name / "coolers_library" / f"{srr}.{SUFFIX}.mcool")


def list_entities(sample_name, srr):
    '''All coolers of one library: the full one plus every subsample.

    Every entity is a dict with entity_id, role, depth, sample_idx and path.
    '''
    entities = [dict(entity_id=f"{srr}_full", role="full",
                     depth=FULL_DEPTH_READS, sample_idx=0,
                     path=full_mcool(sample_name, srr))]

    pattern = str(RAW_ROOT / sample_name / srr / f"{srr}_*_sample*.{SUFFIX}.mcool")
    for path in sorted(glob.glob(pattern)):
        match = re.match(rf"^{srr}_(\d+)_sample(\d+)\.", os.path.basename(path))
        if match is None:
            continue
        depth = int(match.group(1))
        index = int(match.group(2))
        entities.append(dict(entity_id=f"{srr}_{depth}_sample{index}",
                             role="subsample", depth=depth,
                             sample_idx=index, path=path))
    return entities


def eigs_path(sample_name, entity_id):
    '''Where the E1 track of one cooler is cached.'''
    return COMPARTMENT_CACHE / sample_name / f"{entity_id}.eigs.parquet"


def saddle_path(sample_name, entity_id):
    '''Where the saddle of one cooler is cached.'''
    return COMPARTMENT_CACHE / sample_name / f"{entity_id}.saddle.npz"


def insulation_path(sample_name, entity_id):
    '''Where the insulation table of one cooler is cached.'''
    return INSULATION_CACHE / sample_name / f"{entity_id}.ins.parquet"


def is_cached(path):
    '''True if the file exists and is not empty.'''
    return os.path.exists(path) and os.path.getsize(path) > 0


def save_parquet(df, path):
    '''Write a table, first to a temporary file, so a crash cannot leave junk.'''
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = str(path) + ".tmp"
    df.to_parquet(tmp, index=False)
    os.replace(tmp, path)


def save_npz(path, **arrays):
    '''Write arrays the same careful way.'''
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = str(path) + ".tmp"
    np.savez_compressed(tmp, **arrays)
    os.replace(tmp + ".npz", path)


def load_eigs(sample_name, entity_id):
    '''Read a cached E1 track.'''
    return pd.read_parquet(eigs_path(sample_name, entity_id))


def load_insulation(sample_name, entity_id):
    '''Read a cached insulation table.'''
    return pd.read_parquet(insulation_path(sample_name, entity_id))


def load_saddle(sample_name, entity_id):
    '''Read a cached saddle: obs/exp matrix and compartment strength.'''
    with np.load(saddle_path(sample_name, entity_id)) as data:
        interaction_sum = data["interaction_sum"]
        interaction_count = data["interaction_count"]
        strength = float(data["comp_strength"][0])
    with np.errstate(invalid="ignore", divide="ignore"):
        matrix = interaction_sum / interaction_count
    return dict(matrix=matrix, strength=strength)


# ------------------------------------------------------------- computation ---

def open_clr(path, resolution):
    '''Open one resolution of an mcool file.'''
    return cooler.Cooler(f"{path}::resolutions/{resolution}")


def make_view(clr):
    '''Whole-chromosome view frame with autosomes and chrX only.'''
    main_chrom = re.compile(r"^chr([0-9]+|X)$")
    chroms = [c for c in clr.chromnames if main_chrom.match(c)]
    if not chroms:
        chroms = list(clr.chromnames)
    table = pd.DataFrame({"chrom": chroms,
                          "start": 0,
                          "end": [int(clr.chromsizes[c]) for c in chroms],
                          "name": chroms})
    return bioframe.make_viewframe(table, check_bounds=clr.chromsizes)


def ensure_weights(path, resolution):
    '''Balance the cooler if it has no weights yet. True if weights are there.'''
    clr = open_clr(path, resolution)
    if "weight" in clr.bins().columns:
        return True
    bias, stats = cooler.balance_cooler(clr, cis_only=False, store=True,
                                        ignore_diags=2, max_iters=200)
    return bool(stats.get("converged", True))


def load_phasing_track(clr):
    '''GC content on the bins of the cooler; it orients the sign of E1.'''
    gc_track = pd.read_csv(GC_COV_FILE, sep="\t")
    bins = clr.bins()[:][["chrom", "start", "end"]]
    track = bins.merge(gc_track, on=["chrom", "start", "end"], how="left")
    value_column = [c for c in track.columns
                    if c not in ("chrom", "start", "end")][0]
    return track[["chrom", "start", "end", value_column]]


def compute_eigs(clr, phasing, view):
    '''E1 eigenvector track of one cooler.'''
    eigvals, eigvecs = cooltools.eigs_cis(clr, phasing, view_df=view, n_eigs=3)
    return eigvecs[["chrom", "start", "end", "E1"]].copy()


def compute_saddle(clr, eig_track, view, nproc):
    '''Saddle of one cooler plus the compartment strength.

    Strength is (AA + BB) / (AB + BA) taken over the corner groups of the
    observed/expected saddle.
    '''
    expected = cooltools.expected_cis(clr, view_df=view, nproc=nproc)
    interaction_sum, interaction_count = cooltools.saddle(
        clr, expected, eig_track, "cis",
        n_bins=N_GROUPS, qrange=(Q_LO, Q_HI), view_df=view,
    )

    with np.errstate(invalid="ignore", divide="ignore"):
        matrix = interaction_sum / interaction_count
    matrix = matrix[1:-1, 1:-1]           # drop the two outlier flanks

    corner = max(1, int(round(STRENGTH_FRAC * matrix.shape[0])))
    bb = np.nanmean(matrix[:corner, :corner])     # B with B
    aa = np.nanmean(matrix[-corner:, -corner:])   # A with A
    ab = np.nanmean(matrix[:corner, -corner:])
    ba = np.nanmean(matrix[-corner:, :corner])
    strength = (aa + bb) / (ab + ba) if (ab + ba) else np.nan

    return interaction_sum, interaction_count, strength


def compute_insulation(clr, nproc):
    '''Insulation table with Li-thresholded boundaries for both windows.'''
    view = make_view(clr)
    return cooltools.insulation(clr, INSUL_WINDOWS, view_df=view,
                                threshold="Li", nproc=nproc, verbose=False)


# ------------------------------------------------------------- comparisons ---

def aligned_values(test_df, ref_df, column):
    '''Values of one column in two tracks, matched bin by bin, NaNs removed.'''
    merged = test_df[["chrom", "start", "end", column]].merge(
        ref_df[["chrom", "start", "end", column]],
        on=["chrom", "start", "end"], suffixes=("_test", "_ref"))
    test = merged[f"{column}_test"].to_numpy()
    ref = merged[f"{column}_ref"].to_numpy()
    keep = np.isfinite(test) & np.isfinite(ref)
    return test[keep], ref[keep]


def track_correlations(test_df, ref_df, column):
    '''Pearson and Spearman correlation of a track with the reference track.'''
    test, ref = aligned_values(test_df, ref_df, column)
    if test.size < 3:
        return dict(pearson=np.nan, spearman=np.nan, n_bins=int(test.size))
    return dict(pearson=float(pearsonr(test, ref)[0]),
                spearman=float(spearmanr(test, ref)[0]),
                n_bins=int(test.size))


def sign_agreement(test_df, ref_df, column="E1"):
    '''Fraction of bins where both tracks put the bin in the same compartment.'''
    test, ref = aligned_values(test_df, ref_df, column)
    if test.size == 0:
        return np.nan
    return float(np.mean(np.sign(test) == np.sign(ref)))


def bins_by_chrom(positions, resolution):
    '''Boundary positions as sorted bin numbers, one array per chromosome.'''
    out = {}
    for chrom, sub in positions.groupby("chrom"):
        out[chrom] = np.sort((sub["start"].to_numpy() // resolution).astype(int))
    return out


def count_matched(first, second, tol_bins):
    '''How many boundaries of the first set have a partner in the second set.'''
    matched = 0
    for chrom, positions in first.items():
        other = second.get(chrom)
        if other is None or len(other) == 0:
            continue
        for position in positions:
            nearest = np.min(np.abs(other - position))
            if nearest <= tol_bins:
                matched += 1
    return matched


def match_boundary_sets(test_pos, ref_pos, resolution, tol_bins):
    '''Precision, recall and F1 between two boundary sets.

    A boundary counts as found if a boundary of the other set lies within
    tol_bins bins.
    '''
    test = bins_by_chrom(test_pos, resolution)
    ref = bins_by_chrom(ref_pos, resolution)
    n_test = sum(len(v) for v in test.values())
    n_ref = sum(len(v) for v in ref.values())

    precision = count_matched(test, ref, tol_bins) / n_test if n_test else np.nan
    recall = count_matched(ref, test, tol_bins) / n_ref if n_ref else np.nan

    if precision and recall and np.isfinite(precision) and np.isfinite(recall):
        f1 = 2 * precision * recall / (precision + recall)
    else:
        f1 = np.nan

    return dict(precision=precision, recall=recall, f1=f1,
                n_test=n_test, n_ref=n_ref)


# ----------------------------------------------------------------- workers ---

def worker_compartments(task):
    '''Compute and cache E1 and the saddle of one cooler.'''
    sample_name = task["sample_name"]
    entity = task["entity"]
    eig_file = eigs_path(sample_name, entity["entity_id"])
    saddle_file = saddle_path(sample_name, entity["entity_id"])

    if is_cached(eig_file) and is_cached(saddle_file):
        return ("skip", entity["entity_id"])
    if not os.path.exists(entity["path"]):
        return ("missing", entity["entity_id"])

    try:
        if not ensure_weights(entity["path"], RES_COMPARTMENTS):
            return ("no weights", entity["entity_id"])

        clr = open_clr(entity["path"], RES_COMPARTMENTS)
        view = make_view(clr)
        eig = compute_eigs(clr, load_phasing_track(clr), view)
        save_parquet(eig, eig_file)

        inter_sum, inter_count, strength = compute_saddle(clr, eig, view,
                                                          task["nproc"])
        save_npz(saddle_file, interaction_sum=inter_sum,
                 interaction_count=inter_count,
                 comp_strength=np.array([strength]))
        return ("done", entity["entity_id"])
    except Exception as error:
        return ("failed: " + repr(error)[:70], entity["entity_id"])
    finally:
        gc.collect()


def worker_insulation(task):
    '''Compute and cache the insulation table of one cooler.'''
    sample_name = task["sample_name"]
    entity = task["entity"]
    out_file = insulation_path(sample_name, entity["entity_id"])

    if is_cached(out_file):
        return ("skip", entity["entity_id"])
    if not os.path.exists(entity["path"]):
        return ("missing", entity["entity_id"])

    try:
        if not ensure_weights(entity["path"], RES_INSULATION):
            return ("no weights", entity["entity_id"])

        clr = open_clr(entity["path"], RES_INSULATION)
        save_parquet(compute_insulation(clr, task["nproc"]), out_file)
        return ("done", entity["entity_id"])
    except Exception as error:
        return ("failed: " + repr(error)[:70], entity["entity_id"])
    finally:
        gc.collect()


def build_tasks(nproc):
    '''One task per cooler of every library; nproc is used inside a worker.'''
    tasks = []
    for sample_name, srr in all_libraries():
        for entity in list_entities(sample_name, srr):
            tasks.append(dict(sample_name=sample_name, entity=entity, nproc=nproc))
    return tasks


def run_parallel(tasks, worker, n_workers, label=""):
    '''Run a worker over all tasks in separate processes.

    Every task gets a fresh process, so memory is released between coolers.
    A task that fails returns a status string instead of raising, and cached
    coolers are skipped, so the notebook can simply be rerun.
    '''
    results = []
    total = len(tasks)

    if n_workers <= 1:
        for number, task in enumerate(tasks, 1):
            result = worker(task)
            results.append(result)
            print(f"[{number}/{total}] {label} {result[0]:12s} {result[1]}")
    else:
        context = multiprocessing.get_context("fork")
        with context.Pool(processes=n_workers, maxtasksperchild=1) as pool:
            for number, result in enumerate(pool.imap_unordered(worker, tasks), 1):
                results.append(result)
                print(f"[{number}/{total}] {label} {result[0]:12s} {result[1]}")

    print("summary:", dict(Counter(r[0].split(":")[0] for r in results)))
    return results
