//! ckl - CENP-B inter-box gap histograms and histogram divergences. Rendering of centromeres as CENP-B box spacing distributions.
//!
//! `ckl hist`  : box BED (rust_motif_scan --cenpb) -> per-region gap histograms, GCP Model 1 convention:
//!               boxes sorted by start, gap = |start_i - end_{i-1}| (GiuntaLab/GCP-Centeny R/model1.R @3c688ae3).
//!               Three output streams (Rule 5): summary (every requested region), outliers, exceptions.
//! `ckl score` : divergences between histograms. `pairwise` = every query x ref pair;
//!               `loio` = every query vs every chromosome panel built without the query's group
//!               (per-region normalize -> mean within individual -> mean over individuals).

use clap::{Parser, Subcommand};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

mod perturb;
mod xspecies;
mod kin;
mod refdep;

type Res<T> = Result<T, Box<dyn Error>>;

#[derive(Parser)]
#[command(name = "ckl", about = "CENP-B box spacing histograms of centromeres: build, compare, identify, perturb")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Build gap histograms per region.
    Hist {
        /// Box BED(s): chrom start end name score strand ... ('#' lines skipped).
        #[arg(long, required = true)]
        boxes: Vec<String>,
        /// Regions TSV: seq start end region_id (0-based half-open, sequence coordinates).
        /// Several rows may share a region_id (union; gaps never cross intervals). Absent = whole sequences.
        #[arg(long)]
        regions: Option<String>,
        /// Names of every scanned sequence (one per line or .fai). With it, a region on an unknown
        /// sequence is a hard error and a known sequence without boxes is a genuine zero.
        #[arg(long)]
        seqs: Option<String>,
        /// Box count below which a region is `insufficient`.
        #[arg(long, default_value_t = 50)]
        min_boxes: u64,
        /// Gap (bp) above which a within-region gap is written to the exceptions stream.
        #[arg(long, default_value_t = 20000)]
        long_gap: i64,
        /// Output prefix: <out>.hist.tsv, .summary.tsv, .outliers.tsv, .exceptions.tsv
        #[arg(long)]
        out: String,
    },
    /// Divergences between histograms.
    Score {
        #[arg(long)]
        hist: String,
        #[arg(long)]
        summary: String,
        /// TSV with header: region_id individual group hap chrom role  (role: query|ref|panel|both)
        #[arg(long)]
        meta: String,
        /// pairwise | loio
        #[arg(long, default_value = "loio")]
        mode: String,
        /// Epsilon grid for smoothed KL.
        #[arg(long, default_value = "1e-10,1e-6,1e-3")]
        eps: String,
        /// loio only: one row per query (own-chromosome metrics + nearest chromosome by JS and by the first-eps KL(q||panel))
        #[arg(long, default_value_t = false)]
        own_only: bool,
        #[arg(long)]
        out: String,
    },
    /// Error-class battery: edit sequences at HOR-copy boundaries, re-scan, histogram, score vs source.
    Perturb {
        /// FASTA of source centromeres (headers = region ids)
        #[arg(long)]
        sources: String,
        /// TSV seq start end: HOR copies in sequence coordinates
        #[arg(long)]
        units: String,
        /// FASTA of donor active-HOR sequences (other individuals)
        #[arg(long)]
        donors: String,
        /// TSV donor_header chrom
        #[arg(long)]
        donor_meta: String,
        /// TSV chrom family (suprachromosomal family of the dominant active HOR)
        #[arg(long)]
        sf: String,
        #[arg(long, default_value = "1,5,20,100,500")]
        doses: String,
        #[arg(long, default_value_t = 10)]
        placements: u32,
        #[arg(long, default_value_t = 20260928)]
        seed: u64,
        #[arg(long)]
        out: String,
    },
    /// CenSat active_hor BED -> regions / meta / family (labels from the contig suffix hsaN)
    CensatRegions {
        #[arg(long)]
        censat: String,
        #[arg(long, required = true)]
        fai: Vec<String>,
        #[arg(long)]
        species: String,
        #[arg(long)]
        out: String,
    },
    /// Region sequences -> one FASTA per region + ids table (path region_id)
    Extract {
        #[arg(long, required = true)]
        fasta: Vec<String>,
        #[arg(long)]
        regions: String,
        #[arg(long)]
        out_dir: String,
        #[arg(long)]
        ids_out: String,
    },
    /// Rewrite a meta table with one group name; family per chromosome from a map (chrom family)
    MetaSet {
        #[arg(long)]
        meta: String,
        #[arg(long)]
        group: String,
        #[arg(long)]
        family_map: String,
        /// override the role column (e.g. ref)
        #[arg(long)]
        role: Option<String>,
        #[arg(long)]
        out: String,
    },
    /// 1-NN identification for every ordered pair of units (group or individual) from a distance table or size only
    Xident {
        /// distance table: with --metric, a header TSV whose first two columns are query and ref (directional);
        /// without, headerless `a b d` rows treated as symmetric (Mash)
        #[arg(long)]
        dist: Option<String>,
        #[arg(long)]
        metric: Option<String>,
        /// id -> region_id map for the distance table (e.g. Mash file paths)
        #[arg(long)]
        ids: Option<String>,
        /// size-only null: ckl summary files (region_bp, n_boxes); replaces --dist
        #[arg(long)]
        size: Vec<String>,
        #[arg(long, required = true)]
        meta: Vec<String>,
        #[arg(long)]
        family: Vec<String>,
        /// group | individual
        #[arg(long, default_value = "group")]
        by: String,
        /// only regions of these groups are queries (comma list)
        #[arg(long)]
        query_groups: Option<String>,
        /// restrict queries to region ids in this file (first column, or column 2 of a .detail.tsv)
        #[arg(long)]
        queries_from: Option<String>,
        /// labels never scored (comma list)
        #[arg(long, default_value = "")]
        exclude: String,
        #[arg(long, default_value_t = 1000)]
        boot: usize,
        #[arg(long, default_value_t = 20260928)]
        seed: u64,
        #[arg(long)]
        label: String,
        #[arg(long)]
        out: String,
    },
    /// Reference dependence: pool component histograms into consensus recipes and score query assemblies
    Refdep {
        /// TAG=HIST_PREFIX[@REGION_TO_CHROM_MAP], repeatable
        #[arg(long = "component", required = true)]
        components: Vec<String>,
        /// NAME=TAG,TAG,..., repeatable
        #[arg(long = "recipe", required = true)]
        recipes: Vec<String>,
        /// component TAG to score, repeatable
        #[arg(long = "query", required = true)]
        queries: Vec<String>,
        #[arg(long, default_value_t = 1e-12)]
        eps: f64,
        #[arg(long)]
        out: String,
    },
    /// Within-chromosome nearest neighbours: trio identity-by-descent and population-group tests
    Kin {
        #[arg(long)]
        dist: String,
        #[arg(long)]
        metric: Option<String>,
        #[arg(long)]
        ids: Option<String>,
        #[arg(long)]
        meta: String,
        /// TSV child father mother
        #[arg(long)]
        trios: Option<String>,
        /// TSV individual label
        #[arg(long)]
        labels: Option<String>,
        #[arg(long, default_value = "chrY")]
        exclude: String,
        #[arg(long)]
        label: String,
        #[arg(long)]
        out: String,
    },
    /// Held-out nearest neighbour from an all-vs-all distance table (e.g. `mash dist` output):
    /// per query and chromosome, the minimum distance to a reference of another group.
    Nn {
        /// TSV: id_a id_b distance [...] (no header), ids mapped through --ids
        #[arg(long)]
        dist: String,
        /// TSV: id region_id (maps the distance-table ids, e.g. file paths, to meta region ids)
        #[arg(long)]
        ids: String,
        #[arg(long)]
        meta: String,
        #[arg(long)]
        out: String,
    },
}

fn main() {
    let cli = Cli::parse();
    let r = match cli.cmd {
        Cmd::Hist { boxes, regions, seqs, min_boxes, long_gap, out } => {
            run_hist(&boxes, regions.as_deref(), seqs.as_deref(), min_boxes, long_gap, &out)
        }
        Cmd::Score { hist, summary, meta, mode, eps, own_only, out } => run_score(&hist, &summary, &meta, &mode, &eps, own_only, &out),
        Cmd::Perturb { sources, units, donors, donor_meta, sf, doses, placements, seed, out } => {
            perturb::run(&perturb::Args { sources, units, donors, donor_meta, sf, doses, placements, seed, out })
        }
        Cmd::Nn { dist, ids, meta, out } => run_nn(&dist, &ids, &meta, &out),
        Cmd::Refdep { components, recipes, queries, eps, out } => refdep::run(&refdep::RefdepArgs { components, recipes, queries, eps, out }),
        Cmd::Kin { dist, metric, ids, meta, trios, labels, exclude, label, out } => {
            kin::run(&kin::KinArgs { dist, metric, ids, meta, trios, labels, exclude, label, out })
        }
        Cmd::CensatRegions { censat, fai, species, out } => xspecies::run_censat_regions(&censat, &fai, &species, &out),
        Cmd::Extract { fasta, regions, out_dir, ids_out } => xspecies::run_extract(&fasta, &regions, &out_dir, &ids_out),
        Cmd::MetaSet { meta, group, family_map, role, out } => xspecies::run_meta_set(&meta, &group, &family_map, role.as_deref(), &out),
        Cmd::Xident { dist, metric, ids, size, meta, family, by, query_groups, queries_from, exclude, boot, seed, label, out } => {
            xspecies::run_xident(&xspecies::XidentArgs { dist, metric, ids, size, meta, family, by, query_groups, queries_from, exclude, boot, seed, label, out })
        }
    };
    if let Err(e) = r {
        eprintln!("ckl: error: {e}");
        std::process::exit(1);
    }
}

// ───────────────────────────── hist ─────────────────────────────

type Boxes = HashMap<String, Vec<(i64, i64, bool)>>;

fn open(p: &str) -> Res<BufReader<File>> {
    Ok(BufReader::new(File::open(p).map_err(|e| format!("open {p}: {e}"))?))
}

fn create(p: &str) -> Res<BufWriter<File>> {
    Ok(BufWriter::new(File::create(p).map_err(|e| format!("create {p}: {e}"))?))
}

fn read_boxes(paths: &[String]) -> Res<Boxes> {
    let mut m: Boxes = HashMap::new();
    for p in paths {
        for (ln, line) in open(p)?.lines().enumerate() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let t: Vec<&str> = line.split('\t').collect();
            if t.len() < 6 {
                return Err(format!("{p}:{}: expected >= 6 BED columns, got {}", ln + 1, t.len()).into());
            }
            let s: i64 = t[1].parse().map_err(|e| format!("{p}:{}: start {:?}: {e}", ln + 1, t[1]))?;
            let e: i64 = t[2].parse().map_err(|e| format!("{p}:{}: end {:?}: {e}", ln + 1, t[2]))?;
            if e <= s {
                return Err(format!("{p}:{}: end <= start", ln + 1).into());
            }
            let fwd = match t[5] {
                "+" => true,
                "-" => false,
                x => return Err(format!("{p}:{}: strand {x:?}", ln + 1).into()),
            };
            m.entry(t[0].to_string()).or_default().push((s, e, fwd));
        }
    }
    for v in m.values_mut() {
        v.sort_unstable();
    }
    Ok(m)
}

fn read_seq_names(p: &str) -> Res<HashSet<String>> {
    let mut s = HashSet::new();
    for line in open(p)?.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let name = line.split('\t').next().unwrap_or("").trim_start_matches('>').to_string();
        if !name.is_empty() {
            s.insert(name);
        }
    }
    Ok(s)
}

struct Interval {
    rid: String,
    seq: String,
    start: i64,
    end: i64,
}

fn read_regions(p: &str) -> Res<Vec<Interval>> {
    let mut v = Vec::new();
    for (ln, line) in open(p)?.lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() < 4 {
            return Err(format!("{p}:{}: expected seq start end region_id", ln + 1).into());
        }
        let start: i64 = t[1].parse().map_err(|e| format!("{p}:{}: start: {e}", ln + 1))?;
        let end: i64 = t[2].parse().map_err(|e| format!("{p}:{}: end: {e}", ln + 1))?;
        if end <= start {
            return Err(format!("{p}:{}: end <= start", ln + 1).into());
        }
        v.push(Interval { rid: t[3].to_string(), seq: t[0].to_string(), start, end });
    }
    Ok(v)
}

#[derive(Default)]
struct Acc {
    n_int: u64,
    bp: i64,
    n_boxes: u64,
    n_fwd: u64,
    n_rev: u64,
    n_gaps: u64,
    n_overlap: u64,
    max_gap: i64,
    hist: BTreeMap<i64, u64>,
    no_box_seq: bool,
}

/// Gaps between consecutive boxes fully inside [start, end); GCP Model 1 convention.
/// Returns exceptions as (position, gap, reason).
fn accumulate(acc: &mut Acc, boxes: &[(i64, i64, bool)], start: i64, end: i64, long_gap: i64) -> Vec<(i64, i64, &'static str)> {
    let mut exc = Vec::new();
    let first = boxes.partition_point(|b| b.0 < start);
    let mut prev_end: Option<i64> = None;
    let (mut lo, mut hi) = (i64::MAX, i64::MIN);
    for b in &boxes[first..] {
        if b.0 >= end {
            break;
        }
        if b.1 > end {
            continue; // straddles the region end
        }
        acc.n_boxes += 1;
        if b.2 { acc.n_fwd += 1 } else { acc.n_rev += 1 }
        lo = lo.min(b.0);
        hi = hi.max(b.1);
        if let Some(pe) = prev_end {
            let d = b.0 - pe;
            if d < 0 {
                acc.n_overlap += 1;
                exc.push((b.0, d, "overlap_abs_applied"));
            }
            let g = d.abs();
            *acc.hist.entry(g).or_insert(0) += 1;
            acc.n_gaps += 1;
            acc.max_gap = acc.max_gap.max(g);
            if g > long_gap {
                exc.push((b.0, g, "long_gap_inside_region"));
            }
        }
        prev_end = Some(b.1);
    }
    acc.n_int += 1;
    acc.bp += if end == i64::MAX { if hi > lo { hi - lo } else { 0 } } else { end - start };
    exc
}

fn run_hist(boxes_p: &[String], regions_p: Option<&str>, seqs_p: Option<&str>, min_boxes: u64, long_gap: i64, out: &str) -> Res<()> {
    let boxes = read_boxes(boxes_p)?;
    let known = match seqs_p {
        Some(p) => Some(read_seq_names(p)?),
        None => None,
    };
    let intervals: Vec<Interval> = match regions_p {
        Some(p) => read_regions(p)?,
        None => {
            let mut names: Vec<&String> = boxes.keys().collect();
            names.sort();
            names.into_iter().map(|s| Interval { rid: s.clone(), seq: s.clone(), start: 0, end: i64::MAX }).collect()
        }
    };
    let mut order: Vec<String> = Vec::new();
    let mut accs: HashMap<String, Acc> = HashMap::new();
    let mut exc_w = create(&format!("{out}.exceptions.tsv"))?;
    writeln!(exc_w, "region_id\tseq\tpos\tgap\treason\tdefault_applied")?;
    let empty: Vec<(i64, i64, bool)> = Vec::new();
    for iv in &intervals {
        if !accs.contains_key(&iv.rid) {
            order.push(iv.rid.clone());
        }
        let acc = accs.entry(iv.rid.clone()).or_default();
        let bx = match boxes.get(&iv.seq) {
            Some(b) => b,
            None => {
                if let Some(k) = &known {
                    if !k.contains(&iv.seq) {
                        return Err(format!("region {} names sequence {:?} that was not scanned", iv.rid, iv.seq).into());
                    }
                }
                acc.no_box_seq = true;
                &empty
            }
        };
        for (pos, gap, why) in accumulate(acc, bx, iv.start, iv.end, long_gap) {
            let applied = if why == "overlap_abs_applied" { "abs(gap) counted (GCP behaviour)" } else { "counted" };
            writeln!(exc_w, "{}\t{}\t{}\t{}\t{}\t{}", iv.rid, iv.seq, pos, gap, why, applied)?;
        }
    }
    let mut hw = create(&format!("{out}.hist.tsv"))?;
    let mut sw = create(&format!("{out}.summary.tsv"))?;
    let mut ow = create(&format!("{out}.outliers.tsv"))?;
    writeln!(hw, "region_id\tgap\tcount")?;
    writeln!(sw, "region_id\tn_intervals\tregion_bp\tn_boxes\tn_fwd\tn_rev\tn_gaps\tn_overlap\tmax_gap\tstatus")?;
    writeln!(ow, "region_id\treason")?;
    for rid in &order {
        let a = &accs[rid];
        for (g, c) in &a.hist {
            writeln!(hw, "{rid}\t{g}\t{c}")?;
        }
        let status = if a.n_boxes >= min_boxes { "ok" } else { "insufficient" };
        writeln!(sw, "{rid}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{status}", a.n_int, a.bp, a.n_boxes, a.n_fwd, a.n_rev, a.n_gaps, a.n_overlap, a.max_gap)?;
        if a.n_gaps == 0 {
            let why = if a.no_box_seq && a.n_boxes == 0 {
                if known.is_some() { "zero_boxes_in_sequence" } else { "no_boxes_or_name_mismatch(no --seqs given)" }
            } else if a.n_boxes == 0 {
                "zero_boxes_in_region"
            } else {
                "one_box_no_gap"
            };
            writeln!(ow, "{rid}\t{why}")?;
        }
    }
    eprintln!("ckl hist: {} regions, {} sequences with boxes", order.len(), boxes.len());
    Ok(())
}

// ───────────────────────────── score ─────────────────────────────

/// Sparse distribution sorted by gap; probabilities sum to 1; `n` = raw gap count (0 for panels).
#[derive(Clone)]
struct Dist {
    bins: Vec<(i64, f64)>,
    n: f64,
    counts: bool,
}

fn from_counts(c: &[(i64, f64)]) -> Dist {
    let n: f64 = c.iter().map(|x| x.1).sum();
    Dist { bins: c.iter().map(|&(g, v)| (g, v / n)).collect(), n, counts: true }
}

fn mean_of(ds: &[&Dist]) -> Dist {
    let mut m: BTreeMap<i64, f64> = BTreeMap::new();
    let k = ds.len() as f64;
    for d in ds {
        for &(g, p) in &d.bins {
            *m.entry(g).or_insert(0.0) += p / k;
        }
    }
    Dist { bins: m.into_iter().collect(), n: 0.0, counts: false }
}

/// Union support: (gap, p, q).
fn merge(p: &Dist, q: &Dist) -> Vec<(i64, f64, f64)> {
    let (mut i, mut j) = (0, 0);
    let mut v = Vec::with_capacity(p.bins.len() + q.bins.len());
    while i < p.bins.len() || j < q.bins.len() {
        let gp = p.bins.get(i).map(|x| x.0).unwrap_or(i64::MAX);
        let gq = q.bins.get(j).map(|x| x.0).unwrap_or(i64::MAX);
        if gp == gq {
            v.push((gp, p.bins[i].1, q.bins[j].1));
            i += 1;
            j += 1;
        } else if gp < gq {
            v.push((gp, p.bins[i].1, 0.0));
            i += 1;
        } else {
            v.push((gq, 0.0, q.bins[j].1));
            j += 1;
        }
    }
    v
}

/// KL(a||b) after adding eps to every bin of a support of size ns and renormalizing.
/// Bins outside the union carry equal mass in a and b and contribute 0.
fn kl_renorm(u: &[(i64, f64, f64)], eps: f64, ns: f64, a_is_p: bool) -> f64 {
    let z = 1.0 + eps * ns;
    u.iter()
        .map(|&(_, p, q)| {
            let (x, y) = if a_is_p { (p, q) } else { (q, p) };
            let a = (x + eps) / z;
            let b = (y + eps) / z;
            a * (a / b).ln()
        })
        .sum()
}

/// Naive KL: sum over a>0 of a ln(a / (b + eps)), no renormalization.
fn kl_naive(u: &[(i64, f64, f64)], eps: f64, a_is_p: bool) -> f64 {
    u.iter()
        .map(|&(_, p, q)| {
            let (x, y) = if a_is_p { (p, q) } else { (q, p) };
            if x > 0.0 { x * (x / (y + eps)).ln() } else { 0.0 }
        })
        .sum()
}

fn metrics(p: &Dist, q: &Dist, eps: &[f64]) -> Vec<f64> {
    let u = merge(p, q);
    let n_union = u.len() as f64;
    let max_gap = u.last().map(|x| x.0).unwrap_or(0);
    let n_dense = (max_gap + 1) as f64;
    let mut out = Vec::new();
    for &e in eps {
        out.push(kl_renorm(&u, e, n_union, true));
        out.push(kl_renorm(&u, e, n_union, false));
        out.push(kl_renorm(&u, e, n_dense, true));
        out.push(kl_renorm(&u, e, n_dense, false));
        out.push(kl_naive(&u, e, true));
        out.push(kl_naive(&u, e, false));
    }
    // Laplace (+1 count) on the union support: only defined when both sides carry counts.
    if p.counts && q.counts {
        let (np, nq) = (p.n, q.n);
        let lap = |a_is_p: bool| -> f64 {
            u.iter()
                .map(|&(_, pp, qq)| {
                    let a = (pp * np + 1.0) / (np + n_union);
                    let b = (qq * nq + 1.0) / (nq + n_union);
                    if a_is_p { a * (a / b).ln() } else { b * (b / a).ln() }
                })
                .sum()
        };
        out.push(lap(true));
        out.push(lap(false));
    } else {
        out.push(f64::NAN);
        out.push(f64::NAN);
    }
    let ln2 = std::f64::consts::LN_2;
    let (mut js, mut bc, mut l2, mut mn, mut mx, mut w1) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut fp, mut fq) = (0.0, 0.0);
    for k in 0..u.len() {
        let (g, pp, qq) = u[k];
        let m = 0.5 * (pp + qq);
        if pp > 0.0 { js += 0.5 * pp * (pp / m).ln(); }
        if qq > 0.0 { js += 0.5 * qq * (qq / m).ln(); }
        bc += (pp * qq).sqrt();
        l2 += (pp - qq) * (pp - qq);
        mn += pp.min(qq);
        mx += pp.max(qq);
        fp += pp;
        fq += qq;
        if k + 1 < u.len() {
            w1 += (fp - fq).abs() * (u[k + 1].0 - g) as f64;
        }
    }
    out.push((js / ln2).max(0.0).sqrt()); // JS distance, base 2, in [0,1]
    out.push((1.0 - bc).max(0.0).sqrt()); // Hellinger distance
    out.push(l2.sqrt()); // Euclidean
    out.push(if mx > 0.0 { 1.0 - mn / mx } else { f64::NAN }); // weighted Jaccard distance
    out.push(w1); // 1-D Wasserstein, bp
    out
}

fn metric_names(eps: &[f64]) -> Vec<String> {
    let mut v = Vec::new();
    for e in eps {
        let t = format!("{e:e}");
        for s in ["renorm_union", "renorm_dense", "naive"] {
            v.push(format!("kl_{s}_{t}_pq"));
            v.push(format!("kl_{s}_{t}_qp"));
        }
    }
    v.extend(["kl_laplace_union_pq", "kl_laplace_union_qp", "js_dist", "hellinger", "euclid", "wjaccard_dist", "w1_bp"].map(String::from));
    v
}

struct Meta {
    individual: String,
    group: String,
    hap: String,
    chrom: String,
    role: String,
}

fn read_meta(p: &str) -> Res<Vec<(String, Meta)>> {
    let mut lines = open(p)?.lines();
    let header = lines.next().ok_or("meta: empty file")??;
    let cols: Vec<&str> = header.split('\t').collect();
    let idx = |name: &str| -> Res<usize> { cols.iter().position(|c| *c == name).ok_or_else(|| format!("meta: missing column {name}").into()) };
    let (ir, ii, ig, ih, ic, io) = (idx("region_id")?, idx("individual")?, idx("group")?, idx("hap")?, idx("chrom")?, idx("role")?);
    let mut v = Vec::new();
    for (ln, line) in lines.enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() < cols.len() {
            return Err(format!("meta:{}: {} columns, header has {}", ln + 2, t.len(), cols.len()).into());
        }
        v.push((t[ir].to_string(), Meta { individual: t[ii].into(), group: t[ig].into(), hap: t[ih].into(), chrom: t[ic].into(), role: t[io].into() }));
    }
    Ok(v)
}

fn read_hist(p: &str) -> Res<HashMap<String, Vec<(i64, f64)>>> {
    let mut m: HashMap<String, Vec<(i64, f64)>> = HashMap::new();
    for (ln, line) in open(p)?.lines().enumerate().skip(1) {
        let line = line?;
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() != 3 {
            return Err(format!("{p}:{}: expected region_id gap count", ln + 1).into());
        }
        let g: i64 = t[1].parse()?;
        let c: f64 = t[2].parse()?;
        m.entry(t[0].to_string()).or_default().push((g, c));
    }
    for v in m.values_mut() {
        v.sort_by_key(|x| x.0);
    }
    Ok(m)
}

fn read_summary(p: &str) -> Res<HashMap<String, (String, u64)>> {
    let mut m = HashMap::new();
    for line in open(p)?.lines().skip(1) {
        let line = line?;
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() < 10 {
            return Err(format!("{p}: short summary row").into());
        }
        m.insert(t[0].to_string(), (t[9].to_string(), t[3].parse()?));
    }
    Ok(m)
}

#[derive(Default)]
struct ChromPanel {
    total: BTreeMap<i64, f64>,
    n: usize,
    by_group: HashMap<String, (BTreeMap<i64, f64>, usize)>,
}

impl ChromPanel {
    /// Equal-individual mean profile without `group`; None if nothing is left.
    fn without(&self, group: &str) -> Option<(Dist, usize)> {
        let (gs, k) = match self.by_group.get(group) {
            Some((m, k)) => (Some(m), *k),
            None => (None, 0),
        };
        let n = self.n - k;
        if n == 0 {
            return None;
        }
        let mut bins = Vec::with_capacity(self.total.len());
        for (&g, &t) in &self.total {
            let v = (t - gs.and_then(|m| m.get(&g)).copied().unwrap_or(0.0)) / n as f64;
            if v > 1e-12 {
                bins.push((g, v));
            }
        }
        let s: f64 = bins.iter().map(|x| x.1).sum();
        for b in bins.iter_mut() {
            b.1 /= s;
        }
        Some((Dist { bins, n: 0.0, counts: false }, n))
    }
}

fn fmt(x: f64) -> String {
    if x.is_nan() { "NA".into() } else { format!("{x}") }
}

fn run_score(hist_p: &str, sum_p: &str, meta_p: &str, mode: &str, eps_s: &str, own_only: bool, out: &str) -> Res<()> {
    let eps: Vec<f64> = eps_s.split(',').map(|s| s.trim().parse::<f64>()).collect::<Result<_, _>>()?;
    let hist = read_hist(hist_p)?;
    let summ = read_summary(sum_p)?;
    let meta = read_meta(meta_p)?;
    for (rid, _) in &meta {
        if !summ.contains_key(rid) {
            return Err(format!("meta region {rid:?} absent from summary (not histogrammed)").into());
        }
    }
    let dist = |rid: &str| -> Option<Dist> { hist.get(rid).map(|c| from_counts(c)) };
    let names = metric_names(&eps);
    let mut w = create(out)?;
    match mode {
        "pairwise" => {
            writeln!(w, "query\tref\t{}", names.join("\t"))?;
            let qs: Vec<&String> = meta.iter().filter(|(_, m)| m.role == "query" || m.role == "both").map(|(r, _)| r).collect();
            let rs: Vec<&String> = meta.iter().filter(|(_, m)| m.role == "ref" || m.role == "both").map(|(r, _)| r).collect();
            for q in &qs {
                for r in &rs {
                    let row = match (dist(q), dist(r)) {
                        (Some(a), Some(b)) => metrics(&a, &b, &eps).into_iter().map(fmt).collect::<Vec<_>>(),
                        _ => vec!["NA".to_string(); names.len()],
                    };
                    writeln!(w, "{q}\t{r}\t{}", row.join("\t"))?;
                }
            }
        }
        "loio" => {
            // chrom -> individual -> (group, member distributions)
            let mut members: BTreeMap<String, BTreeMap<String, (String, Vec<Dist>)>> = BTreeMap::new();
            for (rid, m) in &meta {
                if !(m.role == "panel" || m.role == "both") || summ[rid].0 != "ok" {
                    continue;
                }
                if let Some(d) = dist(rid) {
                    members.entry(m.chrom.clone()).or_default().entry(m.individual.clone()).or_insert_with(|| (m.group.clone(), Vec::new())).1.push(d);
                }
            }
            // Per chromosome: sum of individual profiles overall and per group; a held-out panel is
            // (total - group) / (n - k), i.e. the equal-individual mean without the query's group.
            let mut panels: BTreeMap<String, ChromPanel> = BTreeMap::new();
            for (c, inds) in &members {
                let cp = panels.entry(c.clone()).or_insert_with(ChromPanel::default);
                for (_ind, (grp, ds)) in inds {
                    let refs: Vec<&Dist> = ds.iter().collect();
                    let prof = mean_of(&refs);
                    cp.n += 1;
                    let gs = cp.by_group.entry(grp.clone()).or_insert_with(|| (BTreeMap::new(), 0));
                    gs.1 += 1;
                    for &(g, p) in &prof.bins {
                        *cp.total.entry(g).or_insert(0.0) += p;
                        *gs.0.entry(g).or_insert(0.0) += p;
                    }
                }
            }
            let mut cache: HashMap<(String, String), Option<(Dist, usize)>> = HashMap::new();
            let (i_kl, i_js) = (4usize, 6 * eps.len() + 2);
            if own_only {
                writeln!(w, "query\tq_individual\tq_hap\tq_chrom\tq_status\tq_n_boxes\tpanel_n_individuals\tnearest_js\tnearest_kl\t{}", names.join("\t"))?;
                for (rid, m) in &meta {
                    if !(m.role == "query" || m.role == "both") {
                        continue;
                    }
                    let (st, nb) = &summ[rid];
                    let qd = match dist(rid) {
                        Some(d) => d,
                        None => {
                            writeln!(w, "{rid}\t{}\t{}\t{}\t{st}\t{nb}\t0\tNA\tNA\t{}", m.individual, m.hap, m.chrom, vec!["NA"; names.len()].join("\t"))?;
                            continue;
                        }
                    };
                    let mut own: Option<(Vec<f64>, usize)> = None;
                    let (mut best_js, mut best_kl) = ((f64::INFINITY, String::from("NA")), (f64::INFINITY, String::from("NA")));
                    for (c, cp) in &panels {
                        let panel = cache.entry((m.group.clone(), c.clone())).or_insert_with(|| cp.without(&m.group));
                        if let Some((pd, k)) = panel {
                            let v = metrics(&qd, pd, &eps);
                            if v[i_js] < best_js.0 { best_js = (v[i_js], c.clone()); }
                            if v[i_kl] < best_kl.0 { best_kl = (v[i_kl], c.clone()); }
                            if *c == m.chrom { own = Some((v, *k)); }
                        }
                    }
                    let (row, k) = match own {
                        Some((v, k)) => (v.into_iter().map(fmt).collect::<Vec<_>>(), k),
                        None => (vec!["NA".to_string(); names.len()], 0),
                    };
                    writeln!(w, "{rid}\t{}\t{}\t{}\t{st}\t{nb}\t{k}\t{}\t{}\t{}", m.individual, m.hap, m.chrom, best_js.1, best_kl.1, row.join("\t"))?;
                }
                return Ok(());
            }
            writeln!(w, "query\tq_individual\tq_hap\tq_chrom\tq_status\tq_n_boxes\tpanel_chrom\tpanel_n_individuals\t{}", names.join("\t"))?;
            for (rid, m) in &meta {
                if !(m.role == "query" || m.role == "both") {
                    continue;
                }
                let (st, nb) = &summ[rid];
                let qd = dist(rid);
                for (c, cp) in &panels {
                    let key = (m.group.clone(), c.clone());
                    let panel = cache.entry(key).or_insert_with(|| cp.without(&m.group));
                    let (row, n_ind) = match (&qd, panel) {
                        (Some(a), Some((pd, k))) => (metrics(a, pd, &eps).into_iter().map(fmt).collect::<Vec<_>>(), *k),
                        (_, Some((_, k))) => (vec!["NA".to_string(); names.len()], *k),
                        _ => (vec!["NA".to_string(); names.len()], 0),
                    };
                    writeln!(w, "{rid}\t{}\t{}\t{}\t{st}\t{nb}\t{c}\t{n_ind}\t{}", m.individual, m.hap, m.chrom, row.join("\t"))?;
                }
            }
        }
        x => return Err(format!("unknown mode {x:?} (pairwise|loio)").into()),
    }
    Ok(())
}

// ───────────────────────────── nn ─────────────────────────────

fn run_nn(dist_p: &str, ids_p: &str, meta_p: &str, out: &str) -> Res<()> {
    let meta: HashMap<String, Meta> = read_meta(meta_p)?.into_iter().collect();
    let mut idmap: HashMap<String, String> = HashMap::new();
    for (ln, line) in open(ids_p)?.lines().enumerate() {
        let line = line?;
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() != 2 {
            return Err(format!("{ids_p}:{}: expected id region_id", ln + 1).into());
        }
        if !meta.contains_key(t[1]) {
            return Err(format!("{ids_p}:{}: region {:?} not in meta", ln + 1, t[1]).into());
        }
        idmap.insert(t[0].to_string(), t[1].to_string());
    }
    // query -> chrom -> min distance to a reference of another group
    let mut best: HashMap<String, BTreeMap<String, f64>> = HashMap::new();
    let mut n_rows: u64 = 0;
    for (ln, line) in open(dist_p)?.lines().enumerate() {
        let line = line?;
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() < 3 {
            return Err(format!("{dist_p}:{}: expected id_a id_b distance", ln + 1).into());
        }
        let a = idmap.get(t[0]).ok_or_else(|| format!("{dist_p}:{}: unknown id {:?}", ln + 1, t[0]))?;
        let b = idmap.get(t[1]).ok_or_else(|| format!("{dist_p}:{}: unknown id {:?}", ln + 1, t[1]))?;
        let d: f64 = t[2].parse().map_err(|e| format!("{dist_p}:{}: distance: {e}", ln + 1))?;
        n_rows += 1;
        for (q, r) in [(a, b), (b, a)] {
            let (mq, mr) = (&meta[q], &meta[r]);
            if mq.group == mr.group {
                continue;
            }
            let e = best.entry(q.clone()).or_default().entry(mr.chrom.clone()).or_insert(f64::INFINITY);
            if d < *e {
                *e = d;
            }
        }
    }
    let chroms: BTreeSet<String> = meta.values().map(|m| m.chrom.clone()).collect();
    let mut w = create(out)?;
    writeln!(w, "query\tq_chrom\tpredicted\tmin_dist\tsecond_chrom\tsecond_dist")?;
    let mut qs: Vec<&String> = best.keys().collect();
    qs.sort();
    for q in qs {
        let mut v: Vec<(f64, &String)> = best[q].iter().map(|(c, d)| (*d, c)).collect();
        v.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
        let (d1, c1) = v[0];
        let (d2, c2) = v.get(1).map(|x| (x.0, x.1.as_str())).unwrap_or((f64::NAN, "NA"));
        writeln!(w, "{q}\t{}\t{c1}\t{d1}\t{c2}\t{}", meta[q].chrom, fmt(d2))?;
    }
    eprintln!("ckl nn: {n_rows} distance rows, {} queries, {} chromosomes", best.len(), chroms.len());
    Ok(())
}

// ───────────────────────────── tests ─────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn hist_of(b: &[(i64, i64, bool)]) -> BTreeMap<i64, u64> {
        let mut a = Acc::default();
        let mut v = b.to_vec();
        v.sort_unstable();
        accumulate(&mut a, &v, 0, i64::MAX, 20000);
        a.hist
    }

    #[test]
    fn gap_is_edge_to_edge() {
        // boxes two monomers apart: gap = 342 - 17 = 325 (GCP Model 1)
        let h = hist_of(&[(0, 17, true), (342, 359, false)]);
        assert_eq!(h.get(&325), Some(&1));
    }

    #[test]
    fn region_excludes_straddling_box_and_never_crosses_intervals() {
        let v = vec![(0, 17, true), (171, 188, true), (342, 359, true), (1000, 1017, true)];
        let mut a = Acc::default();
        accumulate(&mut a, &v, 0, 350, 20000); // box at 342..359 straddles 350
        accumulate(&mut a, &v, 900, 2000, 20000);
        assert_eq!(a.n_boxes, 3);
        assert_eq!(a.n_gaps, 1); // only 0->171 inside the first interval; one box in the second
        assert_eq!(a.hist.get(&154), Some(&1));
    }

    #[test]
    fn composition_channel_is_blind_to_uniform_duplication() {
        let one = [(0, 17, true), (342, 359, true), (513, 530, true)];
        let two = [(0, 17, true), (342, 359, true), (513, 530, true), (855, 872, true), (1197, 1214, true), (1368, 1385, true)];
        let (p, q) = (hist_of(&one), hist_of(&two));
        let d = |h: &BTreeMap<i64, u64>| from_counts(&h.iter().map(|(g, c)| (*g, *c as f64)).collect::<Vec<_>>());
        let m = metrics(&d(&p), &d(&q), &[1e-6]);
        // one junction gap (855-530 = 325) joins an existing bin, so the doubled array keeps a nearly equal shape
        assert!(m[m.len() - 5] < 0.1, "JS distance {}", m[m.len() - 5]);
    }

    #[test]
    fn identical_distributions_score_zero_and_symmetric_kl_is_symmetric() {
        let a = from_counts(&[(154, 3.0), (325, 5.0), (496, 2.0)]);
        let b = from_counts(&[(154, 1.0), (325, 5.0), (667, 4.0)]);
        let aa = metrics(&a, &a, &[1e-6]);
        // naive KL on identical inputs is -eps-order, not exactly 0 (no renormalization)
        assert!(aa.iter().filter(|x| !x.is_nan()).all(|x| x.abs() < 1e-4), "{aa:?}");
        assert!(aa[0].abs() < 1e-12 && aa[2].abs() < 1e-12);
        let ab = metrics(&a, &b, &[1e-6]);
        let ba = metrics(&b, &a, &[1e-6]);
        assert!((ab[0] - ba[1]).abs() < 1e-12); // KL(a||b) computed from either side
        assert!((ab[0] - ab[1]).abs() > 1e-3); // directional KL is asymmetric
    }

    #[test]
    fn wasserstein_moves_mass_by_distance() {
        let a = from_counts(&[(100, 1.0)]);
        let b = from_counts(&[(130, 1.0)]);
        let m = metrics(&a, &b, &[1e-6]);
        assert!((m[m.len() - 1] - 30.0).abs() < 1e-9);
    }
}
