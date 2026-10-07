//! `ckl auroc` - separability of edited / erroneous centromeres from validated ones (reviewer 2, point 10).
//! Input: long TSV with header `set group chrom score`, set = pos | neg. Negatives are validated centromeres scored
//! held-out against their own chromosome; positives are edited or uncorrected arrays scored the same way, grouped
//! (e.g. `collapse|5`). For each positive group, against negatives of the same chromosomes:
//!   auroc_raw         - Mann-Whitney AUROC on the raw score (one scale for all chromosomes)
//!   auroc_chrom_norm  - AUROC after dividing every score by the median negative score of its chromosome
//!   tpr_at_fixed_q95  - fraction of positives above ONE threshold for all chromosomes: the 95th percentile of the
//!                       pooled negatives (false-positive rate 5% by construction)

use crate::{create, open, Res};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, Write};

/// AUROC with mid-ranks for ties: P(pos > neg) + 0.5 P(pos == neg).
pub fn auroc(pos: &[f64], neg: &[f64]) -> f64 {
    if pos.is_empty() || neg.is_empty() {
        return f64::NAN;
    }
    let mut all: Vec<(f64, bool)> = pos.iter().map(|&x| (x, true)).chain(neg.iter().map(|&x| (x, false))).collect();
    all.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut rank_sum_pos = 0.0;
    let mut i = 0;
    while i < all.len() {
        let mut j = i;
        while j + 1 < all.len() && all[j + 1].0 == all[i].0 {
            j += 1;
        }
        let mid = (i + j) as f64 / 2.0 + 1.0;
        for k in i..=j {
            if all[k].1 {
                rank_sum_pos += mid;
            }
        }
        i = j + 1;
    }
    let (np, nn) = (pos.len() as f64, neg.len() as f64);
    (rank_sum_pos - np * (np + 1.0) / 2.0) / (np * nn)
}

fn quantile(v: &mut Vec<f64>, q: f64) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if v.is_empty() {
        return f64::NAN;
    }
    let idx = ((v.len() - 1) as f64 * q).round() as usize;
    v[idx]
}

pub fn run(input: &str, out: &str) -> Res<()> {
    let mut neg: HashMap<String, Vec<f64>> = HashMap::new();
    let mut pos: BTreeMap<String, Vec<(String, f64)>> = BTreeMap::new();
    let mut lines = open(input)?.lines();
    let header = lines.next().ok_or("empty input")??;
    if header.split('\t').collect::<Vec<_>>() != ["set", "group", "chrom", "score"] {
        return Err(format!("header must be: set group chrom score (got {header})").into());
    }
    let mut n_na = 0usize;
    for line in lines {
        let line = line?;
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() != 4 {
            return Err(format!("bad row: {line}").into());
        }
        if t[3] == "NA" {
            n_na += 1;
            continue;
        }
        let s: f64 = t[3].parse()?;
        match t[0] {
            "neg" => neg.entry(t[2].to_string()).or_default().push(s),
            "pos" => pos.entry(t[1].to_string()).or_default().push((t[2].to_string(), s)),
            x => return Err(format!("set must be pos|neg, got {x}").into()),
        }
    }
    let med: HashMap<String, f64> = neg.iter().map(|(c, v)| (c.clone(), quantile(&mut v.clone(), 0.5))).collect();
    let mut w = create(out)?;
    writeln!(w, "group\tn_pos\tn_neg\tauroc_raw\tauroc_chrom_norm\tfixed_threshold_q95\ttpr_at_fixed_q95\tmedian_pos_norm")?;
    for (g, ps) in &pos {
        let chroms: HashSet<&String> = ps.iter().map(|x| &x.0).collect();
        let missing: Vec<&&String> = chroms.iter().filter(|c| !neg.contains_key(**c)).collect();
        if !missing.is_empty() {
            return Err(format!("group {g}: no negatives for {missing:?}").into());
        }
        let negs: Vec<(String, f64)> = chroms.iter().flat_map(|c| neg[*c].iter().map(move |&s| ((*c).clone(), s))).collect();
        let praw: Vec<f64> = ps.iter().map(|x| x.1).collect();
        let nraw: Vec<f64> = negs.iter().map(|x| x.1).collect();
        let pn: Vec<f64> = ps.iter().map(|(c, s)| s / med[c]).collect();
        let nn: Vec<f64> = negs.iter().map(|(c, s)| s / med[c]).collect();
        let thr = quantile(&mut nraw.clone(), 0.95);
        let tpr = praw.iter().filter(|&&s| s > thr).count() as f64 / praw.len() as f64;
        writeln!(w, "{g}\t{}\t{}\t{:.3}\t{:.3}\t{:.4}\t{:.3}\t{:.3}", praw.len(), nraw.len(), auroc(&praw, &nraw), auroc(&pn, &nn), thr, tpr, quantile(&mut pn.clone(), 0.5))?;
    }
    eprintln!("ckl auroc: {} positive groups, {} negative chromosomes, {} NA scores skipped (listed as skipped, not scored)", pos.len(), neg.len(), n_na);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn auroc_extremes_and_ties() {
        assert_eq!(auroc(&[3.0, 4.0], &[1.0, 2.0]), 1.0);
        assert_eq!(auroc(&[1.0, 2.0], &[3.0, 4.0]), 0.0);
        assert_eq!(auroc(&[1.0], &[1.0]), 0.5);
        assert!((auroc(&[2.0, 2.0], &[1.0, 2.0]) - 0.75).abs() < 1e-12);
    }
}
