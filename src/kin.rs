//! `ckl kin` - within-chromosome nearest neighbours for the trio (identity-by-descent) and population-group tests
//!. Reads an all-vs-all distance table (ckl pairwise with --metric, or headerless
//! symmetric rows such as `mash dist` output mapped through --ids), the meta table (individual, group, chrom, hap),
//! optional trios (child father mother) and optional labels (individual label).
//!
//! Candidates for a query are the centromeres of the same chromosome from OTHER individuals. Trio: correct = the
//! nearest candidate belongs to a parent; chance = parent share of the candidates. Population: candidates restricted
//! to other groups (families) with a label; predicted label = the nearest candidate's; chance = the query label's
//! share of those candidates.

use crate::{create, open, read_meta, Meta, Res};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, Write};

pub struct KinArgs {
    pub dist: String,
    pub metric: Option<String>,
    pub ids: Option<String>,
    pub meta: String,
    pub trios: Option<String>,
    pub labels: Option<String>,
    pub exclude: String,
    pub label: String,
    pub out: String,
}

#[derive(Default, Clone, Debug)]
struct Best {
    any: Option<(f64, String)>,       // nearest among other individuals
    labelled: Option<(f64, String)>,  // nearest among other groups whose individual has a label
    parent: Option<(f64, String)>,    // nearest among the query's parents
    unrelated: Option<(f64, String)>, // nearest among other groups
}

fn upd(slot: &mut Option<(f64, String)>, d: f64, r: &str) {
    if slot.as_ref().map(|s| d < s.0).unwrap_or(true) {
        *slot = Some((d, r.to_string()));
    }
}

pub struct Kin {
    meta: HashMap<String, Meta>,
    parents: HashMap<String, (String, String)>,
    labels: HashMap<String, String>,
    exclude: HashSet<String>,
    best: HashMap<String, Best>,
    n_offered: u64,
}

pub struct TrioRow {
    pub child: String,
    pub n: u32,
    pub correct: u32,
    pub chance: f64,
    pub med_parent: f64,
    pub med_unrelated: f64,
    pub med_ratio: f64,
}

pub struct PopRow {
    pub label: String,
    pub n: u32,
    pub correct: u32,
    pub chance: f64,
}

fn median(v: &mut Vec<f64>) -> f64 {
    v.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    if v.is_empty() {
        f64::NAN
    } else if v.len() % 2 == 1 {
        v[v.len() / 2]
    } else {
        (v[v.len() / 2 - 1] + v[v.len() / 2]) / 2.0
    }
}

impl Kin {
    pub fn new(meta: HashMap<String, Meta>, parents: HashMap<String, (String, String)>, labels: HashMap<String, String>, exclude: HashSet<String>) -> Self {
        Kin { meta, parents, labels, exclude, best: HashMap::new(), n_offered: 0 }
    }

    /// One directed distance q -> r.
    pub fn offer(&mut self, q: &str, r: &str, d: f64) {
        let (Some(mq), Some(mr)) = (self.meta.get(q), self.meta.get(r)) else { return };
        if !d.is_finite() || mq.chrom != mr.chrom || mq.individual == mr.individual || self.exclude.contains(&mq.chrom) {
            return;
        }
        self.n_offered += 1;
        let b = self.best.entry(q.to_string()).or_default();
        upd(&mut b.any, d, r);
        if mq.group != mr.group {
            upd(&mut b.unrelated, d, r);
            if self.labels.contains_key(&mr.individual) {
                upd(&mut b.labelled, d, r);
            }
        }
        if let Some((f, m)) = self.parents.get(&mq.individual) {
            if mr.individual == *f || mr.individual == *m {
                upd(&mut b.parent, d, r);
            }
        }
    }

    fn candidates(&self, m: &Meta) -> Vec<&Meta> {
        self.meta.values().filter(|x| x.chrom == m.chrom && x.individual != m.individual).collect()
    }

    fn parent_is_nearest(&self, m: &Meta, b: &Best) -> Option<bool> {
        let (f, mo) = self.parents.get(&m.individual)?;
        let (_, r) = b.any.as_ref()?;
        let ind = &self.meta[r].individual;
        Some(ind == f || ind == mo)
    }

    pub fn trio_rows(&self) -> Vec<TrioRow> {
        let mut acc: BTreeMap<String, (u32, u32, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> = BTreeMap::new();
        for (q, b) in &self.best {
            let m = &self.meta[q];
            let Some(hit) = self.parent_is_nearest(m, b) else { continue };
            let (f, mo) = &self.parents[&m.individual];
            let cands = self.candidates(m);
            let n_par = cands.iter().filter(|x| x.individual == *f || x.individual == *mo).count();
            let e = acc.entry(m.individual.clone()).or_default();
            e.0 += 1;
            e.1 += hit as u32;
            e.2.push(n_par as f64 / cands.len() as f64);
            if let (Some(p), Some(u)) = (&b.parent, &b.unrelated) {
                e.3.push(p.0);
                e.4.push(u.0);
                if u.0 > 0.0 {
                    e.5.push(p.0 / u.0);
                }
            }
        }
        acc.into_iter()
            .map(|(child, (n, c, ch, mut dp, mut du, mut ratio))| TrioRow {
                child,
                n,
                correct: c,
                chance: ch.iter().sum::<f64>() / ch.len() as f64,
                med_parent: median(&mut dp),
                med_unrelated: median(&mut du),
                med_ratio: median(&mut ratio),
            })
            .collect()
    }

    pub fn pop_rows(&self) -> Vec<PopRow> {
        let mut acc: BTreeMap<String, (u32, u32, f64)> = BTreeMap::new();
        for (q, b) in &self.best {
            let m = &self.meta[q];
            let (Some(lab), Some((_, r))) = (self.labels.get(&m.individual), &b.labelled) else { continue };
            let nlab = &self.labels[&self.meta[r].individual];
            let cands: Vec<&Meta> = self.candidates(m).into_iter().filter(|x| x.group != m.group && self.labels.contains_key(&x.individual)).collect();
            let same = cands.iter().filter(|x| self.labels[&x.individual] == *lab).count();
            let e = acc.entry(lab.clone()).or_default();
            e.0 += 1;
            e.1 += (nlab == lab) as u32;
            e.2 += same as f64 / cands.len() as f64;
        }
        acc.into_iter().map(|(label, (n, c, ch))| PopRow { label, n, correct: c, chance: ch / n as f64 }).collect()
    }
}

fn read_two_cols(p: &str) -> Res<Vec<(String, String, Option<String>)>> {
    let mut v = Vec::new();
    for line in open(p)?.lines() {
        let line = line?;
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() >= 2 && !t[0].starts_with('#') {
            v.push((t[0].to_string(), t[1].to_string(), t.get(2).map(|s| s.to_string())));
        }
    }
    Ok(v)
}

pub fn run(a: &KinArgs) -> Res<()> {
    let meta: HashMap<String, Meta> = read_meta(&a.meta)?.into_iter().collect();
    let exclude: HashSet<String> = a.exclude.split(',').filter(|s| !s.is_empty()).map(String::from).collect();
    let mut parents = HashMap::new();
    if let Some(p) = &a.trios {
        for (c, f, m) in read_two_cols(p)? {
            if c == "child" {
                continue;
            }
            let m = m.ok_or_else(|| format!("trio row for {c} lacks a mother column"))?;
            for x in [&c, &f, &m] {
                if !meta.values().any(|y| y.individual == *x) {
                    return Err(format!("trio individual {x} absent from meta").into());
                }
            }
            parents.insert(c, (f, m));
        }
    }
    let mut labels = HashMap::new();
    if let Some(p) = &a.labels {
        for (i, l, _) in read_two_cols(p)? {
            if i != "individual" {
                labels.insert(i, l);
            }
        }
    }
    let mut idmap = HashMap::new();
    if let Some(p) = &a.ids {
        for (path, id, _) in read_two_cols(p)? {
            idmap.insert(path, id);
        }
    }
    let map = |x: &str| -> String { idmap.get(x).cloned().unwrap_or_else(|| x.to_string()) };
    let mut kin = Kin::new(meta, parents, labels, exclude);
    let mut lines = open(&a.dist)?.lines();
    let mut n_rows = 0u64;
    match &a.metric {
        Some(mname) => {
            let header = lines.next().ok_or("empty distance table")??;
            let cols: Vec<&str> = header.split('\t').collect();
            let ci = cols.iter().position(|c| c == mname).ok_or_else(|| format!("metric {mname} not in header: {header}"))?;
            for line in lines {
                let line = line?;
                let t: Vec<&str> = line.split('\t').collect();
                if t.len() <= ci {
                    return Err(format!("short row: {line}").into());
                }
                n_rows += 1;
                if t[ci] == "NA" {
                    continue;
                }
                let d: f64 = t[ci].parse()?;
                kin.offer(&map(t[0]), &map(t[1]), d);
            }
        }
        None => {
            for line in lines {
                let line = line?;
                let t: Vec<&str> = line.split('\t').collect();
                if t.len() < 3 {
                    return Err(format!("short row: {line}").into());
                }
                n_rows += 1;
                let d: f64 = t[2].parse()?;
                let (x, y) = (map(t[0]), map(t[1]));
                kin.offer(&x, &y, d);
                kin.offer(&y, &x, d);
            }
        }
    }
    if kin.best.is_empty() {
        return Err(format!("no same-chromosome pair between different individuals found in {} rows of {} (id mapping?)", n_rows, a.dist).into());
    }
    let mut nn = create(&format!("{}.nn.tsv", a.out))?;
    writeln!(nn, "method\tquery\tindividual\tgroup\tchrom\thap\tnn_any\tnn_any_individual\td_any\tnn_labelled\tnn_labelled_individual\td_labelled\td_best_parent\td_best_unrelated\tparent_is_nearest\tlabel\tnn_label")?;
    let mut qs: Vec<&String> = kin.best.keys().collect();
    qs.sort();
    let show = |s: &Option<(f64, String)>| -> (String, String, String) {
        s.as_ref().map(|(d, r)| (r.clone(), kin.meta[r].individual.clone(), format!("{d}"))).unwrap_or(("NA".into(), "NA".into(), "NA".into()))
    };
    for q in qs {
        let b = &kin.best[q];
        let m = &kin.meta[q];
        let (na, ia, da) = show(&b.any);
        let (nl, il, dl) = show(&b.labelled);
        let dp = b.parent.as_ref().map(|x| format!("{}", x.0)).unwrap_or("NA".into());
        let du = b.unrelated.as_ref().map(|x| format!("{}", x.0)).unwrap_or("NA".into());
        let pin = kin.parent_is_nearest(m, b).map(|x| x.to_string()).unwrap_or("NA".into());
        let lab = kin.labels.get(&m.individual).cloned().unwrap_or("NA".into());
        let nlab = kin.labels.get(&il).cloned().unwrap_or("NA".into());
        writeln!(nn, "{}\t{q}\t{}\t{}\t{}\t{}\t{na}\t{ia}\t{da}\t{nl}\t{il}\t{dl}\t{dp}\t{du}\t{pin}\t{lab}\t{nlab}", a.label, m.individual, m.group, m.chrom, m.hap)?;
    }
    let trio = kin.trio_rows();
    if !trio.is_empty() {
        let mut w = create(&format!("{}.trio.tsv", a.out))?;
        writeln!(w, "method\tchild\tn_centromeres\tparent_is_nearest\trate\tchance\tmedian_d_best_parent\tmedian_d_best_unrelated\tmedian_ratio_parent_over_unrelated")?;
        let (mut tn, mut tc, mut tch) = (0u32, 0u32, 0.0);
        for r in &trio {
            writeln!(w, "{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.6}\t{:.6}\t{:.4}", a.label, r.child, r.n, r.correct, r.correct as f64 / r.n as f64, r.chance, r.med_parent, r.med_unrelated, r.med_ratio)?;
            tn += r.n;
            tc += r.correct;
            tch += r.chance * r.n as f64;
        }
        writeln!(w, "{}\tALL\t{tn}\t{tc}\t{:.4}\t{:.4}\tNA\tNA\tNA", a.label, tc as f64 / tn as f64, tch / tn as f64)?;
    }
    let pop = kin.pop_rows();
    if !pop.is_empty() {
        let mut w = create(&format!("{}.pop.tsv", a.out))?;
        writeln!(w, "method\tlabel\tn\tnn_same_label\taccuracy\tchance")?;
        let (mut tn, mut tc, mut tch) = (0u32, 0u32, 0.0);
        for r in &pop {
            writeln!(w, "{}\t{}\t{}\t{}\t{:.4}\t{:.4}", a.label, r.label, r.n, r.correct, r.correct as f64 / r.n as f64, r.chance)?;
            tn += r.n;
            tc += r.correct;
            tch += r.chance * r.n as f64;
        }
        writeln!(w, "{}\tALL\t{tn}\t{tc}\t{:.4}\t{:.4}", a.label, tc as f64 / tn as f64, tch / tn as f64)?;
    }
    eprintln!(
        "ckl kin {}: {} rows read, {} same-chromosome cross-individual distances, {} queries; trio children {}; labelled queries {}",
        a.label,
        n_rows,
        kin.n_offered,
        kin.best.len(),
        trio.len(),
        pop.iter().map(|r| r.n).sum::<u32>()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(ind: &str, group: &str, chrom: &str) -> Meta {
        Meta { individual: ind.into(), group: group.into(), hap: "h1".into(), chrom: chrom.into(), role: "x".into() }
    }

    fn panel() -> Kin {
        // family F: child C with parents P1, P2; unrelated U1 (AFR), U2 (EUR); C labelled EUR, P1 EUR, P2 EUR
        let meta: HashMap<String, Meta> = [
            ("c1", m("C", "F", "chr1")),
            ("p1", m("P1", "F", "chr1")),
            ("p2", m("P2", "F", "chr1")),
            ("u1", m("U1", "G", "chr1")),
            ("u2", m("U2", "H", "chr1")),
            ("y", m("U2", "H", "chrY")),
        ]
        .into_iter()
        .map(|(k, v)| (k.to_string(), v))
        .collect();
        let parents = [("C".to_string(), ("P1".to_string(), "P2".to_string()))].into_iter().collect();
        let labels = [("C", "EUR"), ("P1", "EUR"), ("P2", "EUR"), ("U1", "AFR"), ("U2", "EUR")].into_iter().map(|(a, b)| (a.to_string(), b.to_string())).collect();
        Kin::new(meta, parents, labels, ["chrY".to_string()].into_iter().collect())
    }

    #[test]
    fn parent_nearest_is_counted_against_parent_share() {
        let mut k = panel();
        k.offer("c1", "p1", 0.01);
        k.offer("c1", "p2", 0.5);
        k.offer("c1", "u1", 0.2);
        k.offer("c1", "u2", 0.3);
        k.offer("c1", "c1", 0.0); // self: ignored
        k.offer("y", "u1", 0.0); // chrY excluded, different chromosome anyway
        let t = k.trio_rows();
        assert_eq!(t.len(), 1);
        assert_eq!((t[0].n, t[0].correct), (1, 1));
        assert!((t[0].chance - 0.5).abs() < 1e-12); // 2 parents of 4 candidates
        assert!((t[0].med_parent - 0.01).abs() < 1e-12);
        assert!((t[0].med_unrelated - 0.2).abs() < 1e-12);
        assert!((t[0].med_ratio - 0.05).abs() < 1e-12);
    }

    #[test]
    fn population_uses_other_labelled_groups_only() {
        let mut k = panel();
        // u1 (AFR): labelled candidates in other groups = c1, p1, p2 (EUR) and u2 (EUR); nearest u2 -> wrong; chance 0
        k.offer("u1", "c1", 0.4);
        k.offer("u1", "p1", 0.5);
        k.offer("u1", "p2", 0.6);
        k.offer("u1", "u2", 0.1);
        // c1 (EUR): parents are the same group and never count; u1 (AFR) nearer than u2 (EUR) -> wrong; chance 0.5
        k.offer("c1", "p1", 0.01);
        k.offer("c1", "u1", 0.2);
        k.offer("c1", "u2", 0.3);
        let p = k.pop_rows();
        let afr = p.iter().find(|r| r.label == "AFR").unwrap();
        let eur = p.iter().find(|r| r.label == "EUR").unwrap();
        assert_eq!((afr.n, afr.correct), (1, 0));
        assert!((afr.chance - 0.0).abs() < 1e-12);
        assert_eq!((eur.n, eur.correct), (1, 0));
        assert!((eur.chance - 0.5).abs() < 1e-12);
    }

    #[test]
    fn median_is_middle_or_mean_of_middles() {
        assert_eq!(median(&mut vec![3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&mut vec![4.0, 1.0, 3.0, 2.0]), 2.5);
        assert!(median(&mut vec![]).is_nan());
    }
}
