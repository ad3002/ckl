//! Cross-species identification. All analysis in Rust:
//! `censat-regions` CenSat active_hor BED -> regions / meta / family tables (labels from the `hsaN` contig suffix)
//! `extract`        region sequences -> one FASTA per region (for Mash), with a balanced ledger
//! `meta-set`       rewrite a meta table (e.g. the Gao human panel as one group `human`) + family by chromosome
//! `xident`         1-NN identification for every ordered pair of units from any distance table, or the size-only null

use crate::{create, open, read_meta, Meta, Res};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::{BufRead, Write};

fn fai_lengths(paths: &[String]) -> Res<HashMap<String, u64>> {
    let mut m = HashMap::new();
    for p in paths {
        for line in open(p)?.lines() {
            let line = line?;
            let t: Vec<&str> = line.split('\t').collect();
            if t.len() >= 2 {
                m.insert(t[0].to_string(), t[1].parse()?);
            }
        }
    }
    Ok(m)
}

/// `chr12_hap1_hsa2a` -> ("hap1", "chr2a")
fn parse_ape_contig(c: &str) -> Option<(String, String)> {
    let parts: Vec<&str> = c.split('_').collect();
    if parts.len() < 3 {
        return None;
    }
    let label = parts[parts.len() - 1].strip_prefix("hsa")?;
    Some((parts[1].to_string(), format!("chr{label}")))
}

pub fn run_censat_regions(censat: &str, fai: &[String], species: &str, out: &str) -> Res<()> {
    let lens = fai_lengths(fai)?;
    let mut regions: Vec<(String, u64, u64, String)> = Vec::new();
    let mut meta: BTreeMap<String, (String, String, String)> = BTreeMap::new(); // region -> (hap, chrom, family)
    let mut skipped: BTreeMap<String, String> = BTreeMap::new();
    for (ln, line) in open(censat)?.lines().enumerate() {
        let line = line?;
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() < 4 || !t[3].starts_with("active_hor(") {
            continue;
        }
        let contig = t[0];
        let (s, e): (u64, u64) = (t[1].parse()?, t[2].parse()?);
        let Some(len) = lens.get(contig) else {
            skipped.insert(contig.to_string(), "contig not in the given .fai (other haplotype file?)".into());
            continue;
        };
        if e > *len || e <= s {
            return Err(format!("{censat}:{}: interval {s}-{e} outside contig {contig} length {len}", ln + 1).into());
        }
        let Some((hap, chrom)) = parse_ape_contig(contig) else {
            skipped.insert(contig.to_string(), "no _hsaN suffix".into());
            continue;
        };
        let fam = t[3].trim_start_matches("active_hor(").trim_end_matches(')').split(',').next().unwrap_or("NA").to_string();
        let rid = format!("{species}|{contig}");
        regions.push((contig.to_string(), s, e, rid.clone()));
        let entry = meta.entry(rid).or_insert((hap, chrom, fam.clone()));
        if entry.2 != fam {
            entry.2 = format!("{}+{}", entry.2, fam);
        }
    }
    let mut w = create(&format!("{out}.regions.tsv"))?;
    for (c, s, e, r) in &regions {
        writeln!(w, "{c}\t{s}\t{e}\t{r}")?;
    }
    let mut m = create(&format!("{out}.meta.tsv"))?;
    writeln!(m, "region_id\tindividual\tgroup\thap\tchrom\trole")?;
    let mut f = create(&format!("{out}.family.tsv"))?;
    for (r, (hap, chrom, fam)) in &meta {
        writeln!(m, "{r}\t{species}_{hap}\t{species}\t{hap}\t{chrom}\tboth")?;
        writeln!(f, "{r}\t{fam}")?;
    }
    let mut sk = create(&format!("{out}.skipped.tsv"))?;
    writeln!(sk, "contig\treason")?;
    for (c, why) in &skipped {
        writeln!(sk, "{c}\t{why}")?;
    }
    eprintln!("ckl censat-regions {species}: {} intervals, {} regions, {} contigs skipped", regions.len(), meta.len(), skipped.len());
    Ok(())
}

pub fn run_extract(fasta: &[String], regions: &str, out_dir: &str, ids_out: &str) -> Res<()> {
    let mut want: BTreeMap<String, Vec<(String, usize, usize)>> = BTreeMap::new();
    let mut contigs: HashSet<String> = HashSet::new();
    for line in open(regions)?.lines() {
        let line = line?;
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() < 4 {
            continue;
        }
        want.entry(t[3].to_string()).or_default().push((t[0].to_string(), t[1].parse()?, t[2].parse()?));
        contigs.insert(t[0].to_string());
    }
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    for p in fasta {
        let mut cur: Option<String> = None;
        for line in open(p)?.lines() {
            let line = line?;
            if let Some(h) = line.strip_prefix('>') {
                let name = h.split_whitespace().next().unwrap_or("").to_string();
                cur = if contigs.contains(&name) { Some(name) } else { None };
                if let Some(n) = &cur {
                    seqs.insert(n.clone(), Vec::new());
                }
            } else if let Some(n) = &cur {
                if let Some(s) = seqs.get_mut(n) {
                    s.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
                }
            }
        }
    }
    std::fs::create_dir_all(out_dir)?;
    let mut ids = create(ids_out)?;
    let mut written = 0usize;
    for (i, (rid, ivs)) in want.iter().enumerate() {
        let path = format!("{out_dir}/{i:05}.fa");
        let mut w = create(&path)?;
        for (c, s, e) in ivs {
            let seq = seqs.get(c).ok_or_else(|| format!("region {rid}: contig {c} not found in the FASTA files"))?;
            if *e > seq.len() {
                return Err(format!("region {rid}: {c}:{s}-{e} beyond length {}", seq.len()).into());
            }
            writeln!(w, ">{rid}:{s}-{e}")?;
            for chunk in seq[*s..*e].chunks(80) {
                w.write_all(chunk)?;
                w.write_all(b"\n")?;
            }
        }
        writeln!(ids, "{path}\t{rid}")?;
        written += 1;
    }
    if written != want.len() {
        return Err(format!("extraction ledger does not balance: requested {} written {written}", want.len()).into());
    }
    eprintln!("ckl extract: {written} regions written");
    Ok(())
}

pub fn run_meta_set(meta: &str, group: &str, family_map: &str, role: Option<&str>, out: &str) -> Res<()> {
    let mut fam: HashMap<String, String> = HashMap::new();
    for line in open(family_map)?.lines() {
        let line = line?;
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() >= 2 {
            fam.insert(t[0].to_string(), t[1].to_string());
        }
    }
    let rows = read_meta(meta)?;
    let mut m = create(&format!("{out}.meta.tsv"))?;
    let mut f = create(&format!("{out}.family.tsv"))?;
    writeln!(m, "region_id\tindividual\tgroup\thap\tchrom\trole")?;
    for (rid, r) in &rows {
        writeln!(m, "{rid}\t{}\t{group}\t{}\t{}\t{}", r.individual, r.hap, r.chrom, role.unwrap_or(&r.role))?;
        writeln!(f, "{rid}\t{}", fam.get(&r.chrom).map(|s| s.as_str()).unwrap_or("NA"))?;
    }
    eprintln!("ckl meta-set: {} rows -> group {group}", rows.len());
    Ok(())
}

struct Rng(u64);
impl Rng {
    fn below(&mut self, n: usize) -> usize {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        (self.0.wrapping_mul(0x2545F4914F6CDD1D) % n.max(1) as u64) as usize
    }
}

pub struct XidentArgs {
    pub dist: Option<String>,
    pub metric: Option<String>,
    pub ids: Option<String>,
    pub size: Vec<String>,
    pub meta: Vec<String>,
    pub family: Vec<String>,
    pub by: String,
    pub query_groups: Option<String>,
    pub queries_from: Option<String>,
    pub exclude: String,
    pub boot: usize,
    pub seed: u64,
    pub label: String,
    pub out: String,
}

pub fn run_xident(a: &XidentArgs) -> Res<()> {
    let mut meta: HashMap<String, Meta> = HashMap::new();
    for p in &a.meta {
        for (r, m) in read_meta(p)? {
            meta.insert(r, m);
        }
    }
    let mut fam: HashMap<String, String> = HashMap::new();
    for p in &a.family {
        for line in open(p)?.lines() {
            let line = line?;
            let t: Vec<&str> = line.split('\t').collect();
            if t.len() >= 2 {
                fam.insert(t[0].to_string(), t[1].to_string());
            }
        }
    }
    let unit = |r: &str| -> Option<String> {
        meta.get(r).map(|m| if a.by == "individual" { m.individual.clone() } else { m.group.clone() })
    };
    let qgroups: Option<HashSet<String>> = a.query_groups.as_ref().map(|s| s.split(',').map(String::from).collect());
    // optional restriction to a query set (first column of a file, e.g. another method's .detail.tsv queries)
    let qset: Option<HashSet<String>> = match &a.queries_from {
        Some(p) => {
            let mut h = HashSet::new();
            for line in open(p)?.lines() {
                let line = line?;
                let t: Vec<&str> = line.split('\t').collect();
                let id = if t.len() > 1 && t[0].starts_with("spacing_") || t.len() > 1 && t[0].starts_with("mash_") || t.len() > 1 && t[0] == "size" { t[1] } else { t[0] };
                h.insert(id.to_string());
            }
            Some(h)
        }
        None => None,
    };
    let is_query = |r: &str| -> bool {
        qset.as_ref().map(|q| q.contains(r)).unwrap_or(true)
            && meta.get(r).map(|m| qgroups.as_ref().map(|q| q.contains(&m.group)).unwrap_or(true)).unwrap_or(false)
    };
    let exclude: HashSet<String> = a.exclude.split(',').filter(|s| !s.is_empty()).map(String::from).collect();
    // best[(query, target unit)] = (distance, ref)
    let mut best: HashMap<(String, String), (f64, String)> = HashMap::new();
    let offer = |q: &str, r: &str, d: f64, best: &mut HashMap<(String, String), (f64, String)>| {
        if !d.is_finite() || !is_query(q) {
            return;
        }
        let (Some(uq), Some(ur)) = (unit(q), unit(r)) else { return };
        if uq == ur {
            return;
        }
        let e = best.entry((q.to_string(), ur)).or_insert((f64::INFINITY, String::new()));
        if d < e.0 {
            *e = (d, r.to_string());
        }
    };
    if !a.size.is_empty() {
        // size-only null: standardized Euclidean on (ln region_bp, ln n_boxes)
        let mut feat: HashMap<String, (f64, f64)> = HashMap::new();
        for p in &a.size {
            for line in open(p)?.lines().skip(1) {
                let line = line?;
                let t: Vec<&str> = line.split('\t').collect();
                if t.len() >= 4 && meta.contains_key(t[0]) {
                    let bp: f64 = t[2].parse()?;
                    let nb: f64 = t[3].parse()?;
                    if bp > 0.0 && nb >= 50.0 {
                        feat.insert(t[0].to_string(), (bp.ln(), nb.ln()));
                    }
                }
            }
        }
        // per-unit descriptive summary (box count, active-HOR bp)
        let mut per: BTreeMap<String, Vec<(f64, f64)>> = BTreeMap::new();
        for (r, v) in &feat {
            if let Some(u) = unit(r) {
                per.entry(u).or_default().push((v.1.exp(), v.0.exp()));
            }
        }
        let mut uw = create(&format!("{}.units.tsv", a.out))?;
        writeln!(uw, "unit\tn_regions\tmedian_boxes\tmedian_bp")?;
        for (u, v) in per.iter_mut() {
            let mut b: Vec<f64> = v.iter().map(|x| x.0).collect();
            let mut l: Vec<f64> = v.iter().map(|x| x.1).collect();
            b.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
            l.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
            writeln!(uw, "{u}\t{}\t{:.0}\t{:.0}", v.len(), b[b.len() / 2], l[l.len() / 2])?;
        }
        let n = feat.len() as f64;
        let (m0, m1) = feat.values().fold((0.0, 0.0), |s, v| (s.0 + v.0 / n, s.1 + v.1 / n));
        let (s0, s1) = feat.values().fold((0.0, 0.0), |s, v| (s.0 + (v.0 - m0).powi(2) / n, s.1 + (v.1 - m1).powi(2) / n));
        let (s0, s1) = (s0.sqrt(), s1.sqrt());
        let qs: Vec<(&String, &(f64, f64))> = feat.iter().filter(|(r, _)| is_query(r)).collect();
        for (q, fq) in &qs {
            for (r, fr) in &feat {
                let d = (((fq.0 - fr.0) / s0).powi(2) + ((fq.1 - fr.1) / s1).powi(2)).sqrt();
                offer(q, r, d, &mut best);
            }
        }
    } else {
        let dist = a.dist.as_ref().ok_or("xident: give --dist or --size")?;
        let mut idmap: HashMap<String, String> = HashMap::new();
        if let Some(p) = &a.ids {
            for line in open(p)?.lines() {
                let line = line?;
                let t: Vec<&str> = line.split('\t').collect();
                if t.len() >= 2 {
                    idmap.insert(t[0].to_string(), t[1].to_string());
                }
            }
        }
        let map = |x: &str| -> String { idmap.get(x).cloned().unwrap_or_else(|| x.to_string()) };
        let mut lines = open(dist)?.lines();
        match &a.metric {
            Some(mname) => {
                let header = lines.next().ok_or("empty distance table")??;
                let cols: Vec<&str> = header.split('\t').collect();
                let ci = cols.iter().position(|c| c == mname).ok_or_else(|| format!("metric {mname} not in header"))?;
                for line in lines {
                    let line = line?;
                    let t: Vec<&str> = line.split('\t').collect();
                    if t.len() <= ci || t[ci] == "NA" {
                        continue;
                    }
                    let d: f64 = t[ci].parse()?;
                    offer(&map(t[0]), &map(t[1]), d, &mut best); // directional: column 1 is the query
                }
            }
            None => {
                for line in lines {
                    let line = line?;
                    let t: Vec<&str> = line.split('\t').collect();
                    if t.len() < 3 {
                        continue;
                    }
                    let d: f64 = t[2].parse()?;
                    let (x, y) = (map(t[0]), map(t[1]));
                    offer(&x, &y, d, &mut best); // symmetric (Mash): both orientations
                    offer(&y, &x, d, &mut best);
                }
            }
        }
    }
    // labels present per target unit
    let mut labels_in: HashMap<String, HashSet<String>> = HashMap::new();
    for (r, m) in &meta {
        if let Some(u) = unit(r) {
            labels_in.entry(u).or_default().insert(m.chrom.clone());
        }
    }
    // (query unit, target unit) -> label -> (n, correct, family agree)
    let mut agg: BTreeMap<(String, String), BTreeMap<String, (u32, u32, u32)>> = BTreeMap::new();
    let mut det = create(&format!("{}.detail.tsv", a.out))?;
    writeln!(det, "method\tquery\tquery_unit\ttarget_unit\tlabel\tpredicted_ref\tpredicted_label\tdistance\tcorrect\tquery_family\tpredicted_family")?;
    let mut keys: Vec<&(String, String)> = best.keys().collect();
    keys.sort();
    for k in keys {
        let (q, tu) = k;
        let (d, r) = &best[k];
        let (mq, mr) = (&meta[q], &meta[r]);
        if exclude.contains(&mq.chrom) || !labels_in.get(tu).map(|s| s.contains(&mq.chrom)).unwrap_or(false) {
            continue;
        }
        let ok = mq.chrom == mr.chrom;
        let fq = fam.get(q).cloned().unwrap_or_else(|| "NA".into());
        let fr = fam.get(r).cloned().unwrap_or_else(|| "NA".into());
        let fa = fq != "NA" && fq == fr;
        let e = agg.entry((unit(q).unwrap_or_default(), tu.clone())).or_default().entry(mq.chrom.clone()).or_insert((0, 0, 0));
        e.0 += 1;
        e.1 += ok as u32;
        e.2 += fa as u32;
        writeln!(det, "{}\t{q}\t{}\t{tu}\t{}\t{r}\t{}\t{d}\t{ok}\t{fq}\t{fr}", a.label, unit(q).unwrap_or_default(), mq.chrom, mr.chrom)?;
    }
    let mut rng = Rng(a.seed | 1);
    let mut w = create(&format!("{}.tsv", a.out))?;
    writeln!(w, "method\tquery_unit\ttarget_unit\tn_queries\tn_labels\ttop1\tci_lo\tci_hi\tfamily_agree")?;
    for ((qu, tu), per) in &agg {
        let labs: Vec<&(u32, u32, u32)> = per.values().collect();
        let n: u32 = labs.iter().map(|x| x.0).sum();
        let c: u32 = labs.iter().map(|x| x.1).sum();
        let f: u32 = labs.iter().map(|x| x.2).sum();
        let mut boots: Vec<f64> = (0..a.boot)
            .map(|_| {
                let (mut bn, mut bc) = (0u32, 0u32);
                for _ in 0..labs.len() {
                    let x = labs[rng.below(labs.len())];
                    bn += x.0;
                    bc += x.1;
                }
                bc as f64 / bn.max(1) as f64
            })
            .collect();
        boots.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        let q = |p: f64| boots.get(((boots.len() as f64 - 1.0) * p).round() as usize).copied().unwrap_or(f64::NAN);
        writeln!(w, "{}\t{qu}\t{tu}\t{n}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}", a.label, labs.len(), c as f64 / n as f64, q(0.025), q(0.975), f as f64 / n as f64)?;
    }
    let units: BTreeSet<&String> = agg.keys().map(|k| &k.0).collect();
    eprintln!("ckl xident {}: {} query units, {} unit pairs", a.label, units.len(), agg.len());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_ape_contig_labels() {
        assert_eq!(parse_ape_contig("chr12_hap1_hsa2a"), Some(("hap1".into(), "chr2a".into())));
        assert_eq!(parse_ape_contig("chr2_mat_hsa3"), Some(("mat".into(), "chr3".into())));
        assert_eq!(parse_ape_contig("chrX_hap2_hsaX"), Some(("hap2".into(), "chrX".into())));
        assert_eq!(parse_ape_contig("chr1_hap1"), None);
    }
}
