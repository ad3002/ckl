//! `ckl refdep` - reference dependence of whole-chromosome scores.
//! Components are per-assembly histogram tables (`ckl hist` output); a recipe pools the raw gap counts of its
//! components per autosome (the manuscript's consensus); every query region is scored against every recipe.
//! Reports per region and per (query, recipe): symmetric KL and KL(query || ref) with an additive eps, JS distance,
//! and the fraction of query gaps that fall in bins the reference lacks (the mass the eps term prices at ~-ln eps).

use crate::{create, from_counts, metric_names, metrics, open, read_hist, Res};
use std::collections::{BTreeMap, HashMap};
use std::io::{BufRead, Write};

pub struct RefdepArgs {
    pub components: Vec<String>,
    pub recipes: Vec<String>,
    pub queries: Vec<String>,
    pub eps: f64,
    pub out: String,
}

/// region id -> autosome, via an explicit map or the `chrN_...` naming of `ckl hist`.
fn chrom_of(region: &str, map: &Option<HashMap<String, String>>) -> Option<String> {
    let c = match map {
        Some(m) => m.get(region)?.clone(),
        None => region.split('_').next()?.to_string(),
    };
    let n: u32 = c.strip_prefix("chr")?.parse().ok()?;
    (1..=22).contains(&n).then_some(c)
}

struct Component {
    /// (region, chrom, counts)
    regions: Vec<(String, String, Vec<(i64, f64)>)>,
}

fn load_component(spec: &str) -> Res<(String, Component)> {
    let (tag, rest) = spec.split_once('=').ok_or_else(|| format!("component {spec}: expected TAG=PREFIX[@MAP]"))?;
    let (prefix, map_p) = match rest.split_once('@') {
        Some((p, m)) => (p, Some(m)),
        None => (rest, None),
    };
    let map = match map_p {
        Some(p) => {
            let mut m = HashMap::new();
            for line in open(p)?.lines() {
                let line = line?;
                let t: Vec<&str> = line.split_whitespace().collect();
                if t.len() >= 2 {
                    m.insert(t[0].to_string(), t[1].to_string());
                }
            }
            Some(m)
        }
        None => None,
    };
    let hist = read_hist(&format!("{prefix}.hist.tsv"))?;
    let mut regions = Vec::new();
    for (r, counts) in hist {
        if let Some(c) = chrom_of(&r, &map) {
            regions.push((r, c, counts));
        }
    }
    regions.sort_by(|a, b| a.0.cmp(&b.0));
    let chroms: std::collections::BTreeSet<&String> = regions.iter().map(|x| &x.1).collect();
    if chroms.len() != 22 {
        return Err(format!("component {tag}: {} autosomes found (expected 22) - check the region map", chroms.len()).into());
    }
    Ok((tag.to_string(), Component { regions }))
}

pub fn run(a: &RefdepArgs) -> Res<()> {
    let mut comps: HashMap<String, Component> = HashMap::new();
    for s in &a.components {
        let (t, c) = load_component(s)?;
        eprintln!("component {t}: {} autosomal regions", c.regions.len());
        comps.insert(t, c);
    }
    // recipe -> chrom -> pooled counts
    let mut refs: Vec<(String, HashMap<String, Vec<(i64, f64)>>)> = Vec::new();
    for s in &a.recipes {
        let (name, parts) = s.split_once('=').ok_or_else(|| format!("recipe {s}: expected NAME=TAG,TAG"))?;
        let mut pooled: HashMap<String, BTreeMap<i64, f64>> = HashMap::new();
        for t in parts.split(',') {
            let c = comps.get(t).ok_or_else(|| format!("recipe {name}: unknown component {t}"))?;
            for (_, chrom, counts) in &c.regions {
                let e = pooled.entry(chrom.clone()).or_default();
                for &(g, v) in counts {
                    *e.entry(g).or_insert(0.0) += v;
                }
            }
        }
        refs.push((name.to_string(), pooled.into_iter().map(|(k, v)| (k, v.into_iter().collect())).collect()));
    }
    let names = metric_names(&[a.eps]);
    let col = |n: &str| names.iter().position(|x| x == n).ok_or_else(|| format!("metric {n} missing"));
    let i_pq = col(&format!("kl_naive_{:e}_pq", a.eps))?;
    let i_qp = col(&format!("kl_naive_{:e}_qp", a.eps))?;
    let i_js = col("js_dist")?;
    let mut per = create(&format!("{}.regions.tsv", a.out))?;
    writeln!(per, "query\tregion\tchrom\treference\tskl\tkl_query_ref\tjs_dist\tfrac_gaps_absent_from_ref")?;
    let mut summ = create(&format!("{}.summary.tsv", a.out))?;
    writeln!(summ, "query\treference\tn_regions\tmean_skl\tmean_kl_query_ref\tmean_js_dist\tfrac_gaps_absent_from_ref_pooled\tmean_region_frac_absent")?;
    for q in &a.queries {
        let c = comps.get(q).ok_or_else(|| format!("unknown query component {q}"))?;
        for (rname, pooled) in &refs {
            let (mut n, mut s_skl, mut s_kl, mut s_js, mut absent, mut total, mut s_frac) = (0usize, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            for (region, chrom, counts) in &c.regions {
                let rc = pooled.get(chrom).ok_or_else(|| format!("{rname} lacks {chrom}"))?;
                let (p, r) = (from_counts(counts), from_counts(rc));
                let m = metrics(&p, &r, &[a.eps]);
                let skl = 0.5 * (m[i_pq] + m[i_qp]);
                let present: std::collections::HashSet<i64> = rc.iter().filter(|x| x.1 > 0.0).map(|x| x.0).collect();
                let (ab, tot): (f64, f64) = counts.iter().fold((0.0, 0.0), |(ab, tot), &(g, v)| (ab + if present.contains(&g) { 0.0 } else { v }, tot + v));
                writeln!(per, "{q}\t{region}\t{chrom}\t{rname}\t{skl:.6}\t{:.6}\t{:.6}\t{:.6}", m[i_pq], m[i_js], ab / tot)?;
                n += 1;
                s_skl += skl;
                s_kl += m[i_pq];
                s_js += m[i_js];
                absent += ab;
                s_frac += ab / tot;
                total += tot;
            }
            writeln!(summ, "{q}\t{rname}\t{n}\t{:.4}\t{:.4}\t{:.4}\t{:.5}\t{:.5}", s_skl / n as f64, s_kl / n as f64, s_js / n as f64, absent / total, s_frac / n as f64)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chrom_of_keeps_autosomes_only() {
        assert_eq!(chrom_of("chr7_MATERNAL", &None), Some("chr7".into()));
        assert_eq!(chrom_of("chrX_PATERNAL", &None), None);
        assert_eq!(chrom_of("chr23_x", &None), None);
        let m: HashMap<String, String> = [("GWH1".to_string(), "chr2".to_string())].into_iter().collect();
        assert_eq!(chrom_of("GWH1", &Some(m.clone())), Some("chr2".into()));
        assert_eq!(chrom_of("GWH9", &Some(m)), None);
    }
}
