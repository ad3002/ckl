//! `ckl perturb` - error-class battery.
//! Edits real centromere sequences in memory at HOR-copy boundaries, re-scans CENP-B boxes with a built-in scanner
//! (NTTCGNNNNANNCGGGN, both strands; count-identity with rust_motif_scan is checked by the `none` class), rebuilds the
//! GCP Model 1 histogram, and scores it against its unedited source. Histograms are written for LOIO scoring.

use crate::{accumulate, create, from_counts, metrics, open, Acc, Res};
use std::collections::{BTreeMap, HashMap};
use std::io::{BufRead, Write};

pub struct Args {
    pub sources: String,
    pub units: String,
    pub donors: String,
    pub donor_meta: String,
    pub sf: String,
    pub doses: String,
    pub placements: u32,
    pub seed: u64,
    pub out: String,
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        // xorshift64*
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % (n.max(1) as u64)) as usize
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn read_fasta(p: &str) -> Res<Vec<(String, Vec<u8>)>> {
    let mut v: Vec<(String, Vec<u8>)> = Vec::new();
    for line in open(p)?.lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            v.push((h.split_whitespace().next().unwrap_or("").to_string(), Vec::new()));
        } else if let Some(last) = v.last_mut() {
            last.1.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        } else if !line.trim().is_empty() {
            return Err(format!("{p}: sequence before first header").into());
        }
    }
    Ok(v)
}

/// Both-strand scan of NTTCGNNNNANNCGGGN; returns (start, end, forward) sorted by start.
pub fn scan_cenpb(s: &[u8]) -> Vec<(i64, i64, bool)> {
    let mut v = Vec::new();
    if s.len() < 17 {
        return v;
    }
    for i in 0..=s.len() - 17 {
        let w = &s[i..i + 17];
        if w[1] == b'T' && w[2] == b'T' && w[3] == b'C' && w[4] == b'G' && w[9] == b'A' && w[12] == b'C' && w[13] == b'G' && w[14] == b'G' && w[15] == b'G' {
            v.push((i as i64, i as i64 + 17, true));
        }
        if w[1] == b'C' && w[2] == b'C' && w[3] == b'C' && w[4] == b'G' && w[7] == b'T' && w[12] == b'C' && w[13] == b'G' && w[14] == b'A' && w[15] == b'A' {
            v.push((i as i64, i as i64 + 17, false));
        }
    }
    v
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        })
        .collect()
}

fn histogram(seq: &[u8]) -> (BTreeMap<i64, u64>, Acc) {
    let boxes = scan_cenpb(seq);
    let mut acc = Acc::default();
    accumulate(&mut acc, &boxes, 0, i64::MAX, 20000);
    (acc.hist.clone(), acc)
}

fn to_dist(h: &BTreeMap<i64, u64>) -> Option<crate::Dist> {
    if h.is_empty() {
        None
    } else {
        Some(from_counts(&h.iter().map(|(g, c)| (*g, *c as f64)).collect::<Vec<_>>()))
    }
}

fn l1_half(a: &BTreeMap<i64, u64>, b: &BTreeMap<i64, u64>) -> u64 {
    let mut keys: Vec<&i64> = a.keys().chain(b.keys()).collect();
    keys.sort();
    keys.dedup();
    let s: u64 = keys.iter().map(|k| (*a.get(k).unwrap_or(&0) as i64 - *b.get(k).unwrap_or(&0) as i64).unsigned_abs()).sum();
    s / 2
}

const CLASSES: [&str; 10] = [
    "none", "collapse", "duplication", "reorder", "inversion", "join_same_sf", "join_cross_sf", "foreign_insertion", "indels", "box_destroy",
];

pub fn run(a: &Args) -> Res<()> {
    let sources = read_fasta(&a.sources)?;
    let donors = read_fasta(&a.donors)?;
    let mut dchrom: HashMap<String, String> = HashMap::new();
    for line in open(&a.donor_meta)?.lines() {
        let line = line?;
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() >= 2 {
            dchrom.insert(t[0].to_string(), t[1].to_string());
        }
    }
    let mut sf: HashMap<String, String> = HashMap::new();
    for line in open(&a.sf)?.lines() {
        let line = line?;
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() >= 2 {
            sf.insert(t[0].to_string(), t[1].to_string());
        }
    }
    let mut units: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    for line in open(&a.units)?.lines() {
        let line = line?;
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() < 3 {
            continue;
        }
        units.entry(t[0].to_string()).or_default().push((t[1].parse()?, t[2].parse()?));
    }
    for v in units.values_mut() {
        v.sort();
    }
    let doses: Vec<usize> = a.doses.split(',').map(|s| s.trim().parse::<usize>()).collect::<Result<_, _>>()?;
    let mut donor_by_chrom: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, (name, _)) in donors.iter().enumerate() {
        let c = dchrom.get(name).ok_or_else(|| format!("donor {name} missing from donor meta"))?;
        donor_by_chrom.entry(c.clone()).or_default().push(i);
    }

    let mut ew = create(&format!("{}.edits.tsv", a.out))?;
    let mut hw = create(&format!("{}.hist.tsv", a.out))?;
    let mut sw = create(&format!("{}.summary.tsv", a.out))?;
    writeln!(ew, "edit_id\tsource\tchrom\tclass\tdose\tplacement\tbp_changed\tlength\tn_boxes\tgaps_changed\tjs_vs_source\tkl_pq_vs_source\tstatus")?;
    writeln!(hw, "region_id\tgap\tcount")?;
    writeln!(sw, "region_id\tn_intervals\tregion_bp\tn_boxes\tn_fwd\tn_rev\tn_gaps\tn_overlap\tmax_gap\tstatus")?;
    let mut rng = Rng(a.seed | 1);
    for (sname, seq) in &sources {
        let chrom = sname.split('_').nth(1).unwrap_or("").trim_start_matches("rc-").to_string();
        let my_sf = sf.get(&chrom).cloned().unwrap_or_default();
        let u = units.get(sname).cloned().unwrap_or_default();
        let (src_h, _) = histogram(seq);
        let src_d = to_dist(&src_h);
        for class in CLASSES {
            let dose_list: Vec<usize> = if class == "none" { vec![0] } else { doses.clone() };
            for &dose in &dose_list {
                let n_place = if class == "none" { 1 } else { a.placements };
                for pl in 0..n_place {
                    let id = format!("{sname}|{class}|{dose}|{pl}");
                    let edited: Result<(Vec<u8>, usize), String> = (|| {
                        let k = dose;
                        let need = match class {
                            "reorder" => 2 * k + 1,
                            "none" | "indels" | "box_destroy" | "foreign_insertion" => 0,
                            _ => k,
                        };
                        if u.len() < need.max(1) && class != "none" && class != "indels" && class != "box_destroy" {
                            return Err(format!("skipped:needs {need} HOR copies, has {}", u.len()));
                        }
                        let pick_block = |rng: &mut Rng, k: usize| -> (usize, usize) {
                            let j = rng.below(u.len() - k + 1);
                            (u[j].0, u[j + k - 1].1)
                        };
                        let donor_seg = |rng: &mut Rng, want_same_sf: Option<bool>, len: usize| -> Result<Vec<u8>, String> {
                            let pool: Vec<usize> = donor_by_chrom
                                .iter()
                                .filter(|(c, _)| {
                                    **c != chrom
                                        && match want_same_sf {
                                            Some(true) => sf.get(*c) == Some(&my_sf),
                                            Some(false) => sf.get(*c).map(|x| x != &my_sf).unwrap_or(false),
                                            None => true,
                                        }
                                })
                                .flat_map(|(_, v)| v.iter().copied())
                                .filter(|&i| donors[i].1.len() > len)
                                .collect();
                            if pool.is_empty() {
                                return Err("skipped:no donor".into());
                            }
                            let d = &donors[pool[rng.below(pool.len())]].1;
                            let off = rng.below(d.len() - len);
                            Ok(d[off..off + len].to_vec())
                        };
                        Ok(match class {
                            "none" => (seq.clone(), 0),
                            "collapse" => {
                                let (x, y) = pick_block(&mut rng, k);
                                ([&seq[..x], &seq[y..]].concat(), y - x)
                            }
                            "duplication" => {
                                let (x, y) = pick_block(&mut rng, k);
                                ([&seq[..y], &seq[x..y], &seq[y..]].concat(), y - x)
                            }
                            "reorder" => {
                                let j = rng.below(u.len() - 2 * k + 1);
                                let (a1, b1) = (u[j].0, u[j + k - 1].1);
                                let (a2, b2) = (u[j + k].0, u[j + 2 * k - 1].1);
                                ([&seq[..a1], &seq[a2..b2], &seq[b1..a2], &seq[a1..b1], &seq[b2..]].concat(), b2 - a1)
                            }
                            "inversion" => {
                                let (x, y) = pick_block(&mut rng, k);
                                ([&seq[..x], &revcomp(&seq[x..y])[..], &seq[y..]].concat(), y - x)
                            }
                            "join_same_sf" | "join_cross_sf" => {
                                let (x, y) = pick_block(&mut rng, k);
                                let seg = donor_seg(&mut rng, Some(class == "join_same_sf"), y - x)?;
                                ([&seq[..x], &seg[..], &seq[y..]].concat(), y - x)
                            }
                            "foreign_insertion" => {
                                let len = 2000 * k;
                                let seg = donor_seg(&mut rng, None, len)?;
                                let (lo, hi) = if u.is_empty() { (0, seq.len()) } else { (u[0].0, u[u.len() - 1].1) };
                                let p = lo + rng.below(hi - lo);
                                ([&seq[..p], &seg[..], &seq[p..]].concat(), len)
                            }
                            "indels" => {
                                // rate = dose x 1e-5 per bp inside the HOR span; 1-50 bp, half insertions
                                let (lo, hi) = if u.is_empty() { (0, seq.len()) } else { (u[0].0, u[u.len() - 1].1) };
                                let rate = k as f64 * 1e-5;
                                let n = ((hi - lo) as f64 * rate).round() as usize;
                                let mut pos: Vec<usize> = (0..n).map(|_| lo + rng.below(hi - lo)).collect();
                                pos.sort();
                                pos.dedup();
                                let mut out = Vec::with_capacity(seq.len());
                                let mut last = 0usize;
                                let mut changed = 0usize;
                                for p in pos {
                                    if p < last {
                                        continue;
                                    }
                                    out.extend_from_slice(&seq[last..p]);
                                    let l = 1 + rng.below(50);
                                    if rng.unit() < 0.5 {
                                        for _ in 0..l {
                                            out.push(b"ACGT"[rng.below(4)]);
                                        }
                                        last = p;
                                    } else {
                                        last = (p + l).min(seq.len());
                                    }
                                    changed += l;
                                }
                                out.extend_from_slice(&seq[last..]);
                                (out, changed)
                            }
                            "box_destroy" => {
                                // fraction = dose / 1000 of boxes lose a core base (C of TTCG, or G of CCCG)
                                let frac = k as f64 / 1000.0;
                                let mut out = seq.clone();
                                let mut changed = 0usize;
                                for (s0, _, fwd) in scan_cenpb(seq) {
                                    if rng.unit() < frac {
                                        let p = s0 as usize + if fwd { 3 } else { 4 };
                                        out[p] = if fwd { b'T' } else { b'A' };
                                        changed += 1;
                                    }
                                }
                                (out, changed)
                            }
                            _ => return Err("unknown class".into()),
                        })
                    })();
                    match edited {
                        Err(why) => {
                            writeln!(ew, "{id}\t{sname}\t{chrom}\t{class}\t{dose}\t{pl}\tNA\tNA\tNA\tNA\tNA\tNA\t{why}")?;
                        }
                        Ok((es, bp)) => {
                            let (h, acc) = histogram(&es);
                            let (js, kl) = match (to_dist(&h), &src_d) {
                                (Some(e), Some(s)) => {
                                    let m = metrics(&e, s, &[1e-12]);
                                    (format!("{}", m[8]), format!("{}", m[4]))
                                }
                                _ => ("NA".into(), "NA".into()),
                            };
                            let gaps_changed = l1_half(&h, &src_h);
                            writeln!(ew, "{id}\t{sname}\t{chrom}\t{class}\t{dose}\t{pl}\t{bp}\t{}\t{}\t{gaps_changed}\t{js}\t{kl}\tok", es.len(), acc.n_boxes)?;
                            for (g, c) in &h {
                                writeln!(hw, "{id}\t{g}\t{c}")?;
                            }
                            let status = if acc.n_boxes >= 50 { "ok" } else { "insufficient" };
                            writeln!(sw, "{id}\t1\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{status}", es.len(), acc.n_boxes, acc.n_fwd, acc.n_rev, acc.n_gaps, acc.n_overlap, acc.max_gap)?;
                        }
                    }
                }
            }
        }
        eprintln!("ckl perturb: {sname} done");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scanner_finds_both_strands() {
        let fwd = b"ATTCGAAAAAAACGGGA".to_vec();
        let rev = revcomp(&fwd);
        assert_eq!(scan_cenpb(&fwd), vec![(0, 17, true)]);
        assert_eq!(scan_cenpb(&rev), vec![(0, 17, false)]);
    }

    #[test]
    fn box_destroy_positions_break_the_match() {
        let mut fwd = b"ATTCGAAAAAAACGGGA".to_vec();
        fwd[3] = b'T';
        assert!(scan_cenpb(&fwd).is_empty());
        let mut rev = revcomp(b"ATTCGAAAAAAACGGGA");
        rev[4] = b'A';
        assert!(scan_cenpb(&rev).is_empty());
    }
}
