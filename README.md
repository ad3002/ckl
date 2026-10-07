# ckl

Centromeres as distributions of CENP-B box spacings: build the histograms, compare them, identify chromosomes and
haplotypes, and measure which assembly errors the comparison can and cannot see.

Written in Rust (single binary, one dependency: `clap`).

## Build

```bash
cargo build --release      # binary: target/release/ckl
cargo test --release
```

## Input

CENP-B boxes as BED (`chrom start end name score strand`), from any exact scanner of `NTTCGNNNNANNCGGGN` on both
strands. The distance between consecutive boxes is `start(i) - end(i-1)` (the GCP Model 1 convention of
Corda & Giunta 2025), so distances one alpha-satellite monomer apart peak at 171 - 17 = 154 bp.

## Commands

| command | does |
|---|---|
| `hist` | gap histograms per region (whole sequences, or a regions TSV such as active HOR arrays; gaps never cross intervals); writes standard / outlier / exception streams |
| `score` | divergences between histograms: KL in several zero-count policies (additive ε, renormalized, Laplace), symmetric KL, Jensen-Shannon, Hellinger, Euclidean, weighted Jaccard, Wasserstein-1. `--mode pairwise` (all pairs) or `--mode loio` (each query against per-chromosome references built without its own group) |
| `nn` | held-out nearest neighbour from any all-vs-all distance table (e.g. `mash dist`) |
| `kin` | within-chromosome nearest neighbours: parent-as-nearest in trios (with chance level and the child's other haplotype), and group-label accuracy |
| `refdep` | reference dependence: pool component histograms into consensus recipes and score assemblies against each, with the fraction of query mass in bins the reference lacks |
| `auroc` | AUROC of edited or erroneous arrays against validated ones, chromosome-normalized, plus detection at one fixed threshold |
| `perturb` | controlled edits of HOR arrays (collapse, duplication, reordering, inversion, joins, foreign insertion, indels, box destruction), rescanned and scored against the source |
| `censat-regions`, `extract`, `meta-set`, `xident` | helpers for annotation-based regions, per-region FASTA for k-mer baselines, and cross-species identification |

Meta tables (`--meta`) have the header `region_id individual group hap chrom role`. The `group` column is the
held-out unit (an individual, or a family when relatives are present).

## Example

```bash
ckl hist --boxes sample.cenpb.bed --regions active_hor.tsv --out sample
ckl score --mode loio --hist panel.hist.tsv --summary panel.summary.tsv --meta panel.meta.tsv --eps 1e-12 --out panel.loio.tsv
ckl kin --dist panel.pairwise.tsv --metric js_dist --meta panel.meta.tsv --trios trios.tsv --label js --out kin.js
```

## Notes on the score

With an additive pseudocount ε, every unit of query probability mass in a bin the reference lacks adds about
`-ln ε` (≈ 28 nats at ε = 1e-12). KL with such ε therefore counts rare distances: it is sensitive to sparse indels,
and a genome scores better against any reference that contains it. `refdep` reports that mass explicitly.

## License

MIT
