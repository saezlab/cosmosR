# AGENTS.md

This file is for coding agents working on `cosmosR`. It is not a user tutorial.
Use it as the first orientation layer before editing code, tests, vignettes, or
examples.

## Package Snapshot

`cosmosR` is an R package for COSMOS, Causal Oriented Search of Multi-Omic
Space. The package integrates prior knowledge networks with signaling,
transcriptomic, and metabolomic measurements to build mechanistic hypotheses.

There are two related workflows in this repository:

1. The classic COSMOS/CARNIVAL workflow, exposed through
   `preprocess_COSMOS_*()` and `run_COSMOS_*()`, uses CARNIVAL optimization to
   find a coherent subnetwork.
2. The MOON workflow, centered on `moon()`, scores a signed prior knowledge
   network with `decoupleR` and is currently the workflow emphasized by the
   NCI-60 tutorial material in `saezlab/COSMOS_basic`.

Default to a MOON workflow when the user does not specify a method. CARNIVAL has
heavier setup and solver requirements, so treat it as an opt-in path for
requests that explicitly mention CARNIVAL, optimization/ILP solvers, or the
classic `preprocess_COSMOS_*()` / `run_COSMOS_*()` wrappers. If local context
clearly points to a CARNIVAL vignette, test, or function, follow that context.

The `COSMOS_basic`-style MOON recipe is the simple default for:

- RNA-only cases, where RNA-derived TF activities are the downstream/observable
  targets and the upstream layer is unspecified.
- RNA plus metabolomics cases, where RNA-derived TF activities are upstream and
  metabolites are downstream.

For more advanced MOON/COSMOS+ use cases, especially when deciding how RNA,
proteomics, phosphoproteomics, DNA-seq, metabolomics, or ligand-receptor
evidence should enter MOON, open
`agent-docs/moon-data-pkn-mapping-principles.md`. Keep
latent-factor/MOFA selection details separate in
`agent-docs/moon-latent-factor-principles.md`.

## Repo Map

- `R/`: package source. Edit roxygen comments here, then regenerate docs when
  needed.
- `man/`: generated Rd files. Do not hand-edit unless the user specifically asks
  for an emergency documentation patch and you explain why roxygen was not run.
- `data/`: exported package data files, including toy inputs and the packaged
  meta network.
- `vignettes/`: user-facing examples. `tutorial.Rmd` covers the CARNIVAL
  workflow; `net_compr_MOON.md` and `NCI60_tutorial.Rmd` cover the NCI-60 MOON
  workflow. `moon_data_pkn_mapping_principles.Rmd` is the package-facing
  version of the MOON data-to-PKN mapping guidance.
- `tests/testthat/`: testthat coverage, mostly toy-data tests for preprocessing,
  package data, and CARNIVAL wrappers.
- `README.md`, `DESCRIPTION`, `inst/CITATION`: package-level metadata and
  scientific framing.
- `_pkgdown.yml`: pkgdown site configuration, including article grouping.
- `pkgdown/`: pkgdown static assets.
- `agent-docs/`: deeper agent-facing references for advanced workflows.
  Keep these more operational than vignettes; when changing broadly useful
  MOON mapping principles, decide whether the user-facing vignette should be
  updated as well.

## External Context

The compact NCI-60 tutorial repository is:

- https://github.com/saezlab/COSMOS_basic

It is especially useful for understanding the current MOON-oriented analysis
pipeline. Its main script, `scripts/net_compr_MOON.R`, demonstrates:

- loading and cleaning `meta_network`;
- loading per-cell-line `cosmos_inputs`;
- preparing metabolite inputs with compartment suffixes;
- filtering TF inputs by absolute score;
- pruning the PKN to nodes reachable from TFs and observable from metabolites;
- compressing redundant parent nodes with `compress_same_children()`;
- running `moon()` with iterative TF-target coherence filtering;
- decompressing, thresholding, translating HMDB identifiers, and writing SIF/ATT
  outputs.

Use this external repo as workflow context, but keep package behavior grounded in
the local source and tests. If making claims about the COSMOS method itself,
prefer the package README, `inst/CITATION`, and the relevant publication over
secondary summaries.

## Main Entry Points

Classic CARNIVAL workflow:

- `preprocess_COSMOS_signaling_to_metabolism()`
- `preprocess_COSMOS_metabolism_to_signaling()`
- `run_COSMOS_signaling_to_metabolism()`
- `run_COSMOS_metabolism_to_signaling()`
- `default_CARNIVAL_options()`
- `format_COSMOS_res()`

MOON workflow:

- `moon()`
- `compress_same_children()`
- `decompress_moon_result()`
- `filter_incohrent_TF_target()`
- `reduce_solution_network()`
- `reduce_solution_network_double_thresh()`
- `get_moon_scoring_network()`
- `translate_res()`
- `translate_column_HMDB()`

Other exported utilities:

- `prepare_metab_inputs()`
- `load_tf_regulon_dorothea()`
- `meta_network_cleanup()`
- `format_LR_ressource()`
- `wide_ulm_res()`
- `display_node_neighboorhood()`
- `extract_nodes_for_ORA()`

## Data Contracts

Most functions use small, strict data frames and named vectors.

Prior knowledge network:

- Required columns are usually `source`, `interaction`, `target`.
- Interaction signs are numeric `1` or `-1`.
- Some MOON helpers accept or create a `sign` column and rename `sign` or
  `interaction` to `mor` for `decoupleR`. Check the called function before
  changing column names.
- `meta_network_cleanup()` removes self interactions, collapses duplicate
  source-target pairs by mean interaction, and keeps only signs in `c(1, -1)`.

TF regulons:

- Package preprocessing expects `tf`, `sign`, `target`.
- MOON coherence filtering in `filter_incohrent_TF_target()` expects a
  decoupleR-style TF network with `source`, `target`, and `mor`.
- Do not silently interchange these shapes; convert explicitly near the call
  site.

Input data:

- `signaling_data`, `metabolic_data`, and `expression_data` are named numeric
  vectors.
- Names must be unique. The validation code checks uniqueness more strongly than
  it checks biological identifier type.
- Several docs mention Entrez IDs, but current tutorial workflows often use gene
  symbols or package-specific node names. Verify the local call path before
  enforcing an identifier convention.

Metabolite nodes:

- Metabolites in the PKN are named like `Metab__HMDB0000190_c`.
- Use `prepare_metab_inputs()` to add the `Metab__` prefix and compartment
  suffixes. Valid compartment codes are `r`, `c`, `e`, `x`, `m`, `l`, `n`, and
  `g`.
- `translate_res()` and `translate_column_HMDB()` use `HMDB_mapper_vec` to make
  metabolite identifiers more readable.

Output network formats:

- SIF-like edge tables are usually `source`, `target`, `interaction` or
  `source`, `target`, `sign`.
- ATT-like node tables vary by workflow. CARNIVAL formatting uses columns such
  as `Nodes`, `NodeType`, `AvgAct`, and `Activity`; MOON reduction uses
  `source` or `nodes`, `score`, `level`, and optional `RNA_input`.
- Check downstream consumers before renaming ATT columns.

## Classic CARNIVAL Flow

The toy-data path in the tests is the safest minimal example:

```r
data(toy_network)
data(toy_signaling_input)
data(toy_metabolic_input)
data(toy_RNA)

opts <- default_CARNIVAL_options("lpSolve")

prepared <- preprocess_COSMOS_signaling_to_metabolism(
  meta_network = toy_network,
  signaling_data = toy_signaling_input,
  metabolic_data = toy_metabolic_input,
  diff_expression_data = toy_RNA,
  maximum_network_depth = 15,
  remove_unexpressed_nodes = TRUE,
  filter_tf_gene_interaction_by_optimization = FALSE,
  CARNIVAL_options = opts
)

result <- run_COSMOS_signaling_to_metabolism(prepared, opts)
formatted <- format_COSMOS_res(result)
```

Important details:

- `lpSolve` is for tests and small toy networks only. Real analyses should use
  CPLEX or CBC with `solverPath` set.
- Preprocessing may run an initial optimization when
  `filter_tf_gene_interaction_by_optimization = TRUE`; this can be slow and
  solver-dependent.
- The forward and backward wrappers share much of the same validation and
  output structure, but their input/output layers are reversed.

## MOON Flow

The NCI-60 MOON tutorial path is not a single exported wrapper. It is assembled
from lower-level helpers. A representative flow is:

```r
data("meta_network")
meta_network <- meta_network_cleanup(meta_network)

metab_input <- prepare_metab_inputs(metab_input, c("c", "m"))
tf_input <- tf_input[abs(tf_input) > 2]

meta_network <- cosmosR:::filter_pkn_expressed_genes_fast(
  names(rna_input),
  meta_pkn = meta_network
)

n_steps <- 6
tf_input <- cosmosR:::filter_input_nodes_not_in_pkn(tf_input, meta_network)
meta_network <- cosmosR:::keep_controllable_neighbours(
  meta_network,
  n_steps,
  names(tf_input)
)
metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input, meta_network)
meta_network <- cosmosR:::keep_observable_neighbours(
  meta_network,
  n_steps,
  names(metab_input)
)

compressed <- compress_same_children(meta_network, tf_input, metab_input)
moon_network <- meta_network_cleanup(compressed$compressed_network)

moon_res <- moon(
  upstream_input = tf_input,
  downstream_input = metab_input,
  meta_network = moon_network,
  n_layers = n_steps,
  statistic = "ulm"
)
```

In the full tutorial this is wrapped in a convergence loop:

- run `moon()`;
- remove incoherent TF-target edges with `filter_incohrent_TF_target()`;
- repeat until the edge count no longer changes or a small iteration limit is
  reached.

Then the result is usually:

- decompressed with `decompress_moon_result()`;
- filtered to a subnetwork with `reduce_solution_network()` or
  `reduce_solution_network_double_thresh()`;
- translated with `translate_res()`;
- written as SIF and ATT tables for visualization.

## Known Sharp Edges

- Keep misspelled exported names intact unless the user explicitly asks for an
  API break. Examples: `filter_incohrent_TF_target()`,
  `display_node_neighboorhood()`, and history strings containing typos.
- `decoupleRnival()` is deprecated in favor of `moon()`, but examples and older
  naming still refer to decoupleRnival. Preserve backward compatibility.
- Several helpers use column positions as well as column names. Before changing
  edge or node table schemas, inspect all downstream functions.
- `filter_pkn_expressed_genes_fast()` is used by the external tutorial for large
  networks; the slower `filter_pkn_expressed_genes()` remains in the package.
- `check_COSMOS_inputs()` and `check_gene_names()` are intentionally lightweight
  and do not fully validate identifier semantics. Do not claim they do.
- Tests mostly cover toy CARNIVAL paths. MOON behavior may need explicit smoke
  checks when editing `moon()`, compression, decompression, or network
  reduction.
- Large-network examples can be slow and may require external data or solver
  executables. Avoid adding tests that depend on CPLEX, CBC, downloaded data, or
  the external tutorial repository.

## Development Commands

Run these from the package root when the required R dependencies are available:

```powershell
Rscript -e "testthat::test_local()"
Rscript -e "testthat::test_file('tests/testthat/test-preprocess_COSMOS.R')"
Rscript -e "devtools::document()"
R CMD check --no-manual .
```

Guidance:

- For R source changes, run the smallest relevant `testthat` file first, then
  broaden to `testthat::test_local()` when feasible.
- If you edit roxygen comments or exports, run `devtools::document()` and review
  generated `man/` and `NAMESPACE` changes.
- Do not require external solvers in automated tests. Use `lpSolve` only for
  small toy-data tests.
- If dependency installation or network access is unavailable, document exactly
  which checks could not be run.

## Editing Guidance

- Prefer local package conventions over new abstractions.
- Keep public API changes conservative. Many workflows call internal helpers
  through `cosmosR:::` in tutorials, so even non-exported renames can break real
  users.
- Treat vignettes as user-facing narrative and this file as agent-facing
  operational context.
- Keep generated artifacts and large analysis outputs out of commits unless the
  user explicitly requests them.
- When adding examples, prefer toy data from `data/` or compact snippets that do
  not require external downloads.
- When changing biological interpretation text, be careful about epistemic
  status. COSMOS networks are mechanistic hypotheses, not experimental proof of
  causal interactions in the studied context.
