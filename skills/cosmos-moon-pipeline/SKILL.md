---
name: cosmos-moon-pipeline
description: Run, troubleshoot, and interpret the operational COSMOS+ MOON workflow in cosmosR. Use for high-level upstream/downstream input decisions via the repository PKN-semantic guide, cleaning and pruning signed PKNs, safe same-children compression, MOON scoring and TF-target coherence loops, score decompression, solution-network extraction, node-name translation, node-specific score explanations, and manuscript-style pathway-control analysis. Confirm a proposed input mapping with the user before finalizing it unless they explicitly request autonomous completion. Do not use for classic CARNIVAL optimization.
---

# COSMOS+ MOON Pipeline

Start with the decision gate below, then run the operational workflow. Treat every network as a mechanistic hypothesis. Preserve internal identifiers and the raw PKN; do not translate labels, claim causal proof, or discard provenance before the final display/export stage.

## Decide and confirm the input mapping

Before choosing MOON directionality, read the package checkout's `AGENTS.md` and `agent-docs/moon-data-pkn-mapping-principles.md`. Also read `agent-docs/moon-latent-factor-principles.md` when input scores derive from latent factors. If no checkout is available, obtain the same current files from `saezlab/cosmosR` before finalizing the mapping. Use these guides to distinguish an activity-compatible PKN input from a measurement that belongs only in a consistency check or annotation.

Propose, in concrete terms:

- the upstream and downstream input vectors, their node semantics, and their identifier domain;
- the network direction and the biological question it answers;
- each measurement retained as a consistency check or annotation rather than a MOON input;
- material assumptions, ambiguities, and alternative valid directions.

Unless the user explicitly asks for autonomous completion, ask for confirmation before treating this decision as final or running a direction-dependent MOON workflow. Restate and check a user-supplied mapping too: it may be a request rather than a biologically settled decision. A request to "run the pipeline" alone is not permission to silently choose the biological mapping. If the user asks for autonomous completion, record the selected mapping and rationale in the run ledger instead of asking.

## Use the workflow in this order

1. Define one identifier domain and make named, unique, finite numeric input vectors. Keep raw measurements, derived TF/kinase activities, and RNA evidence separate.
2. Build a minimal working PKN with exactly `source`, `interaction`, and `target`; clean it, expression-filter it if justified, and record edge/node/input coverage after every filter. Read [PKN preparation](references/pkn-preparation.md).
3. Prune for reachability from upstream inputs and observability from downstream inputs. Re-filter inputs after each prune. Treat the two radii as a permissive search space, not proof that every retained node lies on a seed-to-observation path.
4. Compress only after the uncompressed PKN and input coverage are stable. Preserve the compression mapping and validate that no protected input or TF-regulon node was renamed. Skip compression if that invariant fails.
5. Run `moon()` on the final scoring PKN, inspect the realized levels and upstream coverage, then run a bounded TF-target coherence loop when RNA evidence is appropriate. Read [MOON scoring](references/moon-scoring.md).
6. Keep the final scored PKN, its uncompressed counterpart, score table, compression mapping, thresholds, and convergence log together. Use one matching identifier domain for every downstream topology operation.
7. Decompress scores deliberately, audit a reduced subnetwork independently, and translate only copies used for presentation. Read [network interpretation](references/network-interpretation.md).
8. Use controller-anchored pathway analysis only as a separate ORA workflow with an explicit eligible-gene universe. Read [pathway control](references/pathway-control.md).

## Establish a run ledger

Create a small machine-readable record for each run. Include package commit/version, input construction, PKN source and cleanup counts, lost input IDs, `n_steps`, `n_layers`, statistic, thresholds, compression decision, coherence-loop iterations, and all output filenames. Retain both internal-ID and display-label SIF/ATT tables.

## Apply non-negotiable checks

- Require PKN signs to be numeric `-1` or `1`, with no self loops, missing endpoints, or unrelated metadata columns passed to `meta_network_cleanup()`.
- Require names of every input vector to be unique and present in the current PKN before calling the internal pruning helpers.
- Keep `interaction` as the PKN sign column through pruning, extraction, and score-network inspection. `moon()` accepts `sign` too, but several helpers do not.
- Use the PKN actually supplied to the last successful `moon()` call for score-level explanation. Use a coherently filtered, identifier-matched PKN for extraction; never fall back to the global raw PKN.
- Check reduced edges explicitly: `sign(score[source] * score[target]) == interaction`. Also check that every retained component has an intended upstream-to-level-0 path; neither reducer guarantees this in all cases.
- Treat score cutoffs as analysis choices, not universal biological significance cutoffs. Inspect distributions and perform sensitivity checks when a result matters.

## Resolve common requests

| Request | Read first | Required outcome |
| --- | --- | --- |
| “Prepare my PKN for MOON” | [PKN preparation](references/pkn-preparation.md) | A cleaned, coverage-audited, pruned PKN and matching inputs |
| “Run or debug MOON” | [MOON scoring](references/moon-scoring.md) | A scored table with level, convergence, and coverage checks |
| “Show why this node scored highly” | [network interpretation](references/network-interpretation.md) | A node-specific scoring subnetwork, not an unsupported causal claim |
| “Produce a compact network” | [network interpretation](references/network-interpretation.md) | An audited, score-consistent, input-connected SIF/ATT pair |
| “Make metabolite names readable” | [network interpretation](references/network-interpretation.md) | A separate display copy plus an immutable internal-ID table |
| “Which pathways might this controller affect?” | [pathway control](references/pathway-control.md) | A multiple-testing-aware controller-by-pathway result with explicit background |

## Distinguish the algorithm from the current package

Use the COSMOS+ manuscript’s Algorithms 1–7 as the conceptual workflow, but verify each requested operation against the local cosmosR source. The current package has known operational gaps: same-children compression can rename protected nodes, the RNA coherence helper does not implement protein-aware or pre-score checks, and no exported MOON pathway-control function exists. The references state safe workarounds and checks rather than assuming the manuscript pseudocode is automatically enforced.
