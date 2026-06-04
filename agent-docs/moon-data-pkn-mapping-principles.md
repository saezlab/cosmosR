# MOON Data-to-PKN Mapping Principles

This is the main advanced agent-facing guide for deciding how data should enter
a MOON/COSMOS+ network-scoring workflow. The central question is:

> What biological quantity does this measurement represent, and does that match
> the type of causal interaction encoded by the relevant PKN edges?

Do not map an omic measurement to a node just because the feature name matches.
Only use it as a MOON upstream/downstream node activity when the measurement is
a reasonable proxy for the node state propagated by the signed PKN.

This first draft is based on the COSMOS+ preprint, especially sections 2.1,
2.3, and 4.4 of:

`C:\Users\dugourd\Downloads\2024.07.15.603538v3.full (1).pdf`

It also uses the Factor_COSMOS tutorial code, especially the
`Filtering_RNA_weights` section of:

https://github.com/saezlab/Factor_COSMOS/blob/main/MOFA_to_COSMOS.Rmd

See `agent-docs/moon-latent-factor-principles.md` for the separate question of
how to select and orient latent factors before applying these mapping rules.

## Core Rule

The PKN contains signed, directed mechanistic hypotheses. Its protein/protein
and protein/metabolite edges usually represent activity propagation or causal
regulation, not mere abundance correlation.

Therefore:

- Use footprint-derived activities when the node state is an activity.
- Use direct measurements only when the PKN node semantics match the measured
  quantity.
- Use direct omics measurements as consistency filters or annotations when they
  are informative but do not directly represent propagated activity.
- Keep every transformation explicit: raw measurement, derived activity,
  consistency check, network pruning, and post-hoc annotation are different
  roles.

## RNA

RNA abundance or RNA factor weights should usually not be mapped directly onto
protein signaling nodes as activity.

Preferred roles:

- Estimate transcription factor activity from RNA using TF-target footprints,
  for example CollecTRI with `decoupleR::run_ulm()`.
- Use RNA as a TF-target consistency check: after MOON scores TFs, remove
  TF-target edges whose predicted effect disagrees with the observed transcript
  sign.
- When matched total proteomics is available, build the TF-target consistency
  vector from a protein-gated RNA readout: strong RNA support is required, and
  matched weak or contradictory protein support means the target should not be
  treated as a functional downstream readout of the TF.
- Use RNA to define expressed genes and remove unexpressed protein/gene nodes
  from the PKN when that is appropriate for the analysis.
- Use RNA weights as feature annotations for interpretation, especially when
  comparing a scored node with its transcript-level support.

Avoid:

- Treating RNA abundance of a kinase, receptor, enzyme, or TF as direct protein
  activity.
- Using RNA alone to assert post-translational signaling activity.

Reasoning:

RNA is well matched to transcriptional regulation because TF-target edges
predict transcript changes. RNA is less well matched to protein-protein
signaling edges, which usually encode protein activity, phosphorylation,
complex formation, localization, or regulation.

### Protein-gated RNA for TF-target filtering

The Factor_COSMOS `Filtering_RNA_weights` section provides a concrete pattern:
RNA feature weights are thresholded first; then, when a matched protein weight is
available for the same gene, weak protein support also zeroes the RNA readout
used later in the MOON TF-target consistency loop. Conceptually, this means:

- RNA says whether a TF-target transcript appears regulated in the right
  direction.
- Matched proteomics asks whether that transcript-level regulation plausibly
  translates into a protein-level downstream functional readout.
- If RNA is coherent with the upstream TF but the matched protein abundance does
  not follow, the TF-target edge should not be used as evidence for downstream
  functional regulation.

Implementation note: the current package helper `filter_incohrent_TF_target()`
removes edges with opposite signs via `sign(TF_score * RNA_input * mor) < 0`.
If unsupported targets are encoded as `RNA_input = 0`, that condition alone does
not remove them. For protein-aware workflows, agents should either pre-prune
unsupported TF-target edges or use/update a coherence filter that removes both
sign-incoherent targets and targets with no protein-supported functional
readout when protein evidence is available.

## TF Activities

TF activities are derived variables, not raw RNA measurements.

Preferred roles:

- Use TF footprint scores as MOON downstream/observable targets when asking
  what upstream regulators could explain a transcriptional program.
- Use TF footprint scores as upstream inputs when asking what downstream
  metabolism or ligand/receptor programs could be consistent with TF activity.
- In simple RNA-only cases, TF activities are the downstream/observable layer
  and the upstream layer may be unspecified.
- In simple RNA plus metabolomics cases, TF activities are usually upstream and
  metabolites downstream, following the `COSMOS_basic` pattern.

Check:

- The sign of the TF score should be interpreted as TF activity, not TF mRNA
  abundance.
- If TF-target coherence filtering is used, pass the raw or factor-weight RNA
  vector separately as `RNA_input`; do not replace it with TF activity scores.

## Total Proteomics

Total protein abundance is not the same as protein activity.

Preferred roles:

- Use total protein abundance or protein factor weights as corroborating
  evidence that a protein is present or associated with the signal being
  interpreted.
- Use total proteomics to gate transcript-supported TF-target regulation when
  the analysis requires evidence that RNA changes translate into protein-level
  downstream readouts.
- Use total protein weights as post-hoc node annotations, so a high MOON score
  can be compared with measured protein abundance support.

Avoid:

- Mapping total protein abundance directly onto kinase, receptor, or signaling
  nodes as if abundance were activity.
- Treating total proteomics as a substitute for phosphorylation, kinase
  substrate footprint, or other activity-oriented evidence.

Section 2.3/NCI60 example:

- The NCI60 analysis had total proteomics, not phosphoproteomics.
- Protein weights were used to assess RNA-protein factor consistency and to
  gate RNA-supported TF-target readouts in the `Filtering_RNA_weights` logic.
- Protein weights were not used as generic protein activity constraints on
  signaling edges.
- High protein or transcript weights were used as corroborating evidence for
  nodes with high MOON scores, not as proof of activity by themselves.

## Phosphoproteomics

Phosphoproteomics is better matched to signaling activity than total
proteomics, but it still usually needs a footprint step.

Preferred roles:

- Estimate kinase or phosphatase activities from phosphosite or substrate
  measurements using a kinase-substrate prior.
- Use those derived kinase/phosphatase activity scores as MOON inputs when they
  match the signaling layer being scored.
- Use phosphosite measurements directly only when the PKN explicitly represents
  that modified site or a node whose state is defined by that modification.

Avoid:

- Collapsing phosphoproteomic evidence to total protein abundance semantics.
- Treating every measured phosphosite as a generic activity score for the whole
  protein without a stated aggregation or footprint model.

## DNA Sequencing

DNA sequencing identifies genomic alterations. These alterations are usually
candidate causal perturbations, not direct activity measurements.

Preferred roles:

- Use high-confidence, directionally interpretable alterations as upstream
  perturbation candidates, for example known activating mutations,
  homozygous deletions, truncating loss-of-function events, amplifications with
  expression support, or fusions with known direction.
- Use DNA events as annotations when direction is ambiguous, evidence is weak,
  or the altered gene is not represented in the PKN with compatible semantics.
- Use RNA-seq alongside DNA-seq to infer downstream TF activities and to check
  whether altered genes are expressed in the relevant context.
- Use DNA-derived lesion rules to prune or disable edges whose functional
  assumptions are broken by the alteration.

Avoid:

- Mapping "gene has a mutation" directly to `+1` activity.
- Treating all missense mutations, variants of unknown significance, or
  passenger-like events as causal PKN inputs.
- Treating copy-number gain as activity unless dosage is plausible and
  supported by expression or other context.

Example DNA event roles:

- Known activating `BRAF V600E`: upstream positive candidate if `BRAF` is
  expressed and the relevant PKN node represents BRAF signaling activity.
- Homozygous `PTEN` deletion or strong loss-of-function: negative perturbation
  of PTEN function; downstream consequences may correspond to release of
  inhibition rather than a simple `PTEN = -1` edge readout.
- `ERBB2` amplification: positive upstream candidate when dosage is reflected in
  RNA/protein and the receptor is relevant to the modeled network.
- Ambiguous `TP53` missense: do not force into a sign without deciding whether
  it is loss-of-function, dominant-negative, gain-of-function, or uncertain.
- Low-confidence VUS in a lowly expressed gene: annotation only by default.

### Early stop codons and other likely loss-of-function lesions

An early stop codon should usually be treated as a node-level lesion, not as an
ordinary upstream regulatory signal.

For a PKN edge `A -> B`, where `B` has a likely functional early stop codon:

- Incoming `A -> B` edges should usually not be used to explain loss of B
  activity. The loss is genetically imposed, not necessarily caused by upstream
  regulation.
- Outgoing `B -> C` edges should usually be retained by default and interpreted
  under a forced negative activity for B. In a signed activity-propagation
  network, loss of B function is itself a causal state that should affect
  downstream targets according to the sign of B's outgoing edges.
- Disable outgoing `B -> C` edges only when the edge semantics require a
  retained molecular function that the lesion clearly destroys and this cannot
  be represented as decreased B activity.
- B may itself become a negative upstream perturbation or candidate driver if
  loss of B function is directionally meaningful for the question.

Do not apply this mechanically. Check:

- allele status: heterozygous, homozygous, LOH, clonal fraction;
- nonsense-mediated decay and RNA expression of the mutant gene;
- protein/domain position of the stop codon;
- whether the relevant functional domains and interaction interfaces are lost
  or retained;
- whether the truncated product is unstable, partially functional,
  dominant-negative, or gain-of-function;
- whether the PKN edge represents protein activity, transcript regulation,
  complex membership, or another mechanism.

Default reasoning:

```text
early stop in B -> likely loss of functional B protein
incoming A -> B -> disable as an explanation for B loss
force B activity negative -> retain outgoing B -> C edges by default
outgoing B -> C -> disable only if the edge requires retained function not captured by activity loss
B itself -> possible negative genetic perturbation/candidate driver
```

For DNA+RNA workflows, use RNA-seq to decide whether the lesion is expressed
and to measure downstream consequences. The DNA event supplies a candidate
perturbation; RNA-derived TF activities and TF-target consistency provide the
functional readouts.

## Metabolomics

Metabolite measurements can often be used more directly as downstream MOON
inputs because metabolite nodes in the COSMOS PKN represent metabolite states.

Preferred roles:

- Map metabolite identifiers to PKN metabolite nodes and compartment suffixes
  with `prepare_metab_inputs()`.
- Use metabolite abundances, differential statistics, or factor weights as
  downstream inputs when the question is to infer upstream regulators.
- Keep compartment assumptions explicit; duplicating a metabolite into multiple
  compartments is an analysis choice, not a neutral mapping.

Check:

- Whether the metabolite identifier is HMDB, BiGG, or another namespace.
- Whether the measured metabolite is compatible with the compartments included
  in the PKN.
- Whether the sign represents abundance, differential abundance, or a latent
  factor weight.

## Ligands And Receptors

Ligand-receptor information should be represented according to what is scored.

Preferred roles:

- Estimate LR-pair scores from RNA and/or protein features when the question is
  whether ligand/receptor transcripts or proteins are co-regulated.
- Use receptor or ligand scores as upstream or downstream MOON inputs only
  after defining what the score represents.
- For networks connecting receptors to TFs/metabolites, receptors can be
  upstream candidates and TF activities/metabolites can be downstream targets.
- For TF-to-ligand networks, TF activities can be upstream and ligand scores can
  be downstream, often with a short path length when only direct TF-target
  ligand regulation is intended.

Avoid:

- Treating a ligand-receptor pair score as identical to receptor activity.
- Treating ligand/receptor RNA abundance as proof of extracellular signaling
  activity.

## Derived Activity Scores

Footprint scores are the safest way to bridge between omics measurements and
activity-like PKN nodes.

Common derived scores:

- TF activity from RNA and TF-target regulons.
- Kinase/phosphatase activity from phosphoproteomics and kinase-substrate
  resources.
- LR-pair scores from ligand/receptor RNA and/or protein feature co-regulation.

Agent rule:

If a PKN node represents an activity and the data type measures abundance, first
ask whether a footprint or other transformation can estimate the activity. If
not, keep the abundance as annotation or filtering evidence rather than as a
MOON activity input.

## Worked Reasoning Example: Adding A New Omic Layer

Use this example as a template when a dataset contains an omic layer that is not
already covered by the standard RNA, phosphoproteomics, metabolomics, or LR
recipes.

Example question:

> We have RNA and total proteomics. Can total proteomics improve TF-target
> filtering?

Reasoning:

1. Identify the PKN statement.
   A TF-target edge says that TF activity can regulate a target gene transcript,
   with a sign/mode of regulation.

2. Identify what RNA measures.
   RNA measures the target transcript. This is directly relevant to TF-target
   edges, because those edges predict transcript-level changes.

3. Identify what total proteomics measures.
   Total proteomics measures the abundance of the encoded protein. It does not
   measure TF activity and does not directly measure whether the TF bound or
   regulated the promoter. However, it can test whether the transcript-level
   signal plausibly reaches a downstream protein-level functional readout.

4. Decide each data type's role.
   RNA is the primary target readout for TF-target sign coherence. Total
   proteomics is a functional-readout gate when matched protein evidence exists.

5. Define the transformation.
   Build a protein-gated RNA vector:
   - If RNA is below the selected threshold, set target readout to `0`.
   - If RNA is above threshold and no matched protein measurement exists, keep
     the RNA readout, but mark protein support as unknown.
   - If RNA is above threshold and matched protein is also above threshold in a
     compatible direction, keep the RNA readout.
   - If RNA is above threshold but matched protein is weak or contradictory, set
     target readout to `0` or otherwise flag the TF-target edge as unsupported.

6. Define the edge rule.
   Remove TF-target edges when either:
   - the TF score, target readout, and TF-target sign are incoherent; or
   - the target readout is `0` because the target lacks sufficient functional
     readout support under the chosen RNA/protein rule.

7. State the interpretation.
   Remaining TF-target edges are not experimentally proven. They are edges for
   which the inferred TF activity, transcript response, and available
   protein-level support are mutually compatible with the PKN sign.

Generalized pattern:

```text
data type -> measured quantity -> matching PKN semantics -> permitted role

RNA -> transcript abundance/weight -> TF-target transcript regulation
    -> footprint for TF activity; TF-target consistency/gating

total proteomics -> protein abundance/weight -> protein presence or abundance
    -> protein-level support/gating/annotation, not generic activity

phosphoproteomics -> modified peptide/site abundance -> signaling modification
    -> kinase/phosphatase footprints; direct node input only for explicit site nodes

DNA-seq -> genomic alteration -> genetic perturbation of gene/protein function
    -> upstream candidate, lesion-based edge pruning, or annotation depending on direction/evidence

metabolomics -> metabolite abundance/weight -> metabolite node state
    -> often direct downstream input after identifier and compartment mapping
```

For a new omic layer, repeat the same reasoning before coding:

1. What exact molecular state is measured?
2. Which PKN node or edge type claims something about that state?
3. Is the measurement a direct proxy for propagated activity, a footprint
   substrate, a functional-readout gate, an expression/presence filter, or only
   an annotation?
4. What value means "unsupported" and should remove or ignore an edge?
5. What value means "unknown" and should not be treated as negative evidence?
6. Does the current helper function distinguish unsupported from unknown? If
   not, add that distinction before using the data in a filtering loop.

## Network Direction Choices

Define upstream and downstream layers by matching data semantics to the
question and the PKN, not by omic type alone.

Examples:

- RNA-only perturbation signatures: downstream = TF activities inferred from
  RNA; upstream may be ligands, receptors, drugs, or unspecified PKN sources.
- RNA plus metabolomics: upstream = TF activities; downstream = metabolite
  abundances or metabolite weights.
- Receptor-to-TF/metabolite factor analysis: upstream = receptor or LR-derived
  receptor scores; downstream = TF activity scores plus metabolite weights.
- TF-to-ligand analysis: upstream = TF activities; downstream = ligand scores or
  ligand features, with direct TF-target consistency checks.
- Phosphoproteomics plus RNA: upstream or intermediate = kinase/phosphatase
  activities inferred from phosphosites; downstream = TF activities from RNA.
- DNA-seq plus RNA: upstream candidates = signed, directionally interpretable
  genomic perturbations; downstream = TF activities inferred from RNA, with RNA
  expression used to support or reject DNA-derived lesion assumptions.

When in doubt, write down:

1. What node state does MOON score represent here?
2. Which measurements support that state directly?
3. Which measurements only support consistency, expression, or annotation?
4. Which PKN edge types connect the selected layers?

## Consistency Checks

Use consistency checks to exploit data that is informative but not a direct
activity proxy.

Important checks:

- TF-target coherence: remove TF-target edges where
  `sign(TF_score * RNA_input * mor) < 0`.
- Protein-aware TF-target functionality: when total proteomics is available for
  a target gene, remove or ignore the TF-target edge if the RNA change does not
  have matching protein-level support. This is a stricter check than transcript
  sign coherence alone.
- RNA/protein support: if total protein data is available, use it to assess
  whether RNA-level feature weights are also supported at protein level before
  emphasizing the node in interpretation.
- DNA lesion support: for likely loss-of-function events, treat the altered node
  as a forced negative activity by default and propagate that state through its
  outgoing edges. Disable incoming edges as explanations for the lesion. Disable
  outgoing edges only when the specific edge mechanism requires retained
  function that cannot be modeled as activity loss.
- Edge sign coherence: when reducing a MOON network, keep edges whose signs are
  consistent with the signs of the connected node scores.
- Reachability/observability: prune the PKN to nodes reachable from upstream
  inputs and observable from downstream inputs within the selected number of
  steps.

Do not overinterpret:

- A consistency check removes implausible edges under the chosen assumptions; it
  does not validate the remaining edges experimentally.
- A high MOON score means a node is coherent with the downstream layer and PKN
  signs; it is a mechanistic hypothesis, not direct proof.

## Section 2.3 Heuristics

These are paper-specific starting points, not package defaults:

- In the NCI60 factor-4 workflow, RNA feature weights below a visually selected
  threshold were set to zero; for genes with matched protein weights, weak
  protein support also led to zeroing before downstream mechanistic
  interpretation.
- Metabolite weights with small absolute values were also set to zero.
- Thresholds were selected from visual inspection of the factor-weight
  distributions.
- MOON was run once from receptors to downstream TF activities and metabolites,
  then again from TFs to downstream ligands with a one-step limit.
- The two networks were combined to connect ligands, receptors, TFs, and
  metabolites.
- Absolute MOON score thresholds around `1.5` were used for reduced networks
  and pathway-control analysis.

When reusing these heuristics, state that they are analysis choices and adjust
them to the data distribution, PKN coverage, and biological question.

## Agent Checklist

Before wiring data into `moon()`:

1. List every available data type and whether it measures abundance, activity,
   modification, interaction, or a latent weight.
2. List the PKN edge types involved in the planned analysis.
3. Decide which data become direct MOON inputs and which become footprints,
   filters, consistency checks, or annotations.
4. For RNA, derive TF activity and keep RNA separately for TF-target coherence.
5. For proteomics, distinguish total proteomics from phosphoproteomics.
6. For DNA-seq, classify each alteration by likely functional direction before
   using it as an upstream candidate or edge-pruning rule.
7. For metabolites, map identifiers and compartments explicitly.
8. Define upstream/downstream layers in terms of PKN semantics.
9. Keep output interpretation at the level of mechanistic hypotheses.
