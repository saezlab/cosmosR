# Meta Prior Knowledge Network

Comprehensive Prior Knowledge Network (PKN), which combines signaling
and metabolic interaction networks. The network was constructed using
the Recon3D and STITCH metabolic networks as well as the signaling
network from OmniPath.

## Usage

``` r
data(meta_network)
```

## Format

An object of class “`tibble`” with 117065 rows (interactions) and three
variables:

- `source`:

  Source node, either metabolite or protein

- `interaction`:

  Type of interaction, 1 = Activation, -1 = Inhibition

- `target`:

  Target node, either metabolite or protein

## Source

The network is available in Omnipath:
<https://metapkn.omnipathdb.org/metapkn__20200122.txt>, the scripts used
for the build of the network are available under
<https://github.com/saezlab/Meta_PKN>.

## References

Dugourd, A., Kuppe, C. and Sciacovelli, M. et. al. (2021) *Molecular
Systems Biology*. **17**, e9730.

## Examples

``` r
data(meta_network)
```
