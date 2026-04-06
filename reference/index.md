# Package index

## BLASE Object

- [`BlaseData-class`](https://andrewmccluskey-uog.github.io/blase/reference/BlaseData-class.md)
  [`BlaseData`](https://andrewmccluskey-uog.github.io/blase/reference/BlaseData-class.md)
  : Blase Data Object
- [`as.BlaseData()`](https://andrewmccluskey-uog.github.io/blase/reference/as.BlaseData.md)
  : Conversion to BlaseData
- [`assign_pseudotime_bins()`](https://andrewmccluskey-uog.github.io/blase/reference/assign_pseudotime_bins.md)
  : Assign Pseudotime Bins to a source object's metadata
- [`genes()`](https://andrewmccluskey-uog.github.io/blase/reference/genes-getter.md)
  : Get genes of a BLASE Data object.
- [`` `genes<-`() ``](https://andrewmccluskey-uog.github.io/blase/reference/genes-setter.md)
  : Set genes of a BLASE Data object.
- [`show(`*`<BlaseData>`*`)`](https://andrewmccluskey-uog.github.io/blase/reference/show_blase_object.md)
  : Show an BlaseData object

## BLASE Mapping Result Object

- [`MappingResult()`](https://andrewmccluskey-uog.github.io/blase/reference/MappingResult.md)
  : Blase Mapping Result
- [`` `bulk_name<-`() ``](https://andrewmccluskey-uog.github.io/blase/reference/bulk_name-setter.md)
  : Set name of bulk of a BLASE Mapping Results object.
- [`best_bin()`](https://andrewmccluskey-uog.github.io/blase/reference/mapping-result-best-bin-getter.md)
  : Get best bin of a BLASE Mapping Results object.
- [`best_correlation()`](https://andrewmccluskey-uog.github.io/blase/reference/mapping-result-best-correlation-getter.md)
  : Get best correlation of a BLASE Mapping Results object.
- [`bootstrap_iterations()`](https://andrewmccluskey-uog.github.io/blase/reference/mapping-result-bootstrap-iterations-getter.md)
  : Get the number of bootstrap iterations performed for a BLASE Mapping
  Results object.
- [`bulk_name()`](https://andrewmccluskey-uog.github.io/blase/reference/mapping-result-bulk-name-getter.md)
  : Get name of bulk of a BLASE Mapping Results object.
- [`mapping_history()`](https://andrewmccluskey-uog.github.io/blase/reference/mapping-result-history-getter.md)
  : Get the mapping history for a BLASE Mapping Results object.
- [`metric()`](https://andrewmccluskey-uog.github.io/blase/reference/mapping-result-metric-getter.md)
  : Get the mapping history for a BLASE Mapping Results object.
- [`strong_mapping()`](https://andrewmccluskey-uog.github.io/blase/reference/mapping-result-strong-mapping-getter.md)
  : Get if the result is strong for a BLASE Mapping Results object.
- [`top_2_distance()`](https://andrewmccluskey-uog.github.io/blase/reference/mapping-result-top-2-distance-getter.md)
  : Get the difference in correlation between the top 2 most correlated
  bins for a BLASE Mapping Results object.
- [`show(`*`<MappingResult>`*`)`](https://andrewmccluskey-uog.github.io/blase/reference/show-MappingResult-method.md)
  : Show an MappingResult object

## Mapping Bulks

- [`map_all_best_bins()`](https://andrewmccluskey-uog.github.io/blase/reference/map_all_best_bins.md)
  : Map many bulk samples in the same dataframe
- [`map_best_bin()`](https://andrewmccluskey-uog.github.io/blase/reference/map_best_bin.md)
  : Map the best matching SC bin for a bulk sample

### Plots

- [`plot_bin_population()`](https://andrewmccluskey-uog.github.io/blase/reference/plot_bin_population.md)
  : Plot the populations of a bin
- [`plot_mapping_result()`](https://andrewmccluskey-uog.github.io/blase/reference/plot_mapping_result.md)
  : Plot a summary of the mapping result
- [`plot_mapping_result_corr()`](https://andrewmccluskey-uog.github.io/blase/reference/plot_mapping_result_corr.md)
  : Plot a mapping result's correlation
- [`plot_mapping_result_heatmap()`](https://andrewmccluskey-uog.github.io/blase/reference/plot_mapping_result_heatmap.md)
  : Plot a mapping result heatmap

## Parameter Tuning

- [`evaluate_parameters()`](https://andrewmccluskey-uog.github.io/blase/reference/evaluate_parameters.md)
  : Evaluate n_bins and n_genes for bin mapping
- [`evaluate_top_n_genes()`](https://andrewmccluskey-uog.github.io/blase/reference/evaluate_top_n_genes.md)
  : Evaluate Top Genes
- [`find_best_params()`](https://andrewmccluskey-uog.github.io/blase/reference/find_best_params.md)
  : Identify the Best Parameters For Your Dataset
- [`plot_find_best_params_results()`](https://andrewmccluskey-uog.github.io/blase/reference/plot_find_best_params_results.md)
  : Plot the results of the search for good parameters

## Gene Selection

- [`calculate_gene_peakedness()`](https://andrewmccluskey-uog.github.io/blase/reference/calculate_gene_peakedness.md)
  : calculate_gene_peakedness
- [`gene_peakedness_spread_selection()`](https://andrewmccluskey-uog.github.io/blase/reference/gene_peakedness_spread_selection.md)
  : Gene Peakedness Spread Selection
- [`get_top_n_genes()`](https://andrewmccluskey-uog.github.io/blase/reference/get_top_n_genes.md)
  : Get Top Genes From An AssociationTestResult
- [`plot_gene_peakedness()`](https://andrewmccluskey-uog.github.io/blase/reference/plot_gene_peakedness.md)
  : plot_gene_peakedness
- [`smooth_gene()`](https://andrewmccluskey-uog.github.io/blase/reference/smooth_gene.md)
  : smooth_gene

## Annotation of scRNA-seq

- [`annotate_sce()`](https://andrewmccluskey-uog.github.io/blase/reference/annotate_sce.md)
  : Annotate a SCE with BLASE Mappings

## Utilities

- [`get_bins_as_bulk()`](https://andrewmccluskey-uog.github.io/blase/reference/get_bins_as_bulk.md)
  : Get a pseudobulk of bins with at least 2 replicates

## Data

- [`MCA_PF_SCE`](https://andrewmccluskey-uog.github.io/blase/reference/MCA_PF_SCE.md)
  : Malaria Cell Atlas Plasmodium falciparum for BLASE Vignette
- [`painter_microarray`](https://andrewmccluskey-uog.github.io/blase/reference/painter_microarray.md)
  : Painter 2018 Plasmodium falciparum 48h asexual lifecycle microarray
  data
- [`tradeSeq_BLASE_example_sce`](https://andrewmccluskey-uog.github.io/blase/reference/tradeSeq_BLASE_example_sce.md)
  : TradeSeq Example SCE for BLASE Vignette
- [`zhang_2021_heat_shock_bulk`](https://andrewmccluskey-uog.github.io/blase/reference/zhang_2021_heat_shock_bulk.md)
  : Zhang 2021 Plasmodium falciparum heat shock bulk data
