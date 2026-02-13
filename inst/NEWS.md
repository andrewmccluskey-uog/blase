# blase (development version)

* Initial submission.

## Major changes

* 

## Minor improvements and bug fixes

* 

# blase 1.2.0

## Major changes

* "confident_mapping" has now been renamed as "strong_mapping" to indicate that 
it is not a robust statistical confidence. The `confident_mapping()` getter 
function has been renamed to `strong_mapping()`
* `plot_mapping_result_heatmap()` `annotate_confidence` parameter has been renamed to `annotate_strong`.

## Minor improvements and bug fixes

* Now possible to use BLASE to create pseudotime bins from multiple trajectories at once.
* Speeds up by cell pseudotime bin assignment.
* Possible to apply different metrics for judging the best mapping. Spearman is still recommended.
* It is now possible to rename a bulk sample for a `MappingResult` object. 
* MappingResult heatmaps no longer convert the Y-axis titles (bulk sample name) to 
`character`, allowing numeric style scales.

# blase 1.0.0

* Initial release.
