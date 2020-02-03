# ideal 1.12.0

## New features

* Implemented a series of gentle fail mechanisms to avoid abrupt crashing of the app if the expected objects or formats are not provided

## Other notes

* Alignments of UI elements have been adjusted
* Updated citation info - the preprint for `ideal`, "ideal: an R/Bioconductor package for Interactive Differential Expression Analysis", is now published at bioRxiv, https://doi.org/10.1101/2020.01.10.901652

# ideal 1.10.0

## New features

* `ggplotCounts` gains a new parameter, `labels_repel`, to control the placement of the different labels in the plot - this can be useful when a large number of samples is available

## Bug fixes

* Fixed an error in the initialization of the app due to a new behavior introduced by `shinyAce` in version >= 0.4.0 - occurred in the same way as for `pcaExplorer`

## Other notes

* Better class checks via `is(...)` as per `BiocCheck` suggestion

# ideal 1.8.0

## New features

* `plot_ma` gains an additional parameter, `labels_repel`, for better placing the labels on the features to mark
* The pairwise correlation plots can now use logarithmic scale in the axes, use smaller subsets of the data for quicker inspection, and resizes the correlation info proportionally to its intensity
* The id types can now be chosen among the `keytypes` of the corresponding annotation packages (which still need to be installed when `ideal` is launched). Other input fields that specify id types also behave in a similar manner (e.g. in the Signature Explorer tab). This caused a problem for scenarios where common id types such as `ENSEMBL` are not available, like in Arabidopsis (where ids are provided often as `TAIR`) - thanks to Marc Galland for picking this up in https://github.com/federicomarini/ideal/issues/1
* The Signatures Explorer tab also has a fully fledged tour for first-time users, together with a collapsible help panel to describe its functionality in brief
* The zoomed MA-plot gains a new widget to control the labels size for the names of the genes
* It is possible to export the input data together with the results in a combined `SummarizedExperiment` object, which can be seamlessly fed into `iSEE` (http://bioconductor.org/packages/release/bioc/html/iSEE.html). 
  This leverages a new function, `wrapup_for_iSEE`, which is available and exported from `ideal`.

## Other notes

* An information box has been added to provide detailed information on the required input formats
* Added notification to specify how to install the airway package for demonstration purposes if not already available
* The Signatures Explorer tab is now displayed in a conditional panel, i.e. only when the required inputs are provided


# ideal 1.6.0

## New features

* Added extra diagnostic to results exploration (stratified pvalues histograms, schweder-spjotvoll plot)
* Added functionality for the tab Signatures Explorer: read in gmt files, and explore signatures for all/DE genes visually as heatmaps. The signatures can also be provided in the main call of the function as an argument, on top of uploading at runtime
* Different handling of the selection of ids/gene names, now not requiring anymore that the result is computed, but based on the dds object (and on the presence or not of the annotation object)

## Other notes

* Built project website via pkgdown, with customized reference structure
* Added a `NEWS.md` file to track changes to the package
* Updated the report template as well for including new functionality, and also updated vignette entry
* Correctly adding the resources to shinyBS, loaded via `.onLoad`
* Replaced the first tour structure and call
* New skin for the app, also with logo in the title header
* Instructions start collapsed for more compact main page
* Edited link in button to provide feedback, with subject specified

## Bug fixes

* Fixed behavior of rendering inline the content of the report - did not work properly for server deployed instances

# ideal 1.4.0

## New features

* Specified single go term selection for generating heatmaps of gene signatures
* Added support for logFC shrinkage, following the latest devels of DESeq2

## Bug fixes

* Corrected output for the vignette, as html_document2 is now deprecated
* Menus are back in the expanded form
* Fixed the behavior with addMLE

## Other notes

* Added further progress indicators to give feedback during lengthy steps
* Improved ggplotCounts for better scale display, using exact arg matching, defaulting to the transformed counts

# ideal 0.99.0

## New features

* Ready for Bioc submission
* Completed the news

# ideal 0.9.1

## New features

* Added Instructions fully from rendered version of the vignette to have available at runtime
* Added support for downloading all plots and tables

# ideal 0.9.0

## New features

* Interactive tours are covering now all tabs, with extensive walkthroughs for the user
* Added all screenshots to vignette

# ideal 0.6.2

## New features

* Interactive tours are now available, coded in external files
* Travis-CI is now supported

# ideal 0.6.0

## New features

* Added MA plot with extra custom list to avoid manual selection of many genes
* MA plot function now automatically supports subset of gene to be extra plotted
* Added documentation with roxygen to all functions
* Heatmap functions for genes annotated to a GO term as signature
* Template report also provided
* Full draft of vignette now available, working towards bioc submission
* Added textual help to all sections, with collapsible element
* Added proof of principle to have interactive tours based on rintrojs

# ideal 0.4.0

## New features

* Gene box info added, based on rentrez
* New look for MA plots and volcano plots

# ideal 0.3.0

## New features

* Restructuring of the folders done, package can be correctly installed, loaded - namespace, description are set up
    
# ideal 0.2.0

## New features

* Correct structure of the package

# ideal 0.1.0

## New features

* Package created!
