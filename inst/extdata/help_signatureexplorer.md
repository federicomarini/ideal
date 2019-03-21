## Signatures Explorer

This tab allows you to check the behavior of a number of provided gene signatures in your data at hand, displaying this as a heatmap.

This panel is composed by different well panels:

- in the Setup Options, you can select and upload a gene signature file, in `gmt` format (e.g. like the ones provided in the MSigDB database, or from WikiPathways), and quickly compute the variance stabilized transformed version of your data, which is more amenable for visualization than raw or normalized counts

- in the Conversion options tab, you can create an annotation vector, used to bring the ids from your data and the ids the `gmt` used for encoding the signature elements. 
  This works based on the `org.XX.eg.db` packages.

- the lower well panels control the appearance of the heatmap, also with an option to display all genes annotated in that pathway, or only the ones detected as differentially expressed (for this you need to provide or compute the result object)