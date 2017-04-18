## Functional Analysis

Do you need something more than *just* a list of genes? Looking for insight on the affected biological pathways? Then you can use the modules in the Functional Analysis tab. You can do Gene Ontology overrepresentation analysis based on

- `limma::goana`
- `topGO` (recommended for the somewhat clearer presentation of results, thanks to the algorithms in the topGO package)
- `goseq`

You can analyse

- *all* regulated genes, up- and down- regulated
- up-regulated genes alone
- down-regulated genes only
- two custom lists, that can be uploaded as text files

Once you obtain the interactive table of functions enriched in the gene set of interest, you can click on any row of the DataTable, and this will display a heatmap with the expression values for the genes annotated to that GO Term - as a kind of signature for that function in your data.

You can explore the overlap of the lists as Venn diagrams as well as Upset plots.

Once you are done with this help box, you can close it by clicking on the `Help` text.
