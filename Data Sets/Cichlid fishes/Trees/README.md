# Ancient and recent hybridization in the <i>Oreochromis</i> cichlid fishes

[https://doi.org/10.5061/dryad.p2ngf1w0c](https://doi.org/10.5061/dryad.p2ngf1w0c)

## Description of the data and file structure

All trees inferred, along with the dataset used to generate them, are available in the trees folder, for both M. zebra and O. niloticus mapping dataset.

Mzeb_nc_iqtree.min4.phy.treefile, Onil_nc_iqtree.min4.phy.treefile - The iqtree inferred phylogenetic trees for the nuclear genome (non-coding regions) inferred for the M. zebra (Mzeb) and O. niloticus (Onil) mapping datasets. respectively. Intermediate files generated are also available for the iqtree runs (Mzeb_nc_iqtree.min4.phy.bionj, Mzeb_nc_iqtree.min4.phy.ckp.gz, Mzeb_nc_iqtree.min4.phy.contree, Mzeb_nc_iqtree.min4.phy.log, Mzeb_nc_iqtree.min4.phy.mldist, Mzeb_nc_iqtree.min4.phy.model.gz, Mzeb_nc_iqtree.min4.phy.runtrees, Mzeb_nc_iqtree.min4.phy.splits.nex, Onil_nc_iqtree.min4.phy.bionj, Onil_nc_iqtree.min4.phy.ckp.gz, Onil_nc_iqtree.min4.phy.contree, Onil_nc_iqtree.min4.phy.log, Onil_nc_iqtree.min4.phy.mldist, Onil_nc_iqtree.min4.phy.model.gz, Onil_nc_iqtree.min4.phy.runtrees, Onil_nc_iqtree.min4.phy.splits.nex

Mzeb_bs10_0.5.trees, Onil_bs10_0.5.trees - files containing lists of trees inferred on 10kb recombination-free units as input for ASTRAL phylogenetic analysis.

Mzeb_raw.trees, Onil_raw.trees - the treefiles used prior to collapse of poorly supported nodes and treeshrink outlier removal.

Astral multi-species coalescent trees (in ASTRAL). The treefiles used for these trees are in Mzeb_bs10_0.5.trees and Onil_bs10_0.5.trees. Output trees with all individuals are in Onil_bs10_0.5_astral.out and Mzeb_bs10_0.5_astral.out. Tips are collapsed to the population/ species level in Onil_bs10_0.5_astral_mono.out and Mzeb_bs10_0.5_astral_mono.out. These trees with added information on each quartet frequency are in Onil_bs10_0.5_astral_mono_t8.out and Mzeb_bs10_0.5_astral_mono_t8.out

A neighbour-joining tree from all non-morphological hybrids is in NJ_all_nonmorphhybrid.Nwk. The raw and converted distance matrix are in njtree_dist_raw.mat and njtree_dist_conv.mat, with the R script used to generate it in njtree_dist_conv.mat

A mitochondrial tree for all non-morphological hybrids for which a mitochondrial genome could be assembled is in NoHybridMitoset_longassembly_mars_clean0.5.aln.treefile, with the fasta alignment used to generate this in NoHybridMitoset_longassembly_mars_clean0.5.aln. The iqtree intermediate files are also given (NoHybridMitoset_longassembly_mars_clean0.5.aln.bionj,NoHybridMitoset_longassembly_mars_clean0.5.aln.ckp.gz,NoHybridMitoset_longassembly_mars_clean0.5.aln.contree,NoHybridMitoset_longassembly_mars_clean0.5.aln.iqtree,NoHybridMitoset_longassembly_mars_clean0.5.aln.log,NoHybridMitoset_longassembly_mars_clean0.5.aln.mldist,NoHybridMitoset_longassembly_mars_clean0.5.aln.model.gz,NoHybridMitoset_longassembly_mars_clean0.5.aln.runtrees,NoHybridMitoset_longassembly_mars_clean0.5.aln.splits.nex,NoHybridMitoset_longassembly_mars_clean0.5.aln.uniqueseq.phy).

## Code/Software

Custom python scripts and R scripts are provided in Zenodo. What each python script does is described in Notes.txt and below.

parse_cutadapt.py, parse_cutadapt_sra.py* *runs cutadapt on reads when necessary (i.e. adapters were found)

average_qual_perscaf.py gives average QUAL per scaffold in vcf file

CollapseBSNodes.py collapses nodes with SH-like support below threshold (from iq-tree output trees)

GetTreesForPhylo.py Removes Trees which are adjacent to each other following filtering for recombination (to avoid pseudo-replication)

AlignmentFilter.py filters alignments for missing data per column, removes gappy sequences and short alignments

tree_dists.py gives Robinson-Foulds distances between pairs of trees

rename_chroms_admixture.py - renames .bim file to numeric so can be used for admixture

SpecSpecificSNPs.py - species specific SNPs from pop file and vcf

SpecSNPsperSample.py - which reference species SNPs are present in a test sample

filter_beagle.py - filters beagle output for unlinked sites as part of batch test analysis

Get10kbwindows.py - gives 10kb window bed file to generate trees from

SplitFastasByChromosome.py - splits ANGSD output for multiple individuals into per-scaffold alignments

sig_in_both.py - gives trios significant in both Mzeb and Onil mapping D statistics (Dsuite output)

dtree_addDp.py - add Dp values to Dtree output

get_afs.py - get allele frequency for SNP within specified population (from pop file) - used for getting uro/leu/nil specific SNPs

blt_dct_get_sig.py. blt_dct_sig_pairs.py, blt_dct_sig_pairs_both.py - get significant sets of results from blt and dct mapping, and those which are significant in both Mzeb and Onil mapping.

Dtree_get_sig_results.py - get significant result count from Dtree analysis

get_fb_pval_nodup.py, get_fb_pval_dists_nodup.py  - run the permutation tests to identify if closely/ distantly related species or species sharing drainage basins more likely to have highest fbranch than expected

get_tree_dists.py - get distances between populations on a phylogenetic tree

get_max_fbranch.py - find species pairs with highest fbranch values, and associated information

get_introtopobed.py - find regions where Dwt = 1 over 50kb+ - output bed

plot_Dwt_twisst.R - R script giving the commands used to generate twisst plots, to demonstrate how to use the Dwt test.

plot_twisst_functions.R *- *R script*  *edit of plot_twisst.R from https://github.com/simonhmartin/twisst/blob/master/plot_twisst.R has added function used for this paper - namely plot.twisst.topodiff to give the figure in the main text showing areas with excess of the introgression tree.

plot_map.R - R script used to generate map used in figure 1.
