import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Chromatin Conformation, HiC", divider=True)

    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Look at the files for today's class from /project/60006/hgen_share/lec12")
            
            st.code("data=/project/60006/hgen_share/lec12", language="bash")
        
    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Generate an uniformly binned contact matrix from NPC_1.chr19.pairs.gz")
            
            st.markdown("Check out the first few lines of the file")
            st.code("""
            zcat ${data}/NPC_1.chr19.pairs.gz | head | column -t
            """, language="bash")

            st.markdown("Create matrix from the pairix-indexed and sorted file with cooler cload pairs. You can cancel with `Ctrl + C` because this will take a while")
            st.code("""
            cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 ${data}/mm10.chrom.sizes:10000 ${data}/NPC_1.chr19.pairs.gz NPC_1.chr19.cool
            """, language="bash")

            st.markdown("Examine the output file's metadata using `cooler info`")
            st.code("""
            cooler info NPC_1.chr19.cool
            """, language="bash")

            st.markdown("Perform matrix balancing through `cooler balance`")
            st.code("""
            cooler balance NPC_1.chr19.cool -p 1
            """, language="bash")
            
            st.markdown("Look at the underlying data through `cooler dump`")
            st.code("""
            cooler dump --join --balanced --header NPC_1.chr19.cool | sed -n '1p;12500,12520p;12521q' | column  -t
            """, language="bash")
            
            st.markdown("Coarsen the matrix with `cooler zoomify`")
            st.code("""
            cooler zoomify -p 1 -r 10000,100000 --balance -o NPC_1.chr19.mcool NPC_1.chr19.cool
            """, language="bash")
            
            st.markdown("List the contents inside the mcool file with cooler ls")
            st.code("""
            cooler ls NPC_1.chr19.mcool
            """, language="bash")
            
            st.markdown("Visualize the results using `cooler show`")
            st.code("""
            cooler show NPC_1.chr19.mcool::/resolutions/10000 -b -o 10kb.ice.png chr19
            cooler show NPC_1.chr19.mcool::/resolutions/100000 -b -o 100kb.ice.png chr19
            cooler show NPC_1.chr19.mcool::/resolutions/10000 -o 10kb.raw.png chr19
            cooler show NPC_1.chr19.mcool::/resolutions/100000 -o 100kb.raw.png chr19
            """, language="bash")

            st.markdown("Download results to your local computer, specifically files ending with .png")
            
    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Calculate a contact probability decay curve for the generated mcool file")
            
            st.markdown("Take the average across sub-diagonals with `cooltools compute-expected`")
            st.code("""
            cooltools expected-cis NPC_1.chr19.mcool::/resolutions/10000 -o NPC_1.chr19.cis.tsv -p 1
            """, language="bash")

            st.markdown("Peek at the output")
            st.code("""
            head NPC_1.chr19.cis.tsv | column -t
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Perform eigendecomposition on the same matrix")
            
            st.markdown("Compute the first few eigenvectors using `cooltools call-compartments` with GC content as the reference for sign-flipping")
            st.code("""
            cooltools eigs-cis --n-eigs 3 -o NPC_1.chr19 --phasing-track ${data}/100000.gc NPC_1.chr19.mcool::/resolutions/100000
            """, language="bash")

            # st.markdown("Now compute the first few eigenvectors for trans matrices. See [cooltools](https://cooltools.readthedocs.io/en/latest/cli.html#cooltools-expected-trans)")
            # st.markdown("Look at the outputs")
            # st.code("""
            # cat NPC_1.chr19.cis.lam.txt
            # head NPC_1.chr19.cis.vecs.tsv
            # """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Extra: Visualize eigenvectors using a saddle-plot")
            
            # st.markdown("Use `python` version of cooltools")
            # st.code("""
            # # import standard python libraries
            # import numpy as np
            # import matplotlib.pyplot as plt
            # import pandas as pd
            # import os, subprocess
            # import cooler
            # import cooltools.lib.plotting
            # import cooltools

            # # read in cooler
            # clr = cooler.Cooler('/lustre06/project/6007495/padilr1/projects/HGEN663_W21/pub/lec12_extra/NPC.mcool::/resolutions/100000')
            # # generate dataframe of gc coverage in your organism, here it is mm10 or mouse
            # import bioframe
            # bins = clr.bins()[:]
            # mm10_genome = bioframe.load_fasta('/lustre06/project/6007495/padilr1/pipelines/hic/dcHiC/demo/dcHiC_demo/ESC_NPC_CN_100Kb/mm10_100000_goldenpathData/mm10.fa');
            # gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], mm10_genome)
            # # gc_cov.to_csv('mm10_gc_cov_100kb.tsv',index=False,sep='\t') = this script can export the gc coverage dataframe into a csv file

            # ## Cooltools also allows a view to be passed for eigendecomposition to limit to a certain set of regions. The following code creates the simplest view, of the two chromosomes in this cooler
            # view_df = pd.DataFrame({'chrom': clr.chromnames,
            #                         'start': 0,
            #                         'end': clr.chromsizes.values,
            #                         'name': clr.chromnames}
            #                     )
            # ## obtain first 3 eigenvectors
            # cis_eigs = cooltools.eigs_cis(
            #                         clr,
            #                         gc_cov,
            #                         view_df=view_df,
            #                         n_eigs=3,
            #                         )

            # # cis_eigs[0] returns eigenvalues, here we focus on eigenvectors
            # eigenvector_track = cis_eigs[1][['chrom','start','end','E1']]
            # # calculate the expected cis values
            # cvd = cooltools.expected_cis(
            #         clr=clr,
            #         view_df=view_df,
            # )
            # #
            # Q_LO = 0.025 # ignore 2.5% of genomic bins with the lowest E1 values
            # Q_HI = 0.975 # ignore 2.5% of genomic bins with the highest E1 values
            # N_GROUPS = 38 # divide remaining 95% of the genome into 38 equisized groups, 2.5% each

            # # saddle then plots two matrices: one with the sum for each pair of categories, interaction_sum, and the other with the number of bins for each pair of categories, interaction_count. Typically, interaction_sum/interaction_count is visualized
            # interaction_sum, interaction_count =  cooltools.saddle(
            #         clr,
            #         cvd,
            #         eigenvector_track,
            #         'cis',
            #         n_bins=N_GROUPS,
            #         qrange=(Q_LO,Q_HI),
            #         view_df=view_df
            # )

            # # parameters and functions required to plot the saddle plot
            # import warnings
            # from cytoolz import merge

            # def saddleplot(
            #     track,
            #     saddledata,
            #     n_bins,
            #     vrange=None,
            #     qrange=(0.0, 1.0),
            #     cmap="coolwarm",
            #     scale="log",
            #     vmin=0.5,
            #     vmax=2,
            #     color=None,
            #     title=None,
            #     xlabel=None,
            #     ylabel=None,
            #     clabel=None,
            #     fig=None,
            #     fig_kws=None,
            #     heatmap_kws=None,
            #     margin_kws=None,
            #     cbar_kws=None,
            #     subplot_spec=None,
            # ):
            #     from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
            #     from matplotlib.colors import Normalize, LogNorm
            #     from matplotlib import ticker
            #     import matplotlib.pyplot as plt
            #     class MinOneMaxFormatter(ticker.LogFormatter):
            #         def set_locs(self, locs=None):
            #             self._sublabels = set([vmin % 10 * 10, vmax % 10, 1])
                        
            #         def __call__(self, x, pos=None):
            #             if x not in [vmin, 1, vmax]:
            #                 return ""
            #             else:
            #                 return "{x:g}".format(x=x)
                
            #     track_value_col = track.columns[3]
            #     track_values = track[track_value_col].values
                
            #     digitized_track, binedges = cooltools.digitize(
            #         track, n_bins, vrange=vrange, qrange=qrange
            #     )
            #     x = digitized_track[digitized_track.columns[3]].values.astype(int).copy()
            #     x = x[(x > -1) & (x < len(binedges) + 1)]
                
            #     groupmean = track[track.columns[3]].groupby(digitized_track[digitized_track.columns[3]]).mean()
                
            #     if qrange is not None:
            #         lo, hi = qrange
            #         binedges = np.linspace(lo, hi, n_bins + 1)
                
            #     n = saddledata.shape[0]
            #     X, Y = np.meshgrid(binedges, binedges)
            #     C = saddledata
            #     if (n - n_bins) == 2:
            #         C = C[1:-1, 1:-1]
            #         groupmean = groupmean[1:-1]
                
            #     if subplot_spec is not None:
            #         GridSpec = partial(GridSpecFromSubplotSpec, subplot_spec=subplot_spec)
            #     grid = {}
            #     gs = GridSpec(
            #         nrows=3,
            #         ncols=3,
            #         width_ratios=[0.2, 1, 0.1],
            #         height_ratios=[0.2, 1, 0.1],
            #         wspace=0.05,
            #         hspace=0.05,
            #     )
                
            #     if fig is None:
            #         fig_kws_default = dict(figsize=(5, 5))
            #         fig_kws = merge(fig_kws_default, fig_kws if fig_kws is not None else {})
            #         fig = plt.figure(**fig_kws)
                    
            #     if scale == "log":
            #         norm = LogNorm(vmin=vmin, vmax=vmax)
            #     elif scale == "linear":
            #         norm = Normalize(vmin=vmin, vmax=vmax)
            #     else:
            #         raise ValueError("Only linear and log color scaling is supported")
            #     grid["ax_heatmap"] = ax = plt.subplot(gs[4])
            #     heatmap_kws_default = dict(cmap="coolwarm", rasterized=True)
            #     heatmap_kws = merge(
            #         heatmap_kws_default, heatmap_kws if heatmap_kws is not None else {}
            #     )
            #     img = ax.pcolormesh(X, Y, C, norm=norm, **heatmap_kws)
            #     plt.gca().yaxis.set_visible(False)
                
            #     margin_kws_default = dict(edgecolor="k", facecolor=color, linewidth=1)
            #     margin_kws = merge(margin_kws_default, margin_kws if margin_kws is not None else {})
            #     grid["ax_margin_y"] = plt.subplot(gs[3], sharey=grid["ax_heatmap"])
                
            #     plt.barh(
            #         binedges, height=1/len(binedges), width=groupmean, align="edge", **margin_kws
            #     )
                
            #     plt.xlim(plt.xlim()[1], plt.xlim()[0])  # fliplr
            #     plt.ylim(hi, lo)
            #     plt.gca().spines["top"].set_visible(False)
            #     plt.gca().spines["bottom"].set_visible(False)
            #     plt.gca().spines["left"].set_visible(False)
            #     plt.gca().xaxis.set_visible(False)
            #     grid["ax_margin_x"] = plt.subplot(gs[1], sharex=grid["ax_heatmap"])
                
            #     plt.bar(
            #         binedges, width=1/len(binedges), height=groupmean, align="edge", **margin_kws
            #     )
                
            #     plt.xlim(lo, hi)
            #     plt.gca().spines["top"].set_visible(False)
            #     plt.gca().spines["right"].set_visible(False)
            #     plt.gca().spines["left"].set_visible(False)
            #     plt.gca().xaxis.set_visible(False)
            #     plt.gca().yaxis.set_visible(False)
                
            #     grid["ax_cbar"] = plt.subplot(gs[5])
            #     cbar_kws_default = dict(fraction=0.8, label=clabel or "")
            #     cbar_kws = merge(cbar_kws_default, cbar_kws if cbar_kws is not None else {})
            #     if scale == "linear" and vmin is not None and vmax is not None:
            #         grid["ax_cbar"] = cb = plt.colorbar(img, **cbar_kws)
            #         decimal = 10
            #         nsegments = 5
            #         cd_ticks = np.trunc(np.linspace(vmin, vmax, nsegments) * decimal) / decimal
            #         cb.set_ticks(cd_ticks)
            #     else:
            #         print('cbar')
                    
            #         cb = plt.colorbar(img, format=MinOneMaxFormatter(), cax=grid["ax_cbar"], **cbar_kws)
            #         cb.ax.yaxis.set_minor_formatter(MinOneMaxFormatter())
            #     grid["ax_heatmap"].set_xlim(lo, hi)
            #     grid["ax_heatmap"].set_ylim(hi, lo)
            #     grid['ax_heatmap'].grid(False)
            #     if title is not None:
            #         grid["ax_margin_x"].set_title(title)
            #     if xlabel is not None:
            #         grid["ax_heatmap"].set_xlabel(xlabel)
            #     if ylabel is not None:
            #         grid["ax_margin_y"].set_ylabel(ylabel)
                
            #     return grid

            # ## finally, plot the saddle-plot
            # img = saddleplot(eigenvector_track,
            #         interaction_sum/interaction_count,
            #         N_GROUPS,
            #         qrange=(Q_LO,Q_HI),
            #         cbar_kws={'label':'average observed/expected contact frequency'}
            #         );

            # ## save figure
            # plt.savefig('/lustre06/project/6007495/padilr1/projects/HGEN663_W21/pub/lec12_extra/final_analysis/saddleplot/NPC.saddle.cis.100000.png')
            # """, language="python")

            st.markdown("In ${data}/create.saddle.plot.py we have the code to create a saddleplot. Inspect this file")
            st.code("less ${data}/create.saddle.plot.py", language="bash")

            st.code("python ${data}/create.saddle.plot.py", language="bash")

            st.image("images/lec12.NPC.saddle.png")

    st.divider()
    #############################################################################################

    # with st.container(border=True):
    #         st.markdown("#### Extra: Run differential compartment analysis")
            
    #         st.markdown("Differential compartment analysis using [dcHiC](https://github.com/ay-lab/dcHiC)")
    #         st.code("""
    #         ###                  STEP 1                  ###
    #         ### Step 1: dcHiC will create the raw PC files for each chromosome. ###

    #         Rscript ../scripts/dchicf.r --file input.ES_NPC.txt --pcatype cis --dirovwt T --cthread 2 --pthread 4

    #         ###                  STEP 2                  ###
    #         ### Step 2: dcHiC will select the best pc out of PC1 and PC2 for each chromosome 
    #         ### by comparing each one against GC content and gene density through correlation.

    #         Rscript ../scripts/dchicf.r --file input.ES_NPC_CN.txt --pcatype select --dirovwt T --genome mm10 --cthread 8 --pthread 4 --gfolder /lustre06/project/6007495/padilr1/pipelines/hic/dcHiC/demo/dcHiC_demo/ESC_NPC_CN_100Kb/mm10_100000_goldenpathData

    #         ###                  STEP 3                  ###
    #         ### Step 3: dcHiC will use the selected PC and qunatile normalize them.
    #         ### This will be followed by mahalanobis distance calculation and outlier detection.
    #         ### If replicates are available, dcHiC will apply IHW (independent hypothesis weighting) to adjust the significance. 
    #         ### This will create a directory "DifferentialResult" and few more subfolder. 
    #         ### differential.intra_sample_group.pcOri.bedGraph shows the significance score for all the Hi-C bins.
    #         ### differential.intra_sample_group.Filtered.pcOri.bedGraph shows the filtered (padj < 0.1) the Hi-C bins.
    #         ### The differential.intra_sample_combined.*.bedGraph file shows the compartment scores and significance 
    #         ### of each replicates, differential.intra_sample_group.*.bedGraph files are the subset of these.
    #         ### *.pcOri.* are files with original pc values i.e. before quantile normalized values.
    #         ### *.pcQnm.* are files with quantile normalized values. 
    #         ### Note: All the differential compartment analysis is performed using *.pcQnm.* files. 
    #         ### *.pcOri.* files are generated for vizualization purpose.

    #         Rscript ../scripts/dchicf.r --file input.ES_NPC_CN.txt --pcatype analyze --dirovwt T --diffdir ES_NPC_100Kb


    #         ###                  STEP 4                 ###
    #         ### Step 6: dcHiC will generate the standalone IGV web page.
    #         ### The step reqires "create_datauri" code provided within scripts directory.
    #         ### Please create a soft-link to your ln -s ./scripts/create_datauri ~/.local/bin/

    #         Rscript ../scripts/dchicf.r --file input.ES_NPC.txt --pcatype viz --dirovwt T --diffdir ES_NPC_100Kb --genome mm10 --pcgroup pcOri


    #         ###                  STEP 5              ###
    #         ### Step 7: dcHiC will perform enrichment with the genes overlapping with the 
    #         ### differential A-compartments in each sample. This step will generate a 
    #         ### geneEnrichment directory under DifferentialResult/
    #         ### Under geneEnrichment directory there will sub-directories named as *_geneEnrichment.
    #         ### We recomment using the *_geneList.anchor.txt (Entrez IDs) file to perform the gene-enrichment using
    #         ### https://toppgene.cchmc.org/enrichment.jsp .


    #         Rscript ../scripts/dchicf.r --file input.ES_NPC.txt --pcatype enrich --genome mm10 --diffdir ES_NPC_100Kb --exclA F --pcscore T --region anchor --pcgroup pcOri --gfolder mm10_100000_goldenpathData
    #         """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Identify contact-insulating loci also for the same matrix")
            
            st.markdown("Compute insulation score via cooltools insulation")
            st.code("""
            cooltools insulation NPC_1.chr19.mcool::/resolutions/10000 100000 > NPC_1.chr19.ins.tsv
            """, language="bash")

            st.markdown("Inspect the output")
            st.code("""
            head NPC_1.chr19.ins.tsv
            """, language="bash")

    st.divider()
    #############################################################################################

    # with st.container(border=True):

    #         st.markdown("#### Find dots on again the same matrix")

    #         st.markdown("Use the expected values computed before together with cooltools call-dots")
    #         st.code("""
    #         cooltools dots NPC_1.chr19.mcool::/resolutions/10000 NPC_1.chr19.cis.tsv -o NPC_1.chr19.dots -p 1
    #         """, language="bash")

    #         st.markdown("Look at the outputs")
    #         st.code("""
    #         head NPC_1.chr19.dots
    #         """, language="bash")

    #         st.markdown("Pile up the loops using coolpup.py")
    #         st.code("""
    #         cut -f 1-6 NPC_1.chr19.dots > NPC_1.chr19.dots.bedpe

    #         coolpup.py NPC_1.chr19.mcool::/resolutions/10000 NPC_1.chr19.dots.bedpe \\
    #                     --features_format bedpe \\
    #                     --outname NPC.cis_10000.pileup.clpy --nproc 1 \\
    #                     --expected NPC_1.chr19.cis.tsv
    #         """, language="bash")
            
    #         st.markdown("Plot the loop pile-up using `plotpup.py`")
    #         st.code("""
    #         plotpup.py --input_pups NPC.cis_10000.pileup.clpy --output NPC.cis_10000.pileup.png
    #         """, language="bash")
            
    #         st.markdown("**Extra**: Perform differential loop analysis on union of loops using [pareidolia](https://github.com/koszullab/pareidolia)")
    #         st.code("""
    #         pareidolia -n 8 -b /lb/project/GRID/padilr1/projects/HGEN663/readpileup/union.dots.bedpe -k loops \\
    #                     /lb/project/GRID/padilr1/projects/HGEN663/lec12_extra/ESC_1.mcool::resolutions/10000,/lb/project/GRID/padilr1/projects/HGEN663/lec12_extra/ESC_2.mcool::resolutions/10000,/lb/project/GRID/padilr1/projects/HGEN663/lec12_extra/ESC_3.mcool::resolutions/10000,/lb/project/GRID/padilr1/projects/HGEN663/lec12_extra/ESC_4.mcool::resolutions/10000,/lb/project/GRID/padilr1/projects/HGEN663/lec12_extra/NPC_1.mcool::resolutions/10000,/lb/project/GRID/padilr1/projects/HGEN663/lec12_extra/NPC_2.mcool::resolutions/10000,/lb/project/GRID/padilr1/projects/HGEN663/lec12_extra/NPC_3.mcool::resolutions/10000,/lb/project/GRID/padilr1/projects/HGEN663/lec12_extra/NPC_4.mcool::resolutions/10000 \\
    #                     ESC,ESC,ESC,ESC,NPC,NPC,NPC,NPC \
    #                     output.pareidolia.tsv
    #         """, language="bash")
            
    #         


    with st.container(border=True):
        st.markdown("Lastly, we can view the results of our HiC analysis in [HiGlass](http://206.12.101.70:8888/app)")

    st.divider()

    st.markdown("#### Let's jump to R!!!")

#########################################################################33
    with st.container(border=True):

        st.markdown("#### Set up")

        st.code("""
        library(data.table)
        library(tidyverse)
        library(smacof)
        library(DESeq2)
        library(InteractionSet)
        library(diffHic)
        library(broom)
        library(highcharter)
        library(heatmaply)
        library(plotly)
        library(eulerr)
        library(reactable)
        library(knitr)
        # ("data.table", "tidyverse", "smacof", "DESeq2", "InteractionSet", "diffHic",
        #       "broom", "highcharter", "heatmaply", "plotly", "eulerr")
        # set your working directory
        setwd("~/Documents/HGEN_663/extra/lec12")
        """, language="r")

    st.divider()
#########################################################################33
    with st.container(border=True):

        st.markdown("#### Import")

        st.code("""
        load('lec12.RData')
        """, language="r")

    st.divider()
#########################################################################33
    with st.container(border=True):

        st.markdown("#### Contact matrices")

        st.markdown("There’s a number of ways through which we can examine the correlation matrix")
        st.markdown("##### Hierarchical clustering")

        st.code("""
        tibble(s = names(hr)) %>%
        mutate(type = factor(sub('_.*', '', s))) %>%
        column_to_rownames('s') %>%
        heatmaply(hr, row_side_colors = ., col_side_colors = .,
                    label_names = c("Sample 1", "Sample 2", "SCC"))
        """, language="r")
        st.image("images/lec12.hier.clust.png")

        st.markdown("#### Multidimensional scaling")
        st.code("""
        sim2diss(hr, method = 'corr') %>%
                        mds(ndim = 2) %>%
                        .$conf %>%
                        as.data.frame() %>%
                        rownames_to_column('s') %>%
                        mutate(type = sub('_.*', '', s)) %>%
                        plot_ly(x = ~D1, y = ~D2, color = ~type, text = ~s) %>%
                        add_markers(legendgroup = ~type) %>%
                        add_text(textposition = 'top left', showlegend = F, legendgroup = ~type)
        """, language="r")
        st.image("images/lec12.multiscaling.png")

    st.divider()
#########################################################################33
    with st.container(border=True):

        st.markdown("#### P(s)")

        st.markdown("For the contact probability decay curves and their derivatives that we’ve computed before, we can again just simply visualize them")
        st.code("""
        # colors for each cell type
        hclrs <- c(ESC = "#7cb5ec", NPC = "#434348")

        # plot data
        pd <- ps %>%
        lapply(function(x) {
            d <- rbindlist(x, idcol = 'samp')
            ifelse('slope' %in% names(d), 'slope', 'balanced.avg') %>%
            c('samp', 's_bp', .) %>%
            d[, ., with = F] %>%
            `colnames<-`(c('samp', 'x', 'y'))
        }) %>%
        rbindlist(idcol = 'kind') %>%
        group_by(kind, samp) %>%
        do(data = list_parse2(data.frame(.$x, .$y))) %>%
        ungroup() %>%
        separate(samp, c('ctype', 'rep'), '_', F, fill = 'right') %>%
        #dplyr::filter(is.na(rep)) %>%
        mutate(color = hclrs[ctype],
                ctype = factor(ctype, names(hclrs))) %>%
        arrange(ctype) %>%
        mutate(name = samp,
                id = ifelse(kind != 'log', tolower(name), NA),
                linkedTo = ifelse(kind == 'log', tolower(name), NA),
                yAxis = ifelse(kind == 'log', 0, 2))

        highchart() %>%
        # we use log scale for P(s) and linear for the slope (which was already taken in log space)
        hc_yAxis_multiples(list(title = list(text = "<b>Contact probability</b>, P<sub>c</sub>(s)",
                                            useHTML = TRUE),
                                type = "logarithmic",
                                labels = list(formatter = JS("function(){return this.value.toExponential(0);}")),
                                height = '45%', top = '0%', offset = 0),
                            list(height = '10%', top = '45%',
                                title = NULL,
                                plotLines = list(
                                    list(color = "#a9a9a9", width = 2,
                                        value = .5, zIndex = 1)
                                ),
                                labels = list(enabled = F),
                                gridLineWidth = 0,
                                min = 0, max = 1),
                            list(type = "linear",
                                title = list(
                                    text = "<b>Slope</b>, <sup>d</sup>&frasl;<sub>ds</sub> log P<sub>c</sub>(s)",
                                    useHTML = TRUE
                                ),
                                height = '45%', top = '55%', offset = 0)) %>%
        # grab the data
        hc_add_series_list(pd) %>%
        # a bit of formatting
        hc_xAxis(type = "logarithmic",
                title = list(text = "<b>Genomic separation</b> (bp), s",
                                useHTML = T),
                minorTickInterval = 'auto',
                min = 1e4, max = 1e8) %>%
        hc_tooltip(headerFormat = '<span style="font-size: 10px">{point.key:,.0f} bp</span><br/>') %>%
        hc_chart(zoomType = "xy") %>%
        hc_plotOptions(line = list(marker = list(enabled = F)))
        """, language="r")
        st.image("images/lec12.Ps.png")

    st.divider()
#########################################################################33
    with st.container(border=True):

        st.markdown("#### Compartments")

        st.markdown("Next we’ll visualize the first eigenvector (i.e., “compartment score”)")
        st.code("""
        imageList <- list.files("~/Documents/HGEN_663/extra/lec12", pattern= "saddleplot.png", full.names=T)
        include_graphics(imageList,dpi = 600)
        """, language="r")
        st.image("images/lec12.saddle.png")

    st.divider()
#########################################################################33
    with st.container(border=True):

        st.markdown("#### Insulation")

        st.markdown("For the simplest comparison we can just count the number of shared boundaries")
        st.code("""
        bdrs <- ins[c('ESC', 'NPC')] %>%
        lapply(function(x) {
            na.omit(x) %>%
            dplyr::filter(boundary_strength_100000 > .5) %>%
            mutate(start = start + 1) %>%
            makeGRangesFromDataFrame()
        })

        Reduce(c, bdrs) %>% 
        unique() %>%
        {lapply(bdrs, function(x) overlapsAny(., x))} %>%
        bind_cols() %>%
        euler() %>%
        plot(quantities = T)
        """, language="r")
        st.image("images/lec12.insulation.png")

    st.divider()
#########################################################################33
    with st.container(border=True):

        st.markdown("#### Loops")
        st.markdown("#### Called Loops")

        st.markdown("These are loops called separately in each sample. There are specific classes in R that handles paired ranges, one of which is GInteractions")
        st.code("""
        lps <- dots[c('ESC', 'NPC')] %>%
        lapply(function(x) {
            list(x[,1:3], x[,4:6]) %>%
            lapply(function(y) {
                y %>%
                `colnames<-`(c('chr', 'start', 'end')) %>%
                mutate(start = start + 1) %>%
                makeGRangesFromDataFrame()
            }) %>%
            {GInteractions(.[[1]], .[[2]], mode = 'reverse')}
        })
        """, language="r")

        st.markdown("ESC")
        st.code("""
        lps$ESC 
        """, language="r")
        st.image("images/lec12.called.loops.png")

        st.markdown("NPC")
        st.code("""
        lps$NPC
        """, language="r")
        st.image("images/lec12.called.loops2.png")
        
    st.divider()
#########################################################################33
    with st.container(border=True):

        st.markdown("#### Loop pile-up")

        st.code("""
        imageList <- list.files("~/Documents/HGEN_663/extra/lec12", pattern= "pileup.png", full.names=T)
        include_graphics(imageList,dpi = 600)
        """, language="r")
        st.image("images/lec12.loop.pile.up.png")

    st.divider()
#########################################################################33
    with st.container(border=True):

        st.markdown("#### Differential loop calling")

        st.code("""
        diff_loop <- fread("output.pareidolia.tsv")
        diff_loop %>% reactable()
        """, language="r")
        st.image("images/lec12.diff.loops.png")

    st.divider()