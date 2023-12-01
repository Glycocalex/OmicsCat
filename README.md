# OmicsCat
**Concatenation integration for generic multi-omics datasets.**


![image5](https://github.com/Alex-RW-Bennett/OmicsCat/assets/131264603/444e5063-b630-4efa-b80d-84eae3d0d353)


Hi, thanks for having a look at OmicsCat!

OC is very much a WIP, however I am working to develop this into a simple package for temporal analysis of biological data. Though robust tools for the concatenation of generic -omics data are available, as far as I am aware no such tools have been developed with a focus on biological rhythms. As well as good python-practice for myself, I am hoping this will eventually be a useful, quick-and-easy tool for basic multi-omics integration for data collected over the course of circadian or other biological rhythms.

As for the definition of 'Generic-omics' in terms of this script, I am referring to any dataset that describes changes in the abundance (be that: expression, reads, counts, ect.), of unique species (RNA, proteins, metabolites, glycans, ect.).

At present, the script produces a network graph of co-expression relationships between 'molecular species' among different -omics datasets. It is capable of handling time-course datasets. I have developed this for examination of co-regulatory relationships during circadian rhythms, though any regular-interval timeseries is suitable. The approach for this is:

After z-score normalisation, euclidean distances are determined between all molecular species at all timepoints. (This can be very intensive- I would reccomend only inputting species-of-interest at this stage)

The average distance among all timepoints is calculated for every pair of species. (This gives a similarity-over-time)

The average distances are used as edge weights in the construction of a network graph, species with similar expression are grouped closely together within the network.

Additionally, I have included a function which performs dynamic time warping on elements of a dataset. In theory this allows for the identification of molecular species with time-delayed or idiosyncatic manifestations of the same underlying rhytmic trends. In practice, I am unsure if it is suitable for analysing data where a majority of rhythms are sinusoidal. Nevertheless, I have included it as more utility is not a bad thing!

At this stage the script does not include functionality for detection of rhythmicity, I am testing for rhythmicity using other tools and importing rhythmic species manually.

Finally, I have updated the script to use pyvis for the visualisation of network graphs. Pyvis creates a html file which allows interactive exploration of the network graphs generated from OmicsCat.

My next steps will include:

    More robust normalisation of different omics varieties, before concatenation. One thing I have noticed is that the majority of low-weight (closely co-regulated) interactions are occuring within members of different 'omics-levels'. While this seems biologically intuative (i.e: it's no surpirse similar metabolites are often found together), I do wonder if extra steps can be taken to standardise the different levels. I have no specific strategies in mind at the moment but I will consider this an ongoing area for improvement.

    Following the implementation of proper visualisation tools, the next step is to develop the script into functions that can be utilised as a package. In doing so, I hope the expand on the visualisation to incorporate metadata (colour nodes by omics-level, for example).

    Implementation of algorithms for detecting rhythmicity- this may be an inevitable inclusion.

The current and proposed workflows are illustrated in the figure below.

![image](https://github.com/Alex-RW-Bennett/OmicsCat/assets/131264603/90c3f085-3320-4014-9228-50f96665837b)

As a final note, I should mention I have used data published by Bignon _et al._ during my development/testing. Thanks! (PMID: 36862511)
 
