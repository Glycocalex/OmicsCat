# OmicsCat
Concatenation integration for generic multi-omics datasets.

Currently this is just a wip script, however I intend to develop this into a simple package.
Though robust tools for the concatenation of generic -omics data are available, as far as I am aware no such tools have been developed with a focus on circadian rhythms. As well as good python-practice for myself, I am hoping this will eventually be a useful, quick-and-easy tool for basic multi-omics integration for data collected over the course of biological rhythms.

As for the definition of 'Generic-omics' in terms of this script, I am referring to any dataset that describes changes in the abundance (be that: expression, reads, counts, ect.), of unique species (RNA, proteins, metabolites, glycans, ect.).

At present, the script produces a network graph of co-expression relationships between 'molecular species' among different -omics datasets.
It is capable of handling time-course datasets. I have developed this for examination of co-regulatory relationships during circadian rhythms.
The approach for this is:
1) After z-score normalisation, euclidean distances are determined between all molecular species at all timepoints. (This can be very intensive- I would reccomend only inputting species-of-interest at this stage)
2) The average distance among all timepoints is calculated for every pair of species. (This gives a similarity-over-time)
3) The average distances are used as edge weights in the construction of a network graph, species with similar expression are grouped closely together within the network.

At this stage the script does not include functionality for detection of rhythmicity, I am testing for rhythmicity using other tools and importing rhythmic species manually.
Incorporating rhytmicity detection may be beyond the scope of this script, we will see how things develop.

My next steps will include:
More robust normalisation of different omics varieties, before concatenation.
Implementation of improved network visualisation. (Currently the base visualisation features of networkX are being used)
