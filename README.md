# OmicsCat
Concatenation integration for generic multi-omics datasets.

Currently this is just a wip script, however I intend to develop this into a simple package.
Though robust tools for the concatenation of generic -omics data is available, as far as I am aware no such tools have been developed with a particular focus on circadian rhythms. As well as good python-practice for myself I am hoping this will eventually be a useful, quick-and-easy tool for basic multi-omics integration for data collected over the course of biological rhythms.

At present, the script produces a network graph of co-expression relationships between 'molecular species' among different -omics datasets.
It is capable of handling time-course datasets. I have developed this for examination of co-regulatory relationships during circadian rhythms.

At this stage the script does not include functionality for detecion of rhythmicity, I am testing for rhythmicity using other tools and importing rhythmic species manually.

Incorporating rhytmicity detetcion may be beyond the scope of this script, we will see how things develop.
