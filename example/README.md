# dotPlotly example
Compare diploid [Brassica rapa](https://genomevolution.org/coge/GenomeInfo.pl?gid=28890) and allo-tetraploid [Brassica napus](https://genomevolution.org/coge/GenomeInfo.pl?gid=25695)
(hybrid of B. rapa and B. oleracea)

Mummer commands:

```
nucmer --maxmatch -l 100 -c 500 Brassica_rapa.faa Brassica_napus_rape.faa -p Brapa_Bnapus.nucmer
delta-filter -r Brapa_Bnapus.nucmer.delta > Brapa_Bnapus.nucmer.delta.filter
show-coords -c Brapa_Bnapus.nucmer.delta.filter > Brapa_Bnapus.nucmer.delta.filter.coords
```

Make dot plot:


```
../mummerCoordsDotPlotly.R -i Brapa_Bnapus.nucmer.delta.filter.coords -o Brapa_Bnapus -p 8 -q 20000 -m 0 -l -s -k 10
```
