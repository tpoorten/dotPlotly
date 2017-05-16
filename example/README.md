# dotPlotly example
Compare diploid [Brassica rapa](https://genomevolution.org/coge/GenomeInfo.pl?gid=28890) and allo-tetraploid [Brassica napus](https://genomevolution.org/coge/GenomeInfo.pl?gid=25695)
(hybrid of B. rapa and B. oleracea)

Mummer commands:

```
nucmer --maxmatch -l 100 -c 500 Brassica_rapa.faa Brassica_napus_rape.faa -p Brapa_Bnapus.nucmer
show-coords -c Brapa_Bnapus.nucmer.delta > Brapa_Bnapus.nucmer.coords
```

Make dot plot:


```
../mummerCoordsDotPlotly.R -i Brapa_Bnapus.nucmer.coords -o Brapa_Bnapus.nucmer -p 8 -q 100000 -m 500 -l -s -k 10
```
