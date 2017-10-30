# dotPlotly example
Compare diploid [Brassica rapa](https://genomevolution.org/coge/GenomeInfo.pl?gid=28890) and allo-tetraploid [Brassica napus](https://genomevolution.org/coge/GenomeInfo.pl?gid=25695)
(hybrid of B. rapa and B. oleracea)

### 1. Mummer example

```
nucmer --maxmatch -l 80 -c 100 Brassica_rapa.faa Brassica_napus_rape.faa -p Brapa_Bnapus.nucmer
delta-filter -r Brapa_Bnapus.nucmer.delta > Brapa_Bnapus.nucmer.delta.filter
show-coords -c Brapa_Bnapus.nucmer.delta.filter > Brapa_Bnapus.nucmer.delta.filter.coords
```

Make dot plot:


```
../mummerCoordsDotPlotly.R -i Brapa_Bnapus.nucmer.delta.filter.coords -o Brapa_Bnapus.nucmer.plot -s -m 500 -q 500000 -k 10 -x -t -l
```

![Dotplot](https://github.com/tpoorten/dotPlotly/tree/master/example/Brapa_Bnapus.nucmer.plot.png)

### 2. Minimap2 example 
(16 seconds on a Macbook Pro)

```
minimap2 -x asm10 Brassica_rapa.faa Brassica_napus_rape.faa > Brapa_Bnapus.minimap2.paf
```

Make dot plot:

```
../pafCoordsDotPlotly.R -i Brapa_Bnapus.minimap2.paf -o Brapa_Bnapus.minimap2.plot -s -m 500 -q 500000 -k 10 -x -t -l
```

![Dotplot](https://github.com/tpoorten/dotPlotly/tree/master/example/Brapa_Bnapus.minimap2.plot.png)
