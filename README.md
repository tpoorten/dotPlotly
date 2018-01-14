# dotPlotly
Create an interactive dot plot from mummer output OR PAF format

R script that makes a plotly interactive and/or static (png/pdf) dot plot.

[Shiny app available for testing](https://tom-poorten.shinyapps.io/dotplotly_shiny/)

### Example (more [here](https://github.com/tpoorten/dotPlotly/tree/master/example))

For mummer (nucmer -> show-coords) output:

```
show-coords -c example.delta > example.coords
./mummerCoordsDotPlotly.R -i example.coords -o out -s -t -m 500 -q 500000 -k 7 -l
```

For PAF format (e.g. [minimap2](https://github.com/lh3/minimap2)):

```
./pafCoordsDotPlotly.R -i example.paf -o out -s -t -m 500 -q 500000 -k 7 -l
```

### Updates (10/29/17)

  + Added script for [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md)
  + Fixes for Mummer format script: add -t param, add handling of one ref chrom 

### Dependencies

The script requires three R packages: `install.packages(c("optparse", "ggplot2", "plotly"))`. 

### Script parameters

Use `./mummerCoordsDotPlotly.R -h` to see options.

```
	-i INPUT, --input=INPUT
		coords file from mummer program 'show.coords' [default NULL]

	-o OUTPUT, --output=OUTPUT
		output filename prefix [default out]

	-v, --verbose
		Print out all parameter settings [default]

	-q MIN-QUERY-LENGTH, --min-query-length=MIN-QUERY-LENGTH
		filter queries with total alignments less than cutoff X bp [default 4e+05]

	-m MIN-ALIGNMENT-LENGTH, --min-alignment-length=MIN-ALIGNMENT-LENGTH
		filter alignments less than cutoff X bp [default 10000]

	-p PLOT-SIZE, --plot-size=PLOT-SIZE
		plot size X by X inches [default 15]

	-l, --show-horizontal-lines
		turn on horizontal lines on plot for separating scaffolds  [default FALSE]

	-k NUMBER-REF-CHROMOSOMES, --number-ref-chromosomes=NUMBER-REF-CHROMOSOMES
		number of sorted reference chromosomes to keep [default all chromosmes]

	-s, --similarity
		turn on color alignments by percent similarity [default FALSE]

	-t, --identity-on-target
		turn on calculation of % identity for on-target alignments only [default FALSE]

	-x, --interactive-plot-off
		turn off production of interactive plotly [default TRUE]

	-r REFERENCE-IDS, --reference-ids=REFERENCE-IDS
		comma-separated list of reference IDs to keep [default NULL]

	-h, --help
		Show this help message and exit
```
