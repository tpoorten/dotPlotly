# dotPlotly
Create an interactive dot plot from mummer output

R script that makes a plotly interactive dot plot. The input is the coords file from the MUMmer utility program `show-coords -c`.

Example:
```
show-coords -c example.delta > example.coords
./mummerCoordsDotPlotly.R -i example.coords -o out -s -t -m 500 -q 500000 -k 7 -l
```

Updated (10/29/17) :  add -t param, add handling of one ref chrom


The script requires three R packages: `install.packages(c("optparse", "ggplot2", "plotly"))`. 
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

	-h, --help
		Show this help message and exit
```