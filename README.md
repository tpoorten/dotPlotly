# dotPlotly
Create an interactive dot plot from mummer output

R script that makes a plotly interactive dot plot. The input is the coords file from the MUMmer utility program `show-coords -c`.

Example:
```
show-coords -c example.delta > example.coords
./mummerCoordsDotPlotly.R -i example.coords -k 4 -m 10000 -o out
```

The script requires three R packages: `install.packages(c("optparse", "ggplot2", "plotly"))`. 
Use `./mummerCoordsDotPlotly.R -h` to see options.

```
	-i INPUT, --input=INPUT
		coords file from mummer program 'show.coords' [default NULL]

	-o OUTPUT, --output=OUTPUT
		output filename prefix [default out]

	-m MIN-ALIGNMENT-LENGTH, --min-alignment-length=MIN-ALIGNMENT-LENGTH
		filter alignments less than cutoff X bp [default 10000]

	-k NUMBER-REF-CHROMOSOMES, --number-ref-chromosomes=NUMBER-REF-CHROMOSOMES
		number of sorted reference chromosomes to keep [default all chromosmes]

	-s, --similarity
		turn on color alignments by percent similarity [default FALSE]

	-h, --help
		Show this help message and exit
```