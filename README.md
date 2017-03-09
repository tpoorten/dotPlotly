# dotPlotly
Create an interactive dot plot from mummer output

R script make a dot plot from the coords file from MUMmer utility program output of `show-coords -c`.

Example:
```
show-coords -c example.delta > example.coords
./mummerCoordsDotPlotly.R -i example.coords -k 4 -m 10000 -o out
```

