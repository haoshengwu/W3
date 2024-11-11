in gunplot

set view map
set xrage[1.8:3.2]
splot 'Bfield_in_Pol' using 1:2:3 with image
replot 'sep1' using ($1/1000):($2:1000):(0)

