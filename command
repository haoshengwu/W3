in gunplot

set view map
set xrange[1.2:3.2]
splot 'Bfield_in_Pol' using 1:2:3 with image
replot 'sep1' using ($1/1000):($2/1000):(0)

#plot line
plot 'euler_line_tracing' using ($1):($3)
splot 'euler_line_tracing_xyz'
replot 'sep1' using (-$1/1000):(0):($2/1000)
replot 'sep1' using ($1/1000):(0):($2/1000)