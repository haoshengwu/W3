in gunplot

set view map
set xrange[1.2:3.2]
splot 'Bfield_in_Pol' using 1:2:3 with image
replot 'sep1' using ($1/1000):($2/1000):(0)
replot 'sep1' using ($1/1000):($2/1000):(0
#plot line
plot 'euler_line_tracing' using ($1):($3)
splot 'euler_line_tracing_xyz'
replot 'sep1' using (-$1/1000):(0):($2/1000)
replot 'sep1' using ($1/1000):(0):($2/1000)

#plot 3D lines

set vime map
splot 'sep1' using ($1/1000):(0):($2/1000) lt rgb "red" title 'sep1'
replot 'sep3' using ($1/1000):(0):($2/1000) lt rgb "red" title 'sep1'
replot 'sep1' using (-$1/1000):(0):($2/1000) lt rgb "blue" title 'sep2
replot 'sep3' using (-$1/1000):(0):($2/1000) lt rgb "blue" title 'sep2'
replot 'euler_line_tracing_xyz' title 'line1' lt rgb "black"
replot 'euler_line_tracing_xyz_line2' title 'line2' lt rgb "green"
