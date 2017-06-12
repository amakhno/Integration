reset

# eps
set encoding koi8r
set terminal postscript eps size 3.5,2.62 enhanced color font 'Helvetica,15' lw 1
set output 'roots1.eps'
set xzeroaxis
# Line styles
set border linewidth 1.5
# Legend
#set key at 6,1.3
# Axes label 
set xlabel '{/Helvetica-Oblique ๗าลอั, t}'
set ylabel "{/Helvetica-Oblique ๗าลอั, t\'}"
# Axis ranges
set xrange[0:30]
set yrange [-8: 4.2]
# Axis labels
#set ytics 1
#set tics scale 0.75
# Functions to plot
a = 0.9
f(x) = a * sin(x)
g(x) = a * cos(x)
# Plot
plot 	'roots.txt' using ($1):($2) notitle w p pt 7 ps 0.5\


