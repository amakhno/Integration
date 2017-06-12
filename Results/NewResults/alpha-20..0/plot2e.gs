reset

# eps
set encoding koi8r
set terminal postscript eps size 3.5,2.62 enhanced color font 'Helvetica,15' lw 1
set output 'alpha-20..0.eps'
set xzeroaxis
# Line styles
set border linewidth 1.5
set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 2  dt 1 # blue
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 2  dt 2 # red
set style line 3 linecolor rgb '#000000' linetype 1 linewidth 2  dt 4 # black
# Legend
set key at 19,-1
# Axes label 
set xlabel "{/Helvetica-Oblique ๗าลอั, t\'}"
set ylabel "{/Helvetica-Oblique Alpha (t\', t)}"
# Axis ranges
#set xrange[20:70]
set yrange[-1.5:1.5]
# Axis labels
#set ytics 1
set tics scale 0.75
# Functions to plot
a = 0.9
f(x) = 0.6
g(x) = -0.6
# Plot
plot 	'analytical.txt' using ($1):($2) title 'Alpha' with lines ls 1, \
		g(x) title '-sqrt(2 eps)' with lines ls 2, \
		f(x) title 'sqrt(2 eps)' with lines ls 3


