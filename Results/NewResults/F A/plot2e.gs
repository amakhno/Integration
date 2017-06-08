reset

# eps
set encoding koi8r
set terminal postscript eps size 3.5,2.62 enhanced color font 'Helvetica,15' lw 0.5
set output 'F_A.eps'
set xzeroaxis
set yzeroaxis
# Line styles
set border linewidth 1.5
set style line 1 linecolor rgb '#000000' linetype 0 linewidth 2  # blue
set style line 2 linecolor rgb '#000000' linetype 1 linewidth 2  # red
# Legend
#set key at 6,1.3
# Axes label 
set xlabel '{/Helvetica-Oblique ๗าลอั, t}'
# Axis ranges
#set xrange[30:50]
#set yrange[-1.5:1.5]
# Axis labels
#set ytics 1
set tics scale 0.75
# Functions to plot
a = 0.9
f(x) = a * sin(x)
g(x) = a * cos(x)
# Plot
plot 	'f.txt' using ($1):($2) title 'F(t)' with lines ls 2, \
		'a.txt' using ($1):($2) title 'A(t)' with lines ls 1


