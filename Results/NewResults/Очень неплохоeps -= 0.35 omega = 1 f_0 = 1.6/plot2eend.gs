reset

# eps
set encoding koi8r
set terminal postscript eps size 3.5,2.62 enhanced color font 'Helvetica,15' lw 1
set output 'end2.eps'
set xzeroaxis
# Line styles
set border linewidth 1.5
set style line 1 linecolor rgb '#0060ad' linetype 0 linewidth 2   # blue
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 2  # red
# Legend
#set key at 6,1.3
# Axes label 
set xlabel '{/Helvetica-Oblique �����, t}'
set ylabel '{M (eps, t)}'
# Axis ranges
set xrange[20:70]
# Axis labels
#set ytics 1
set tics scale 0.75
# Functions to plot
a = 0.9
f(x) = a * sin(x)
g(x) = a * cos(x)
# Plot
plot 	'analytical.txt' using ($1):($2) title '������������� �������' with lines ls 1, \
		'numerical.txt' using ($1):($2) title '��������� �������' with lines ls 2


