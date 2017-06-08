reset

# eps
set encoding koi8r
set terminal postscript eps size 3.5,2.62 enhanced color font 'Helvetica,15' lw 1
set output 'terminal_roots_eps.eps'
set title 'Сравнение'
set xzeroaxis lw 1.5
# Line styles
set border linewidth 1.5
set style line 1 lc rgb '#0060ad' pt 3
# Legend
#set key at 6,1.3
# Axes label 
set xlabel '{/Helvetica-Oblique Время, t}'
set ylabel '{/Helvetica-Oblique M (eps, t)}'
# Axis ranges
set xrange[0:35]
set yrange[-5:5]
# Axis labels
#set ytics 1
set tics scale 0.75
# Functions to plot
a = 0.9
f(x) = a * sin(x)
g(x) = a * cos(x)
# Plot
plot 'roots.txt' notitle 'Корень' ls 1

