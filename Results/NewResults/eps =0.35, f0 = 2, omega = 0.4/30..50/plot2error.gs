reset

# eps
set encoding koi8r
set terminal postscript eps size 3.5,2.62 enhanced color font 'Helvetica,15' lw 1
set output 'terminal_error_eps.eps'
set title 'Относительная ошибка'
# Line styles
set border linewidth 1.5
set style line 1 linecolor rgb '#0060ad' linetype 0 linewidth 2   # blue
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 2  # red
# Legend
#set key at 6,1.3
# Axes label 
set xlabel '{/Helvetica-Oblique Время, t}'
set ylabel '{/Helvetica-Oblique M (eps, t)}'
# Axis ranges
set xrange[5:35]
#set yrange[-1.5:1.5]
# Axis labels
#set ytics 1
#set tics scale 0.75
# Plot
plot 	'error.txt' using ($1):(abs($2-$3)/$3) notitle with lines ls 2

