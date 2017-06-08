reset

# eps
set encoding koi8r
set terminal postscript eps size 3.5,2.62 enhanced color font 'Helvetica,15' lw 1
set output 'terminal_eps.eps'
set xzeroaxis lw 1.5
# Line styles
set border linewidth 1.5
set style line 1 lc rgb '#0060ad' pt 5 ps 0.5
set style line 2 lc rgb '#dd181f' pt 7 ps 0.5
set style line 3 lc rgb '#000000' pt 3	 ps 0.5
# Legend
set key at 6,1.3
# Axes label 
set xlabel '{/Helvetica-Oblique Время, t}'
# Axis ranges
set xrange[0:20]
set yrange[-5:5]
# Axis labels
#set ytics 1
set tics scale 0.75
# Functions to plot
# Plot
plot 'roots1.txt' notitle 'Корень' ls 1, \
		'roots2.txt' notitle 'Корень' ls 2, \
		'roots3.txt' notitle 'Корень' ls 3, \

