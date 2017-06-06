reset

# eps
set encoding koi8r
set terminal postscript eps size 3.5,2.62 enhanced color font 'Helvetica,15' lw 1
set output 'terminal_roots_eps.eps'
set title "Корни уравнения"
set xzeroaxis lw 1.5
# Line styles
set border linewidth 1.5
# Legend
#set key at 6,1.3
# Axes label 
set xlabel '{/Helvetica-Oblique Время, t}'
# Axis ranges
set xrange[0:70]
set yrange[-5:5]
# Axis labels
#set ytics 1
# Plot
plot 'roots.txt' notitle w p pt 7 ps 0.5

