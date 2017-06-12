reset

# eps
set encoding koi8r
set terminal postscript eps size 3.5,2.62 enhanced color font 'Helvetica,15' lw 0.5
set output 'I_0.eps'
set xzeroaxis
# Line styles
set termoption dash
set border linewidth 1.5
set style line 1 dt 1 linecolor rgb '#0060ad' # blue
set style line 2 dt 2 linecolor rgb '#dd181f' linetype 1 linewidth 2  # red
set style line 3 dt 4 linecolor rgb '#000000' linetype 1 linewidth 2  # red
# Legend
#set key at 15, -140
# Axes label 
set xlabel "๗าลอั, t"
set ylabel "|M_0|^2 (eps, t)"
# Axis ranges
set xrange[0:20]
set yrange[-0.1:0.3]
# Plot
plot	'I_0.txt' title '|M_0|^2 (eps, t)' with line ls 1, \


