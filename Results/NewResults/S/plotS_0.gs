reset

# eps
set encoding koi8r
set terminal postscript eps size 3.5,2.62 enhanced color font 'Helvetica,15' lw 0.5
set output 'S.eps'
# Line styles
set termoption dash
set border linewidth 1.5
set style line 1 dt 1 linecolor rgb '#0060ad' # blue
set style line 2 dt 2 linecolor rgb '#dd181f' linetype 1 linewidth 2  # red
set style line 3 dt 4 linecolor rgb '#000000' linetype 1 linewidth 2  # red
# Legend
set key at 15, -140
# Axes label 
set xlabel "{/Helvetica-Oblique ๗าลอั, t\'}"
set ylabel "S{/Helvetica-Oblique (t\', t)}"
# Axis ranges
#set xrange[20:70]
#set yrange[-1.5:1.5]
# Plot
plot	'S-20.txt' title "S{/Helvetica-Oblique (t\', -20)}" with line ls 2, \
		'S0.txt' title "S{/Helvetica-Oblique (t\', 0)}" with line ls 1, \
		'S20.txt' title "S{/Helvetica-Oblique (t\', 20)}" with line ls 3 \


