set terminal wxt size 800,800 enhanced font "Helvetica,20" persist
set output 'output.png'
set format y "%.12f"
set ytics format "%.12f"

plot 'analytical.txt' using ($1):($2*2.5) w l ls 1, 'numerical.txt' w l
