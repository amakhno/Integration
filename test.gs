set terminal wxt size 800,800 enhanced font "Helvetica,20" persist

plot 'analytical.txt' using ($1):($2) w l ls 1, 'numerical.txt' using ($1):($2) w l
