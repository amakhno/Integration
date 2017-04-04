set terminal wxt size 800,800 enhanced font "Helvetica,20" persist
set output 'output.png'
set style line 1 lt 2 lc rgb '#007dd2bd' # blue


plot 'F1.txt' ls 1 with lines, 'S.txt'  
