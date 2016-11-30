set terminal png size 800,800 enhanced font "Helvetica,20"
set output 'output.png'
set style line 1 lt 1 lc rgb '#007dd2bd' # blue

plot 'S.txt' ls 1 with lines
