set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'T7.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 0.1
set ytics 0.1
plot "T7.dat" using 1:2 with lines,\
"T7_zeros.dat" using 1:2,\
"T7_extremas.dat" using 1:2
