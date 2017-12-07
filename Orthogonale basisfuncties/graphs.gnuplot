set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/original.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/original.dat" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/interpolated.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/interpolated.dat" using 1:2 with lines, \
"data/points.dat" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/combined.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/interpolated.dat" using 1:2 with lines, \
"data/original.dat" using 1:2 with lines, \
"data/points.dat" using 1:2 with points