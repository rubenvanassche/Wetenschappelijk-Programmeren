set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'Polynomial.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 10
set ytics 1
plot "Polynomial.dat" using 1:2 with lines, \
"Points.dat" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'naturalSpline.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 10
set ytics 1
plot "naturalSpline.dat" using 1:2 with lines, \
"Points.dat" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'CubicSpline.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 10
set ytics 1
plot "CubicSpline.dat" using 1:2 with lines, \
"Points.dat" using 1:2 with points
