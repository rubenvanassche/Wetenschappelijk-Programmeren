set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'fit.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "fittedPointsN6.dat" using 1:2 with lines, \
"fittedPointsN7.dat" using 1:2 with lines, \
"fittedPointsN8.dat" using 1:2 with lines, \
"fittedPointsN11.dat" using 1:2 with lines, \
"fittedPointsN13.dat" using 1:2 with lines, \
"fittedPointsN15.dat" using 1:2 with lines, \
"points.dat" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'fitlow.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "fittedPointsN6.dat" using 1:2 with lines, \
"fittedPointsN7.dat" using 1:2 with lines, \
"fittedPointsN8.dat" using 1:2 with lines, \
"points.dat" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'fithigh.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "fittedPointsN11.dat" using 1:2 with lines, \
"fittedPointsN10.dat" using 1:2 with lines, \
"fittedPointsN12.dat" using 1:2 with lines, \
"fittedPointsN14.dat" using 1:2 with lines, \
"fittedPointsN13.dat" using 1:2 with lines, \
"fittedPointsN15.dat" using 1:2 with lines, \
"points.dat" using 1:2 with points
