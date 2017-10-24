set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'BasicPoints1.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "BasicPoints1.dat" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'BasicPoints2.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "BasicPoints2.dat" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'ChebyChevLeastSquares.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "ChebyChevLeastSquares.dat" using 1:2 with lines, \
"BasicPoints1.dat" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'TrioLeastSquares.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "TrioLeastSquares.dat" using 1:2 with lines, \
"BasicPoints2.dat" using 1:2 with lines