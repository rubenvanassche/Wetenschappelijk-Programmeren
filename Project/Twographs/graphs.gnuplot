set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'trioLS-interpolated.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 0.1
plot "trioLS-interpolated.dat" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'trioI-interpolated.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 0.1
plot "trioI-interpolated.dat" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'PolynomialRandom-interpolated.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 0.1
plot "PolynomialRandom-interpolated.dat" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'PolynomialEquidistant-interpolated.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 0.1
set ytics 0.1
plot "PolynomialEquidistant-interpolated.dat" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'PolynomialLS-interpolated.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 10
plot "PolynomialLS-interpolated.dat" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 2000, 2000
set output 'All.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 0.1
set key autotitle columnhead
plot "trioLS-interpolated.dat" using 1:2 with lines title "Triogoniometric Least Squares", \
"trioI-interpolated.dat" using 1:2 with lines title "Polynomial Triogoniometric Interpolation", \
"PolynomialRandom-interpolated.dat" using 1:2 with lines  title "Polynomial Random Interpolation", \
"PolynomialEquidistant-interpolated.dat" using 1:2 with lines title "Polynomial Equidistant Interpolation"