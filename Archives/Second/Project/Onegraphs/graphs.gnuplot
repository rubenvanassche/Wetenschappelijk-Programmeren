set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'Sin-Plot.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 0.1
set ytics 0.1
plot "Sin-Plot.dat" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'Polynomial-Interpolated.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 0.1
set ytics 0.1
set key autotitle columnhead
plot "Polynomial-interpolated.dat" using 1:2 with lines


set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'Polynomial-Error.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 0.1
set ytics 0.000000000001
plot "Polynomial-error.dat" using 1:2 with lines


set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'Spline-Interpolated.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 0.1
set ytics 0.1
set key autotitle columnhead
plot "Spline-interpolated.dat" using 1:2 with lines


set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'Spline-Error.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 0.1
set ytics 0.00000000001
plot "Spline-error.dat" using 1:2 with lines


set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'LeastSquares-Interpolated.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 0.1
set ytics 0.1
set key autotitle columnhead
plot "LeastSquares-interpolated.dat" using 1:2 with lines


set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'LeastSquares-Error.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 0.1
set ytics 0.001
plot "LeastSquares-error.dat" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'All.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 0.1
set ytics 0.1
set key autotitle columnhead
plot "LeastSquares-interpolated.dat" using 1:2 with lines title "Least Squares", \
"Polynomial-interpolated.dat" using 1:2 with lines title "Polynomial", \
"Sin-Plot.dat" using 1:2 with lines  title "sin(x)", \
"Spline-interpolated.dat" using 1:2 with lines title "Natural Cubic Spline"