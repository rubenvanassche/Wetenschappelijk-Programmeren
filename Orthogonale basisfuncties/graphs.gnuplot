set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/polynomial.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/polynomial-function.data" using 1:2 with lines, \
"data/polynomial-interpolation.data" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/polynomial-error.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1e-3
plot "data/polynomial-error.data" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/polynomial-error2.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1e-7
plot "data/polynomial-error2.data" using 1:2 with lines


set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/polynomial-approx.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set boxwidth 0.1
set style fill solid
set xtics 1e-1
set ytics 1e-16
plot "data/polynomial-approx.data" using 1:2:xtic(1) with boxes

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/mockcheb.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/mockcheb-function.data" using 1:2 with lines, \
"data/mockcheb-interpolation.data" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/mockcheb-error.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1e-3
plot "data/mockcheb-error.data" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/mockcheb-error2.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1e-7
plot "data/mockcheb-error2.data" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/mockcheb-approx.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set boxwidth 0.1
set style fill solid
set xtics 1
set ytics 1e-16
plot "data/mockcheb-approx.data" using 1:2:xtic(1) with boxes

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/mockcheb-gridPoints.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set style fill solid
set yrange [-1: 1]
plot "data/mockcheb-grid.data" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/mockcheb-chebPoints.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set style fill solid
set yrange [-1: 1]
plot "data/mockcheb-cheb.data" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/spline.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/spline-function.data" using 1:2 with lines, \
"data/spline-interpolation.data" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/spline-error.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/spline-error.data" using 1:2 with lines

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/spline-approx.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set boxwidth 0.1
set style fill solid
set xtics 1e-1
set ytics 1
set yrange[-1:1]
plot "data/spline-approx.data" using 1:2:xtic(1) with boxes