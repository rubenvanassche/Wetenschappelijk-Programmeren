set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/combined.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/original.dat" using 1:2 with points, \
"data/noisy.dat" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/original.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/original.dat" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/noisy.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/noisy.dat" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/combinedLS.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/OLS.dat" using 1:2 with lines, \
"data/MnoisyLS.dat" using 1:2 with lines, \
"data/LnoisyLS.dat" using 1:2 with lines, \
"data/CnoisyLS.dat" using 1:2 with lines, \
"data/noisy.dat" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/monomialsCombinedLS.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/OLS.dat" using 1:2 with lines, \
"data/MnoisyLS.dat" using 1:2 with lines, \
"data/noisy.dat" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/OLS.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/OLS.dat" using 1:2 with lines, \
"data/original.dat" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/MnoisyLS.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/MnoisyLS.dat" using 1:2 with lines, \
"data/noisy.dat" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/LnoisyLS.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/LnoisyLS.dat" using 1:2 with lines, \
"data/noisy.dat" using 1:2 with points

set terminal png  transparent enhanced font "arial,10" fontscale 1.0 size 800, 800
set output 'plots/CnoisyLS.png'
set key bmargin left horizontal Right noreverse enhanced autotitles box linetype -1 linewidth 1.000
set samples 10000
set xtics 1
set ytics 1
plot "data/CnoisyLS.dat" using 1:2 with lines, \
"data/noisy.dat" using 1:2 with points