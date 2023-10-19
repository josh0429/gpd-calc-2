set logscale x
set logscale y

set format x "10^{%T}"
set format y "10^{%T}"

set yrange[0.001:10000]
set xrange[1:250]

plot "cff_zp07lnq2.dat" using 1:($3*4) w l lw 2 lt 2 notitle,\
     "cff_zp15lnq2.dat" using 1:($3*.8) w l lw 2 lt 6 notitle,\
     "cff_zp25lnq2.dat" using 1:($3*.25) w l lw 2 lt 7 notitle,\
     "cff_zp35lnq2.dat" using 1:($3*.15) w l lw 2 lt 1 notitle,\
     "cff_zp001lnq2.dat" using 1:($3*7.5) w l lw 2 lt 4 notitle,\
     "cff_zp57lnq2.dat" using 1:($3*.05) w l lw 2 lt 8 notitle
     