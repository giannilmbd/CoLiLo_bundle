# ogpf libray
# Rev. 0.22 of March 9th, 2018
# Licence: MIT

# gnuplot global setting
set term X11 size 640,480 enhanced font "verdana,10" title "ogpf libray: Rev. 0.22 of March 9th, 2018"

# ogpf extra configuration
# -------------------------------------------
# color definitions
set style line 1 lc rgb "#800000" lt 1 lw 2
set style line 2 lc rgb "#ff0000" lt 1 lw 2
set style line 3 lc rgb "#ff4500" lt 1 lw 2
set style line 4 lc rgb "#ffa500" lt 1 lw 2
set style line 5 lc rgb "#006400" lt 1 lw 2
set style line 6 lc rgb "#0000ff" lt 1 lw 2
set style line 7 lc rgb "#9400d3" lt 1 lw 2

# Axes
set border linewidth 1.15
set tics nomirror

# grid
# Add light grid to plot
set style line 102 lc rgb "#d6d7d9" lt 0 lw 1
set grid back ls 102

# plot style
set style data linespoints

# -------------------------------------------

 
# options
set style data linespoints


 
# plot scale
 
# Annotation: title and labels
set title "pf \& XX Residual" font ","
set xlabel "points" font ","
set ylabel "ftest" font ","
 
# axes setting

plot "-" notitle
 -0.500000000000000       0.500000000000000     
 -0.388888888888889       0.302469135802469     
 -0.277777777777778       0.154320987654321     
 -0.166666666666667       5.555555555555557E-002
 -5.555555555555558E-002  6.172839506172845E-003
  5.555555555555552E-002  6.172839506172832E-003
  0.166666666666667       5.555555555555553E-002
  0.277777777777778       0.154320987654321     
  0.388888888888889       0.302469135802469     
  0.500000000000000       0.500000000000000     
e
