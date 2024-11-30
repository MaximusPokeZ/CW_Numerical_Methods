set terminal

set multiplot layout 1,2 title "Graphs: Velocity and Distance vs Time"

set title "Velocity vs Time"
set xlabel "Time (t)"
set ylabel "Velocity (u)"
plot "result.txt" using 2:3 with lines title "u(t)"

set title "Distance vs Time"
set xlabel "Time (t)"
set ylabel "Distance (x)"
plot "result.txt" using 2:4 with lines title "x(t)"
