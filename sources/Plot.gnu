# ===========================================================================================================
# réinitialisation des paramètres
reset
# fichier de sortie
set term postscript eps size 3.5,2.62 enhanced color font 'Helvetica,12'
set output "../graphes/pave.eps"
set encoding utf8

# paramètres
set title "Test\nparticules"
set grid
set xlabel "x"
set ylabel "y"
set size ratio -1 # repère orthonormé
set xrange[-0.2:1.2]
set yrange[-0.2:1.2]

# tracé
plot "../sorties/x.dat" u 1:2 lc rgb "#008000" lw 3 title "particules",\
     "../sorties/nvec.dat" u 1:2:3:4 with vectors head filled lt rgb "red"
#    "../sorties/grad.dat" u 1:2:3:4 with vectors head filled lt rgb "green"



# ===========================================================================================================
reset
set term postscript eps size 3.5,2.62 enhanced color font 'Helvetica,12'
set output "../graphes/circle.eps"
set encoding utf8

# paramètres
set title "Test\nparticules"
set grid
set xlabel "x"
set ylabel "y"
set size ratio -1 # repère orthonormé
#set xrange[-0.2:1.2]
#set yrange[-0.2:1.2]
set border 16

# tracé
plot "../sorties/x_circle.dat" u 1:2 lc rgb "#008000" lw 1 title "particules"



# ===========================================================================================================
# affichage écran
set term x11 enhanced
replot
