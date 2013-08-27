reset

file="xxx"

##x=m0^2, y=lambda
Nf=1
const(x,y)=-0.5*(8.0+x)/(8.0*Nf*y)
root(x,y)=sqrt(const(x,y)*const(x,y) + 1.0/(8.0*Nf*y))
easykappa(x)=1/(x+8.0)
kappa_min(x,y) = (y==0)?easykappa(x): const(x,y) - root(x,y)
kappa_pl(x,y)  = (y==0)?easykappa(x): const(x,y) + root(x,y)


set xlabel "cutoff in GeV"
set ylabel "Higgs boson mass in GeV"

set key top left reverse
set pointsize 0.5
set title "CEP in broken phase"

set xrange[0:5000]
set yrange[0:200]

# # # # # using (246/$6):(sqrt($7)*246/$6)
plot file using (246/$6):(sqrt($7)*246/$6) with points pt 6 lc rgbcolor "blue" title "naive"
     
     