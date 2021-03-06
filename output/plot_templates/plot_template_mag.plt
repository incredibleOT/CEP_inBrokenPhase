reset

file="xxx"


##x=m0^2, y=lambda
Nf=1
const(x,y)=-0.5*(8.0+x)/(8.0*Nf*y)
root(x,y)=sqrt(const(x,y)*const(x,y) + 1.0/(8.0*Nf*y))
easykappa(x)=1/(x+8.0)
kappa_min(x,y) = (y==0)?easykappa(x): const(x,y) - root(x,y)
kappa_pl(x,y)  = (y==0)?easykappa(x): const(x,y) + root(x,y)


set xlabel "kappa"
set ylabel "magnetization"

set key top left reverse
set pointsize 0.5
set title "CEP in broken phase"

plot file using (kappa_min($1,$2)):($6/sqrt(2.0*kappa_min($1,$2))) with points pt 6 lc rgbcolor "blue" title "l=0.0" 