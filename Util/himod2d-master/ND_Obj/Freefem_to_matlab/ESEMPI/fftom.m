!FreeFem++ lapl.edp
clear
addpath ./../
[points seg tri]=importfilemesh('griglia.msh');
sol=importfiledata('sol');

pdemesh(points,seg,tri);
pause
pdeplot(points,seg,tri,'xydata',sol,'zdata',sol,'mesh','on','colormap','jet');
