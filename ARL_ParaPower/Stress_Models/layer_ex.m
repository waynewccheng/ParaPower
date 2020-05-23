function sigma=layer_ex(row,col,lay,dT,dz,cte,E,nu,nlsub,Mat,NL)
% This function calculates the thermal stress in the film layers
% Based on paper by C. H. Hsueh, Thin Solid Films, Vol 418, 2002
% global dz cte E nu nlsub Mat NL
if lay == nlsub+1
    z=dz(lay)/2;
else
    z=sum(dz(nlsub+1:lay-1))+dz(lay)/2;
end
%fprintf(' MAT(%.0f, %.0f, %.0f) = %.0f\n',row,col,nlsub, Mat(row,col,nlsub))

Ebs=E(Mat(row,col,nlsub))/(1-nu(Mat(row,col,nlsub)));
dTs=mean(dT(row,col,1:nlsub));
ts=sum(dz(1:nlsub));
Ebl=E(Mat(row,col,lay))/(1-nu(Mat(row,col,lay)));
% Calculate sums in the equations for c, tb and r
sumcn=0;
sumd=0;
sumtbn=0;
sumrn=0;
sumrd=0;
for i=nlsub+1:NL
    if Mat(row,col,i) == 0
        sumcn=sumcn+0;
        sumd=sumd+0;
        sumtbn=sumtbn+0;
    else
        Ebi=E(Mat(row,col,i))/(1-nu(Mat(row,col,i)));
        sumcn=sumcn+Ebi*dz(i)*cte(Mat(row,col,i))*dT(row,col,i);
        sumd=sumd+Ebi*dz(i);
        if i == nlsub+1
            him1=0;
        else
            him1=sum(dz(nlsub+1:i-1));
        end
        sumtbn=sumtbn+Ebi*dz(i)*(2*him1+dz(i));
    end
end
c=(Ebs*ts*cte(Mat(row,col,nlsub))*dTs+sumcn)/(Ebs*ts+sumd);
% TC eq. 6 (c = uniform strain component)
tb=((-Ebs*ts^2)+sumtbn)/(2*(Ebs*ts+sumd));
% TC eq. 8 (tb = location of the bending axis)
for i=nlsub+1:NL
    if Mat(row,col,i) == 0
        sumrn=sumrn+0;
        sumrd=sumrd+0;
    else
        Ebi=E(Mat(row,col,i))/(1-nu(Mat(row,col,i)));
        if i == nlsub+1
            him1=0;
        else
            him1=sum(dz(nlsub+1:i-1));
        end
        sumrd=sumrd+Ebi*dz(i)*(c-cte(Mat(row,col,i))*dT(row,col,i))*(2*him1+dz(i));
        sumrn=sumrn+Ebi*dz(i)*(6*(him1^2+him1*dz(i))+2*dz(i)^2-3*tb*(2*him1+dz(i)));
    end
end
r=((Ebs*ts^2)*(2*ts+3*tb)+sumrn)/(3*(Ebs*(c-cte(Mat(row,col,nlsub))*dTs)*ts^2-sumrd));
% TC eq. 10 (r = radius of curvature of system)
% found r by inversing equation for curvature (1/r)
eps=c+(z-tb)/r;
% TC eq. 2 (total strain in multilayer)
sigma=Ebl*(eps-cte(Mat(row,col,lay))*dT(row,col,lay));
% TC eq. 3 (normal stresses in substrate - also found in ppt. p. 4)
% TC Ebl = modulus of elasticity (Young's modulus)
% TC cte = coefficients of thermal expansion
% TC dt = cooling temperature range
end
