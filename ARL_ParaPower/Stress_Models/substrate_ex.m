function sigma=substrate_ex(row,col,lay,dT,dz,cte,E,nu,nlsub,Mat,NL)
% This function calculates the thermal stress in the top layer of the
% substrate material
% global dz cte E nu nlsub Mat NL
% TC calculate z
z=-sum(dz(1:nlsub))+sum(dz(1:lay))-dz(lay)/2;
Ebs=E(Mat(row,col,lay))/(1-nu(Mat(row,col,lay)));
dTs=mean(dT(row,col,1:nlsub));
ts=sum(dz(1:nlsub));
% Calculate sums in the equations for c, tb and r
sumcn=0;
% TC sumcn = summation in numerator of eq. 6 (for c)
sumd=0;
% TC sumd = summation in denominator of eq. 8 (tb) and eq. 6 (c)
sumtbn=0;
% TC sumtbn = summation in numerator of eq. 8 (for tb)
sumrn=0;
% TC sumrn = summation in numerator of eq. 10 (for r)
sumrd=0;
% TC sumrd = summation in denominator of eq. 10 (for r)

% TC The following loop is calculates sumcn, sumd, and sumtbn.
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

% TC go from 2nd layer to n
% TC This loop calculates sumrd and sumrn. 
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
r=((Ebs*ts^2)*(2*ts+3*tb)+sumrn)/(3*(Ebs*(c-cte(Mat(row,col,lay))*dTs)*ts^2-sumrd));
% TC eq. 10 (r = radius of curvature of system)
% found r by inversing equation for curvature (1/r)
eps=c+(z-tb)/r;
% TC eq. 2 (total strain in multilayer)
sigma=Ebs*(eps-cte(Mat(row,col,lay))*dT(row,col,lay));
% TC eq. 3
% TC Ebs = modulus of elasticity (Young's modulus)
% TC cte = coefficients of thermal expansion
% TC dt = cooling temperature range
end
