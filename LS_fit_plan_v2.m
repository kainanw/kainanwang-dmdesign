function [coefs]=LS_fit_plan_v2(x,y,W)
%
%   To be used to avoid problems with NaNs.
%
%   cfr Eberly_2008
%   ( http://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf )
%
%
%
%   Example of reconstruction:
%
%   x_temp=[-R:R/10:R];
%   y_temp=[-R:R/10:R];
%   [x_plan,y_plan]=meshgrid(x_temp,y_temp);
%   W_plan=zeros(size(x_plan));
%   for tt=1:1:size(x_plan,1)
%       for kk=1:1:size(x_plan,2)
%           W_plan(tt,kk)=coefs(1)+coefs(2)*x_plan(tt,kk)+coefs(3)*y_plan(tt,kk);
%       end
%   end
%
%   Or more simply: W_plan=coefs(1)+coefs(2)*x_plan+coefs(3)*y_plan;
%
% 
% Transformation to Zernike coefficients:
%
%   - Standard convention:
%       coef_piston=coefs(1);
%       coef_tilt_x=coefs(2)*R/2;
%       coef_tilt_y=coefs(3)*R/2;
%
%   - Modified Malacara's convention:
%       coef_piston=coefs(1);
%       coef_tilt_x=coefs(3)*R/2;
%       coef_tilt_y=coefs(2)*R/2;
%
%   - Malacara's convention:
%       coef_piston=coefs(1);
%       coef_tilt_x=coefs(3)*R;
%       coef_tilt_y=coefs(2)*R;
%
% Last update: 02/03/2012

flag_NaN=sum(sum(isnan(W)));

%%% The algorithm needs x,y and W to be vectors

[s1,s2]=size(x);
if(flag_NaN==0) % No NaN

    if(s1>1 & s2>1)
        x=reshape(x,s1*s2,1);
        y=reshape(y,s1*s2,1);
        W=reshape(W,s1*s2,1);
    end

elseif(flag_NaN>0) % There are NaNs

    x_temp=zeros(s1*s2,1);
    y_temp=zeros(s1*s2,1);
    W_temp=zeros(s1*s2,1);
    compteur=0;
    for tt=1:1:s1
        for kk=1:1:s2
            if(isnan(W(tt,kk))==0)
                compteur=compteur+1;
                x_temp(compteur)=x(tt,kk);
                y_temp(compteur)=y(tt,kk);
                W_temp(compteur)=W(tt,kk);
            end
        end
    end
    
    x=x_temp(1:compteur);
    y=y_temp(1:compteur);
    W=W_temp(1:compteur);

end

J=[ sum(x.*x) sum(x.*y) sum(x)
    sum(x.*y) sum(y.*y) sum(y)
    sum(x)    sum(y)    length(x) ];

p=[ sum(x.*W)
    sum(y.*W)
    sum(W)    ];

coefs=inv(J)*p;

coefs=[coefs(3) coefs(1) coefs(2)]'; % Piston, Tilt x, Tilt y