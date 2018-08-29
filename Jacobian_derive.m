function Jac=Jacobian_derive(Jacobian,coord_x,coord_y,Sensor_Parameter,Sensor_matrix)

Pupil=Sensor_Parameter.pupil;

%% Jacobian for the voltage computation

% Slope sensor: Use sensor matrix multiplied by node-displacement Jacobian
% Displacement sensor: Use manual mesh projected onto node-displacement Jacobian

N=200;
X=[-1.05*Pupil/2:1.05*Pupil/(N-1):1.05*Pupil/2];
Lx=length(X);
Y=[-1.05*Pupil/2:1.05*Pupil/(N-1):1.05*Pupil/2];
[xgrid,ygrid]=meshgrid(X,Y);

x=reshape(xgrid,length(xgrid)^2,1);
y=reshape(ygrid,length(ygrid)^2,1);
    
mask_NaN=ones(size(x));
for tt=1:length(x)
    distance_2=(x(tt)^2)+(y(tt)^2);
    if(distance_2>(Pupil/2)^2) % coef security
        mask_NaN(tt)=NaN;
    end
end
x=x.*mask_NaN;
y=y.*mask_NaN;
    
% coordinate
x_notNaN=x(isfinite(x));
y_notNaN=y(isfinite(y));


W_jacobian_def=zeros(length(x_notNaN),size(Jacobian,2));
    
for ii=1:size(Jacobian,2)
    F=scatteredInterpolant(coord_x,coord_y,Jacobian(:,ii),'natural','nearest');
    W_jacobian_def(:,ii)=F(x_notNaN,y_notNaN);
end
    
coord_x_cal=x_notNaN;
coord_y_cal=y_notNaN;
    
    
if(strcmp(Sensor_Parameter.type,'Displacement_def'))
    
    W_jacobian=W_jacobian_def;
    W_jacobian=[W_jacobian,zeros(size(W_jacobian,1),3)];
    
    max_W=max(max(abs(W_jacobian)));
    
    W_jacobian(:,end-2)=max_W*ones(size(W_jacobian,1),1);           % Piston
    W_jacobian(:,end-1)=coord_x_cal/max(abs(coord_x_cal))*max_W;            % Tilt x
    W_jacobian(:,end)=coord_y_cal/max(abs(coord_x_cal))*max_W;              % Tilt y
    
    W_jacobian_NoPTT=W_jacobian(:,1:end-3);
    
    
elseif(strcmp(Sensor_Parameter.type,'Local_slope'))
    
    W_jacobian=Sensor_matrix*Jacobian;
    
    nb_active_lenslet=size(W_jacobian,1)/2;
    
    W_jacobian=[W_jacobian,zeros(size(W_jacobian,1),2)];
    
    max_slope_x=max(max(abs(W_jacobian(1:nb_active_lenslet,:))));
    max_slope_y=max(max(abs(W_jacobian(nb_active_lenslet+1:2*nb_active_lenslet,:))));
    
    W_jacobian(:,end-1)=max_slope_x*[ones(size(W_jacobian,1)/2,1);zeros(size(W_jacobian,1)/2,1)];   % Tilt x
    W_jacobian(:,end)=max_slope_y*[zeros(size(W_jacobian,1)/2,1);ones(size(W_jacobian,1)/2,1)];     % Tilt y
    
    W_jacobian_NoPTT=W_jacobian(:,1:end-2);
    
end

% W_jacobian_inv=(W_jacobian'*W_jacobian)\W_jacobian';

Jac.N=N;
Jac.xgrid=xgrid;
Jac.ygrid=ygrid;
Jac.x=x;
Jac.y=y;
Jac.coord_x_cal=coord_x_cal;
Jac.coord_y_cal=coord_y_cal;
Jac.W_jacobian=W_jacobian;
Jac.W_jacobian_NoPTT=W_jacobian_NoPTT;
