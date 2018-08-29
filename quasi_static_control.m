function [Err_RMS,Err_RMS_ratio,Voltage]=quasi_static_control(Jacobian,coord_x,coord_y,Sensor_Parameter,Sensor_matrix,Nb_Zernike,text_options,plot_options,Geo_Parameter,Nb_electrode,G_matrix)

Pupil=Sensor_Parameter.pupil;

if(text_options==1)
    
    try
        Word=actxGetRunningServer('Word.Application');
    catch
        Word=actxserver('Word.Application'); 
    end

    set(Word,'Visible',1);

    Document=Word.Documents.Add;

    Document.PageSetup.TopMargin=10;
    Document.PageSetup.BottomMargin=10;
    Document.PageSetup.LeftMargin=10;
    Document.PageSetup.RightMargin=10;

    Content=Document.Content; 

    Content.Start=0;

    Tab=Document.Tables.Add(Content,2*(Nb_Zernike-3+1),4);

    Tab.Borders.InsideLineStyle='wdLineStyleSingle';
    Tab.Borders.InsideLineWidth='wdLineWidth025pt';
    Tab.Borders.OutsideLineStyle='wdLineStyleSingle';
    Tab.Borders.OutsideLineWidth='wdLineWidth150pt';
    
end

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
    
elseif(strcmp(Sensor_Parameter.type,'Local_slope'))
    
    W_jacobian=Sensor_matrix*Jacobian;
    
    nb_active_lenslet=size(W_jacobian,1)/2;
    
    W_jacobian=[W_jacobian,zeros(size(W_jacobian,1),2)];
    
    max_slope_x=max(max(abs(W_jacobian(1:nb_active_lenslet,:))));
    max_slope_y=max(max(abs(W_jacobian(nb_active_lenslet+1:2*nb_active_lenslet,:))));
    
    W_jacobian(:,end-1)=max_slope_x*[ones(size(W_jacobian,1)/2,1);zeros(size(W_jacobian,1)/2,1)];   % Tilt x
    W_jacobian(:,end)=max_slope_y*[zeros(size(W_jacobian,1)/2,1);ones(size(W_jacobian,1)/2,1)];     % Tilt y
    
end

W_jacobian_inv=(W_jacobian'*W_jacobian)\W_jacobian';
% W_jacobian_inv_1=(Jacobian'*Jacobian)\Jacobian';
%% Deposit the Zernike mode onto the coordinate mesh

dimension=length(coord_x);

rho_c=zeros(dimension,1);
theta_c=zeros(dimension,1);
for ii=1:dimension
    rho_c(ii)=sqrt(coord_x(ii)^2+coord_y(ii)^2);
    theta_c(ii)=atan2(coord_y(ii),coord_x(ii));
end
rho_c=rho_c/(Pupil/2);

for j=2:Nb_Zernike+1
    Z_temp=sqrt(pi)*eva_pol_zer(j-1,rho_c,theta_c);
    Wtarget(:,j-1)=Z_temp;
end
clear Z_temp;
clear Z_origin_height;

for j=1:Nb_Zernike
    Wtarget(:,j)=1e-6*(Wtarget(:,j)/(max(Wtarget(:,j))-min(Wtarget(:,j))));
end

Wtarget_cal_def=zeros(length(x_notNaN),size(Wtarget,2));
for ii=1:size(Wtarget,2)
    F=scatteredInterpolant(coord_x,coord_y,Wtarget(:,ii),'natural','nearest');
    Wtarget_cal_def(:,ii)=F(x_notNaN,y_notNaN);
end

if(strcmp(Sensor_Parameter.type,'Displacement_def'))
    Wtarget_cal=Wtarget_cal_def;
elseif(strcmp(Sensor_Parameter.type,'Local_slope'))
    Wtarget_cal=Sensor_matrix*Wtarget;
end

voltages=W_jacobian_inv*Wtarget_cal;

if(strcmp(Sensor_Parameter.type,'Displacement_def'))
    voltages=voltages(1:end-3,:);
elseif(strcmp(Sensor_Parameter.type,'Local_slope'))
    voltages=voltages(1:end-2,:);
end

Voltage=zeros(Nb_electrode,Nb_Zernike);

Err_RMS=zeros(Nb_Zernike,1);
Err_RMS_ratio=zeros(Nb_Zernike,1);

for ii=3:Nb_Zernike
    
    if(text_options==1)
        Tab.Cell(2*(ii-2)-1,1).Merge(Tab.Cell(2*(ii-2)-1,2));
        Tab.Cell(2*(ii-2)-1,1).Merge(Tab.Cell(2*(ii-2)-1,2));
        Tab.Cell(2*(ii-2)-1,1).Merge(Tab.Cell(2*(ii-2)-1,2));
    end
        
        
    if(plot_options)
        figure; 
    end
        
	[PV_target,RMS_target]=surface_plot(coord_x_cal,coord_y_cal,Wtarget_cal_def(:,ii),Pupil,plot_options,1,[0,90]);
	colorbar;
	box on;
    
    
	if(text_options==1)
        hgexport(gcf,'-clipboard');
        Tab.Cell(2*(ii-2),1).Range.Paste;
    end
           
	if(plot_options)
        figure; 
    end
        
	[PV_result,RMS_result]=surface_plot(coord_x_cal,coord_y_cal,W_jacobian_def*voltages(:,ii),Pupil,plot_options,1,[0,90]);
	colorbar;
	box on;

	if(text_options==1)
        hgexport(gcf,'-clipboard');
        Tab.Cell(2*(ii-2),2).Range.Paste;
    end

	if(plot_options)
        figure; 
    end
    
	[PV_error,RMS_error]=surface_plot(coord_x_cal,coord_y_cal,Wtarget_cal_def(:,ii)-W_jacobian_def*voltages(:,ii),Pupil,plot_options,1,[0,90]);
	colorbar;
	box on;

	if(text_options==1)
        hgexport(gcf,'-clipboard');
        Tab.Cell(2*(ii-2),3).Range.Paste;
	end

  
    voltages_zernike=voltages(:,ii)/PV_target*1e-6;
    
    voltage_range=max(voltages_zernike)-min(voltages_zernike);
    
    if(Geo_Parameter.Ring==1)
        voltages_zernike=[voltages_zernike;zeros(Geo_Parameter.Ring_dim(2,1),1)];
    end
    
    v_min=min(voltages_zernike);
    v_max=max(voltages_zernike);
    
    if(plot_options)
        figure;
        voltage_map(voltages_zernike,Nb_electrode,v_min,v_max,G_matrix);
        axis equal;
    end
    
    if(text_options==1)
        hgexport(gcf, '-clipboard');
        Tab.Cell(2*(ii-2),4).Range.Paste;    
        Tab.Cell(2*(ii-2)-1,1).Range.Text=['Zernike mode ',num2str(ii),': Err. RMS Error=',num2str(RMS_error/PV_target*1e-6*1e9,3),'nm, Voltage range=',num2str(voltage_range,4),'V'];
    end
    
    Err_RMS(ii)=RMS_error/PV_target*1e-6*1e9;
    Err_RMS_ratio(ii)=RMS_error/RMS_target;
    
    Voltage(:,ii)=voltages_zernike;
    
    close all;

end