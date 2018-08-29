clear;
clc;
close all;

evalc(['!rm *.dat* *.log* *.dia* *.run* *.adb* *.sdb* *.sam *.spy *base* *#* *.asv* *.u18* *.in* *.out* *.res*  *.des* *.fac*']);%% Problem independents for central patch actuation test
input_material;

%% Computation setting

Filename=['DM_keystone'];
Task='IF';

Filename=[Filename,'_',Task];

%% Pattern (Geometry) parameters

Geo_Parameter.layout='Keystone';% Name of the layout
Geo_Parameter.D=0.1;% Diameter of the wafer [m]
Geo_Parameter.D_cc_region=0.08;% Diameter of the circumcircle of the electrode region [m]
Geo_Parameter.Nb_R=5;% Number of the radial order [/]
Geo_Parameter.Nb_e_R=[1,8,8,8,12];% Number of electrodes in each crown [/]
Geo_Parameter.Angle_shift_e_R=[0,20,0,0,0];% Angle of the electrode shift in each crown [/]
Geo_Parameter.W_gap=1e-3;% Width of the gap [m]
Geo_Parameter.R_e_R=[0.008,0.016,0.024,0.032,0.04];% Radius of each crown [m]
Geo_Parameter.Support_surface=[0,0,0,0,0,0,0,0,0,0];
Geo_Parameter.Sparse_factor=[0.3,0.3];

% TriSupport config
Geo_Parameter.TriSupport=0;
Geo_Parameter.TriFeet_radius=0.8*Geo_Parameter.D/2;
Geo_Parameter.TriFeet_angle=90;
Geo_Parameter.TriFeet_diameter=0.005;
Geo_Parameter.TriFeet_piezo_intersection=1;
Geo_Parameter.TriFeet_location_domain=[129,131,133];

% Ring config
Geo_Parameter.Ring=0;
    Ring_piezo_1=0.066;
    Ring_piezo_2=Geo_Parameter.D/2-Geo_Parameter.W_gap;
Geo_Parameter.Ring_dim=[Ring_piezo_1,Ring_piezo_2;6,0];

%% Mesh parameters

Mesh_Parameter.degree=1;

%% Substrate parameters

Pro_Parameter.hypothesis='MINDLIN';
Pro_Parameter.material={'SiC_alpha','PZT_PIC255','Dummy_Non'};
Pro_Parameter.sequence=[1,2,3]; % Default all the substrate
Pro_Parameter.thickness=[1000e-6,200e-6,200e-6];
Pro_Parameter.sequence_act=[3,1,2];
Pro_Parameter.sequence_deact=[1];
Pro_Parameter.sequence_shunt=[2,1,3];

%% Optical sensor parameters

% Sensor_Parameter.type='Local_slope';
Sensor_Parameter.type='Displacement_def';
Sensor_Parameter.pupil=0.05;
Sensor_Parameter.sampling=100;
Sensor_Parameter.Nb_lenslet=17;

%% Boundary condition and load parameters

Load_Parameter.BC='Isostatic';
Load_Parameter.Load=Task;

%% Samcef file generation

% Geometry matrix generation (2D)
[G_matrix,Nb_electrode]=geometry_matrix_generation_v00(Geo_Parameter);

figure;
pdegplot(G_matrix,'EdgeLabels','off','SubdomainLabels','off');
xlim([-Geo_Parameter.D/2,Geo_Parameter.D/2]);
ylim([-Geo_Parameter.D/2,Geo_Parameter.D/2]);
xlabel('X [m]');
ylabel('Y [m]');
axis square;

if(Geo_Parameter.TriSupport==1)
    
    hold on;
    TriFeet_angle=Geo_Parameter.TriFeet_angle;
    TriFeet_radius=Geo_Parameter.TriFeet_radius;
    TriFeet_diameter=Geo_Parameter.TriFeet_diameter;
    
    Trifeet_1=[TriFeet_radius*cos((TriFeet_angle)/180*pi),TriFeet_radius*sin((TriFeet_angle)/180*pi)];
    Trifeet_2=[TriFeet_radius*cos((TriFeet_angle+120)/180*pi),TriFeet_radius*sin((TriFeet_angle+120)/180*pi)];
    Trifeet_3=[TriFeet_radius*cos((TriFeet_angle+240)/180*pi),TriFeet_radius*sin((TriFeet_angle+240)/180*pi)];

    draw_feet(Trifeet_1(1),Trifeet_1(2),TriFeet_diameter/2);
    draw_feet(Trifeet_2(1),Trifeet_2(2),TriFeet_diameter/2);
    draw_feet(Trifeet_3(1),Trifeet_3(2),TriFeet_diameter/2);
    
end

% Geometry file generation (2D)
line_nb_start=500;
contour_nb_start=300;
domain_nb_start=360;

[pre_mesh_matrix,nb_structure]=geometry_matrix_translator_2d_samcef_v00([Filename,'_geo'],Geo_Parameter,G_matrix,line_nb_start,contour_nb_start,domain_nb_start,1);

% Mesh file generation (2D)
[nb_piezo_index]=mesh_generation_samcef_v00([Filename,'_mesh'],Geo_Parameter,Mesh_Parameter,Load_Parameter,pre_mesh_matrix,nb_structure);

% Property file generation
property_generation_samcef_v00([Filename,'_pro'],Geo_Parameter,Pro_Parameter,Load_Parameter,Nb_electrode,nb_piezo_index);

% Sensor file generation
active_index=sensor_generation_samcef_v00([Filename,'_sen'],Sensor_Parameter,Geo_Parameter,G_matrix,Nb_electrode,nb_piezo_index);

% Load file generation
load_generation_samcef_v00([Filename,'_load'],Geo_Parameter,Load_Parameter,Nb_electrode,nb_piezo_index);

% Samcef file assamble
preparation_samcef_v00(Filename,Geo_Parameter,G_matrix,Mesh_Parameter,Pro_Parameter,Load_Parameter);
% return;

if(Geo_Parameter.Ring==1)
    Nb_Actuator=Nb_electrode-Geo_Parameter.Ring_dim(2,1);
else
    Nb_Actuator=Nb_electrode;
end
Nb_Signal=2*length(active_index);

disp(['The number of the actuator is ',int2str(Nb_Actuator)]);
disp(['The number of the reconstruction signal is ',int2str(Nb_Signal)]);
disp(['The SAR is ',num2str(Nb_Signal/Nb_Actuator,3)]);

% clc;
% answer=questdlg('Satisfied with the configuration?', ...
% 	'DM design tool', ...
% 	'Yes','No','No');
% % Handle response
% switch answer
%     case 'Yes'
%         disp([answer ' Configuration formed'])
%         dec=1;
%     case 'No'
%         disp([answer ' Restart the design'])
%         dec=2;
% end
% 
% if(dec==2)
%     return;
% end
