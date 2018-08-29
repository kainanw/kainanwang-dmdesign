clear;
clc;
close all;

save('MATERIAL_DATABASE.mat');

%% Reflector material
input_material_property_v00('Silicon',130e9,0.275,2.59e-6,2330,0,[],[],[],[]);
input_material_property_v00('Kapton',2.8e9,0.34,20e-6,1420,0,[],[],[],[]);
input_material_property_v00('Peek',3.6e9,0.38,30e-6,1320,0,[],[],[],[]);
input_material_property_v00('SiC',410e9,0.14,4e-6,3100,0,[],[],[],'Accuratus Corp.');
input_material_property_v00('SiC_alpha',476e9,0.19,5.12e-6,3210,0,[],[],[],'Ferro-Ceramic Grinding Inc.');
% input_material_property_v00('SiC_highPR',476e9,0.34,5.12e-6,3210,0,[],[],[],'Ferro-Ceramic Grinding Inc.');
% input_material_property_v00('SiC_lowPR',476e9,0.01,5.12e-6,3210,0,[],[],[],'Ferro-Ceramic Grinding Inc.');
% input_material_property_v00('SiC_thermal',476e9,0.19,3.8e-6,3210,0,[],[],[],'Properties of Advanced Semiconductor Materials: GaN, AIN, InN, BN, SiC, SiGe');
% input_material_property_v00('SiC_3C',476e9,0.19,3.8e-6,3210,0,[],[],[],'Properties of Advanced Semiconductor Materials: GaN, AIN, InN, BN, SiC, SiGe');

%% Metal material (electrode)
input_material_property_v00('Au',78e9,0.43,14.2e-6,19300,0,[],[],[],[]);
input_material_property_v00('Steel',210e9,0.3,13.7e-6,7800,0,[],[],[],[]);
input_material_property_v00('Aluminium',69e9,0.32,23.1e-6,2700,0,[],[],[],[]);

%% Piezo material
input_material_property_v00('PZT_PIC255',62e9,0.34,4e-6,7800,0,1.55e-8,-180e-12,-180e-12,'PIC255');
input_material_property_v00('PVDF',2.5e9,0.34,140e-6,1750,0,8.854e-11,3e-12,3e-12,[]);

% Thermal research
input_material_property_v00('PZT_PIC255_m1',62e9,0.34,4.12e-6,7800,0,1.55e-8,-180e-12,-180e-12,'PIC255');
input_material_property_v00('PZT_PIC255_p0',62e9,0.34,5.12e-6,7800,0,1.55e-8,-180e-12,-180e-12,'PIC255');
input_material_property_v00('PZT_PIC255_p1',62e9,0.34,6.12e-6,7800,0,1.55e-8,-180e-12,-180e-12,'PIC255');
input_material_property_v00('PZT_PIC255_p2',62e9,0.34,7.12e-6,7800,0,1.55e-8,-180e-12,-180e-12,'PIC255');
input_material_property_v00('PZT_PIC255_p3',62e9,0.34,8.12e-6,7800,0,1.55e-8,-180e-12,-180e-12,'PIC255');

input_material_property_v00('PZT_PIC255_Dummy_m1',62e9,0.34,4.12e-6,7800,0,[],[],[],'PIC255');
input_material_property_v00('PZT_PIC255_Dummy_p0',62e9,0.34,5.12e-6,7800,0,[],[],[],'PIC255');
input_material_property_v00('PZT_PIC255_Dummy_p1',62e9,0.34,6.12e-6,7800,0,[],[],[],'PIC255');
input_material_property_v00('PZT_PIC255_Dummy_p2',62e9,0.34,7.12e-6,7800,0,[],[],[],'PIC255');

%% Piezo material (dummy)
input_material_property_v00('PZT_PIC255_Dummy',62e9,0.34,4e-6,7800,0,[],[],[],[]);
input_material_property_v00('PVDF_Dummy',2.5e9,0.34,140e-6,1750,0,[],[],[],[]);

%% Dummy material (light and soft)
input_material_property_v00('Dummy_Non',1,0.19,5.12e-6,1,0,[],[],[],[]);