function []=input_material_property_v00(Name_str,YM,PR,TEC,D,TREF,Permittivity,d31,d32,Comment_str)

%% Input new matrial

Material_name=Name_str;

% Material property
Material.YM=YM; % Young's Modulus [Pa]
Material.PR=PR; % Poisson Ratio [/]
Material.TEC=TEC; % Thermal expansion coefficient [/]
Material.D=D; % Mass density [kg/m^3]
Material.TREF=TREF; % Reference temperature [K]

% Piezo
Material.PZEE=Permittivity; % Permittivity
Material.PZUE_1=d31; % d31
Material.PZUE_2=d32; % d32

Material.N=Material_name; % Name

time_ind=clock;% Time record
Material.TIME=[num2str(time_ind(1)),'-',...
    num2str(time_ind(2),2),'-',...
    num2str(time_ind(3),2),'-',...
    num2str(time_ind(4),2),'-',...
    num2str(time_ind(5),2)]; % Documentation time

Comment=Comment_str; % Comments
Material.Comment=Comment; % Comments

str_1=[Material_name,'=Material;'];
evalc(str_1);
% save('MATERIAL_DATABASE',Material,'-append');
save('MATERIAL_DATABASE',Material_name,'-append');
%% Load the material database

% load('MATERIAL_DATABASE');

% switch Force_mode
%     case 0
%         if (exist(Material_name)==1)
%             str_2=['Material_repeat=',Material_name,';'];
%             evalc(str_2);
%             disp('The material exists in the database, please check it again');
%             disp('');
%             disp(['Name:                          ',num2str(Material_repeat.N,5)]);
%             disp(['Young modulus:                 ',num2str(Material_repeat.YM,5)]);
%             disp(['Poisson ratio:                 ',num2str(Material_repeat.PR,5)]);
%             disp(['Density:                       ',num2str(Material_repeat.D,5)]);
%             disp(['Thermal expansion Coefficient: ',num2str(Material_repeat.TEC,5)]);
%         else
%             save('MATERIAL_DATABASE',Material_name,'-append');
%         end
%     case 1
%         Material_name=Material;
%         str_3=[Material_name,'=Material;'];
%         evalc(str_3);
        save('MATERIAL_DATABASE',Material_name,'-append');
