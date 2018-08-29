function []=property_generation_samcef_v00...
    (pro_file_name,Geo_Parameter,Pro_Parameter,Load_Parameter,nb_electrode,nb_piezo_index)

if(Geo_Parameter.Ring==1)
    nb_actuator=nb_electrode-Geo_Parameter.Ring_dim(2,1);
else
    nb_actuator=nb_electrode;
end

Hypothesis=Pro_Parameter.hypothesis;

Material=Pro_Parameter.material;
nb_material=length(Material);

% Default all the substrate
Substrate_material=Pro_Parameter.sequence;
Substrate_thickness=Pro_Parameter.thickness;
nb_substrate=length(Substrate_material);

% (De)Active part sequence
Sequence_act=Pro_Parameter.sequence_act;
Sequence_act_char=[];
for tt=1:length(Sequence_act)
    Sequence_act_char=[Sequence_act_char,int2str(Sequence_act(tt)),' '];
end

Sequence_deact=Pro_Parameter.sequence_deact;
Sequence_deact_char=[];
for tt=1:length(Sequence_deact)
    Sequence_deact_char=[Sequence_deact_char,int2str(Sequence_deact(tt)),' '];
end

if(Geo_Parameter.Ring==1)
    Sequence_shunt=Pro_Parameter.sequence_shunt;
    Sequence_shunt_char=[];
    for tt=1:length(Sequence_shunt)
        Sequence_shunt_char=[Sequence_shunt_char,int2str(Sequence_shunt(tt)),' '];
    end
end

h_pro=fopen([pro_file_name,'.dat'],'wt');

fprintf(h_pro,['! PROPERTY FILE\n']);
fprintf(h_pro,['! =================================================================================================\n']);

% Define hypothesis

fprintf(h_pro,['.HYP ',Hypothesis,'\n']);

% Define material

fprintf(h_pro,['.MAT\n']);
for ii=1:nb_material
    
    Material_name=Material(ii);
    Material_name=char(Material_name);
    
    load('MATERIAL_DATABASE.mat',Material_name);
    
    str_1=['Material_input=',Material_name,';'];
    evalc(str_1);
    
    N=Material_input.N; % String
    YM=Material_input.YM;
    PR=Material_input.PR;
    TEC=Material_input.TEC;
    D=Material_input.D;
    TREF=Material_input.TREF;
    
    % Piezo material check
    if (strcmp(char(Material_input.PZEE),''))
        piezo_flag=0;
    else
        piezo_flag=1;
    end
    
    if (piezo_flag==1)
        PZEE=Material_input.PZEE;
        PZUE_1=Material_input.PZUE_1;
        PZUE_2=Material_input.PZUE_2;
    end
    
    
    if (piezo_flag==1)
        fprintf(h_pro,['	I %i NOM "',Material_name,'" YT %g NT %g A %g M %g TREF %g PZEE %g PZUE %g %g\n'],...
            ii,YM,PR,TEC,D,TREF,PZEE,PZUE_1,PZUE_2);
    else
        fprintf(h_pro,['	I %i NOM "',Material_name,'" YT %g NT %g A %g M %g TREF %g\n'],...
            ii,YM,PR,TEC,D,TREF);
    end
    
    clear Material_name
        
end

% % Dummy material
% fprintf(h_pro,['.MAT\n']);
% fprintf(h_pro,['	I %i NOM "Dummy" YT %g NT %g A %g M %g TREF %g\n'],...
%             nb_material+1,1,0.14,4e-6,1,0);

% Define the ply
fprintf(h_pro,['.PLI\n']);
Material_angle=zeros(size(Substrate_material)); % Temporary setting
for ii=1:nb_substrate
    fprintf(h_pro,['	AN %g T %g MAT %i PLI %i\n'],...
        Material_angle(ii),Substrate_thickness(ii),Substrate_material(ii),ii);
end

% % Dummy layer
% fprintf(h_pro,['.PLI\n']);
% fprintf(h_pro,['	AN %g T %g MAT %i PLI %i\n'],...
%         0,sum(Substrate_thickness),nb_material+1,nb_substrate+1);
    
% Define the laminate

fprintf(h_pro,['.LAM\n']);
% Active part
fprintf(h_pro,['	LAM 1 PLI ',Sequence_act_char,'\n']);
% Deactive part
fprintf(h_pro,['	LAM 2 PLI ',Sequence_deact_char,'\n']);
if (Geo_Parameter.Ring==1)
    % Shunt part
    fprintf(h_pro,['	LAM 3 PLI ',Sequence_shunt_char,'\n']);
end


% Apply the property onto the substrate
Substrate_angle=zeros(nb_actuator,1); % Temporary setting
fprintf(h_pro,['.ETA\n']);
for ii=1:nb_actuator
    fprintf(h_pro,['	GROUP "electrode_%i" LAM 1 ANGLE %g\n'],ii,Substrate_angle(ii));
end

if (Geo_Parameter.Ring==1)
    Substrate_angle_shunt=zeros(nb_electrode-nb_actuator+1,1); % Temporary setting
    for ii=nb_actuator+1:nb_electrode
        fprintf(h_pro,['	GROUP "electrode_%i" LAM 3 ANGLE %g\n'],ii,Substrate_angle_shunt(ii-nb_actuator));
    end
end

fprintf(h_pro,['	GROUP "cells_gap" LAM 2 ANGLE %g\n'],0);

if (~strcmp(Load_Parameter.Load,'Thermal'))
    % Associate the electrical nodes to the electrodes
    fprintf(h_pro,['.AEL\n']);
    for ii=1:nb_electrode
        fprintf(h_pro,['	GROUP "electrode_%i" EN %i\n'],ii,nb_piezo_index(ii));
    end
end

fprintf(h_pro,['RETURN']);
fclose(h_pro);