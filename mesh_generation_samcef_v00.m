function [Nb_piezo_index]=mesh_generation_samcef_v00...
    (mesh_file_name,Geo_Parameter,Mesh_parameter,Load_Parameter,...
    Pre_mesh_matrix_1,Nb_structure_1,...
    Pre_mesh_matrix_2,Nb_structure_2)

% Piezo node start index
Nb_piezo_start_index=1e9; % Check the number of the nodes and make sure it is enough big!

switch nargin
    
    case 6
        Pre_mesh_matrix=Pre_mesh_matrix_1;
        Nb_list.line_nb_list=Nb_structure_1.line_nb_list;
        Nb_list.contour_nb_list=Nb_structure_1.contour_nb_list;
        Nb_list.domain_nb_list=Nb_structure_1.domain_nb_list;
        line_segmentation_marker=length(Nb_structure_1.line_nb_list);
        contour_segmentation_marker=length(Nb_structure_1.contour_nb_list);
        domain_segmentation_marker=length(Nb_structure_1.domain_nb_list);
        edge_index=Nb_structure_1.edge_index;
        Surf_coeff=Geo_Parameter.Support_surface;
    case 8
        Pre_mesh_matrix=[Pre_mesh_matrix_1,Pre_mesh_matrix_2];
        Nb_list.line_nb_list=[Nb_structure_1.line_nb_list;...
            Nb_structure_2.line_nb_list];
        Nb_list.contour_nb_list=[Nb_structure_1.contour_nb_list;...
            Nb_structure_2.contour_nb_list];
        Nb_list.domain_nb_list=[Nb_structure_1.domain_nb_list;...
            Nb_structure_2.domain_nb_list];
        line_segmentation_marker=length(Nb_structure_1.line_nb_list);
        contour_segmentation_marker=length(Nb_structure_1.contour_nb_list);
        domain_segmentation_marker=length(Nb_structure_1.domain_nb_list);
        edge_index=Nb_structure_1.edge_index;
        Surf_coeff=Geo_Parameter.Support_surface;
end

degree=Mesh_parameter.degree;

h_mesh=fopen([mesh_file_name,'.dat'],'wt');

fprintf(h_mesh,['! MESH FILE\n']);
fprintf(h_mesh,['! =================================================================================================\n']);
fprintf(h_mesh,['! Generation of the mesh\n']);
fprintf(h_mesh,['! =================================================================================================\n']);

% Define the degree of the element

domain_nb_list=Nb_list.domain_nb_list;
line_nb_list=Nb_list.line_nb_list;

for ii=1:length(domain_nb_list)
    fprintf(h_mesh,['.GEN DEGRE %i DOMAINES %i ATTRIBUT %i\n'],degree,domain_nb_list(ii),domain_nb_list(ii));
end

% if(Geo_Parameter.TriSupport==1)
%     
%     fprintf(h_mesh,['.GEN DEGRE %i DOMAINES %i ATTRIBUT %i\n'],degree,1,1);
%     fprintf(h_mesh,['.GEN DEGRE %i DOMAINES %i ATTRIBUT %i\n'],degree,2,2);
%     fprintf(h_mesh,['.GEN DEGRE %i DOMAINES %i ATTRIBUT %i\n'],degree,3,3);
%     
%     fprintf(h_mesh,['.GEN MODIFIE LIGNE %i ELEMENT %i\n'],1,6);
%     fprintf(h_mesh,['.GEN MODIFIE LIGNE %i ELEMENT %i\n'],2,6);
%     fprintf(h_mesh,['.GEN MODIFIE LIGNE %i ELEMENT %i\n'],3,6);
%     fprintf(h_mesh,['.GEN MODIFIE LIGNE %i ELEMENT %i\n'],4,6);
%     fprintf(h_mesh,['.GEN MODIFIE LIGNE %i ELEMENT %i\n'],5,6);
%     fprintf(h_mesh,['.GEN MODIFIE LIGNE %i ELEMENT %i\n'],6,6);
%     fprintf(h_mesh,['.GEN MODIFIE LIGNE %i ELEMENT %i\n'],7,6);
%     fprintf(h_mesh,['.GEN MODIFIE LIGNE %i ELEMENT %i\n'],8,6);
%     fprintf(h_mesh,['.GEN MODIFIE LIGNE %i ELEMENT %i\n'],9,6);
%     fprintf(h_mesh,['.GEN MODIFIE LIGNE %i ELEMENT %i\n'],10,6);
%     fprintf(h_mesh,['.GEN MODIFIE LIGNE %i ELEMENT %i\n'],11,6);
%     fprintf(h_mesh,['.GEN MODIFIE LIGNE %i ELEMENT %i\n'],12,6);
%     
% end



for jj=1:length(line_nb_list)
    fprintf(h_mesh,['.GEN MODIFIE LIGNE %i ELEMENT %i\n'],Pre_mesh_matrix(16,jj),Pre_mesh_matrix(15,jj));
end

fprintf(h_mesh,['.GEN MAILLE %i A %i TRIANGULATION DELAUNAY\n'],domain_nb_list(1),domain_nb_list(end));
% if(Geo_Parameter.TriSupport==1)
%     fprintf(h_mesh,['.GEN MAILLE %i A %i TRIANGULATION DELAUNAY\n'],1,3);
% end

fprintf(h_mesh,['.COL NOEUDS\n']);
fprintf(h_mesh,['      EXECUTE\n']);
fprintf(h_mesh,['.PURGE PURGE NOEUDS\n']);
fprintf(h_mesh,['.REN NOEUD\n']);
fprintf(h_mesh,['      EXECUTE\n']);
fprintf(h_mesh,['.REN MAILLE\n']);
fprintf(h_mesh,['      EXECUTE\n']);

% Selection of all the nodes and cells
fprintf(h_mesh,['.SEL\n']);
fprintf(h_mesh,['   GROUP NOM "nodes_reflector" NOEUDS TOUT\n']);
fprintf(h_mesh,['   GROUP NOM "cells_reflector" MAILLE TOUT\n']);

% Selection of the electrode
fprintf(h_mesh,['.SEL\n']);
for ii=1:length(domain_nb_list)
    if(ii==length(domain_nb_list))
        fprintf(h_mesh,['   GROUP NOM "cells_gap" MAILLE AT %i\n'],domain_nb_list(ii));
    else
        fprintf(h_mesh,['   GROUP NOM "electrode_%i" MAILLE AT %i\n'],ii,domain_nb_list(ii));
    end
end

% Selection of the edge nodes

edgeline_nb_id=[];
for ii=1:length(edge_index)
    edgeline_nb_id=[edgeline_nb_id,' ',int2str(edge_index(ii))];
end

fprintf(h_mesh,['.SEL\n']);
fprintf(h_mesh,['   GROUP NOM "nodes_edge" NOEUDS LIGNES ',edgeline_nb_id,'\n']);

if(Geo_Parameter.TriSupport==1)
    
    TriFeet_angle=Geo_Parameter.TriFeet_angle;
    TriFeet_radius=Geo_Parameter.TriFeet_radius;
    TriFeet_diameter=Geo_Parameter.TriFeet_diameter;
    
    Trifeet_1=[TriFeet_radius*cos((TriFeet_angle)/180*pi),TriFeet_radius*sin((TriFeet_angle)/180*pi)];
    Trifeet_2=[TriFeet_radius*cos((TriFeet_angle+120)/180*pi),TriFeet_radius*sin((TriFeet_angle+120)/180*pi)];
    Trifeet_3=[TriFeet_radius*cos((TriFeet_angle+240)/180*pi),TriFeet_radius*sin((TriFeet_angle+240)/180*pi)];

    fprintf(h_mesh,['.SEL\n']);
    fprintf(h_mesh,['   GROUP NOM "nodes_trifeet" NOEUDS POINTS 1 2 3\n']);    
end

% Project the nodes
% fprintf(h_mesh,['.MOD NOEUD\n']);
% fprintf(h_mesh,['	LIEU 90 COEFFICIENTS %g %g %g %g %g %g %g %g %g %g DIRECTION 0 0 1\n'],...
%     Surf_coeff(1),Surf_coeff(2),Surf_coeff(3),Surf_coeff(4),Surf_coeff(5),...
%     Surf_coeff(6),Surf_coeff(7),Surf_coeff(8),Surf_coeff(9),Surf_coeff(10));
% fprintf(h_mesh,['	DEPLACE GROUP "nodes_reflector"\n']);

Nb_piezo_index=[];
if (~strcmp(Load_Parameter.Load,'Thermal'))
    fprintf(h_mesh,['.NOE\n']);
    for ii=1:length(domain_nb_list)-1
        fprintf(h_mesh,['     I %i X 0 0 %g\n'],Nb_piezo_start_index+ii,-1);
        Nb_piezo_index=[Nb_piezo_index,Nb_piezo_start_index+ii];
    end
    fprintf(h_mesh,['\n']);
    fprintf(h_mesh,['.SEL\n']);
    fprintf(h_mesh,['   GROUP NOM "nodes_electric" NOEUDS\n']);
    fprintf(h_mesh,['   I %i J %i\n'],Nb_piezo_index(1),Nb_piezo_index(end));
end


fprintf(h_mesh,['RETURN']);
fclose(h_mesh);