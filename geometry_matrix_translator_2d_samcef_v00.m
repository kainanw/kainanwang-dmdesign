function [Pre_mesh_matrix,Nb_structure]=...
    geometry_matrix_translator_2d_samcef_v00...
    (geometry_file_name,...
    Geo_Parameter,...
    G_matrix,...
    line_nb_start,contour_nb_start,domain_nb_start,...
    stage_marker)

% The geometry matrix should be ordered as: 1)electrodes, and 2)wafer
% The subdomain number should be ordered as: 1)electrodes, and 2)wafer

% Geo_Parameter: the structure contains basic information of the geometry

% The numbering of the line/arc, contour, and domain starts from
% line_nb_start, contour_nb_start, domain_nb_start

% stage_marker: the stage of the finite element formulaiton
%               1: the generation of the line/contour/domain of electrodes
%               and wafer(external boundary)
%               2: the generation the line/contour/domain of the other
%               objectives

%% Internal parameter

Sparse_mesh=Geo_Parameter.Sparse_factor;

Sparse_mesh_factor=Sparse_mesh(1);
Sparse_mesh_factor_boundary=Sparse_mesh(2);

% Width of the gap

W_gap=Geo_Parameter.W_gap;

% Segmentation of the geometry matrix

region_index_record(1)=0;
region_index_record(2)=0;

seg_index=0;
seg_nb=zeros(size(G_matrix,2),1);
seg_subdomain=zeros(size(G_matrix,2),1);

for ii=1:size(G_matrix,2)
    
   region_index(1)=G_matrix(6,ii);
   region_index(2)=G_matrix(7,ii);
   
   if ((region_index(1)~=region_index_record(1))||...
           (region_index(2)~=region_index_record(2)))&&...
           ((region_index(1)~=region_index_record(2))||...
           (region_index(2)~=region_index_record(1)))
       seg_index=seg_index+1;
       seg_nb(seg_index)=seg_nb(seg_index)+1;
       
   elseif ((region_index(1)==region_index_record(1))&&...
           (region_index(2)==region_index_record(2)))||...
           ((region_index(1)==region_index_record(2))&&...
           (region_index(2)==region_index_record(1)))
       seg_nb(seg_index)=seg_nb(seg_index)+1;
       
       if min(region_index(1),region_index(2))~=0 
           seg_subdomain(seg_index)=min(region_index(1),region_index(2));
       else % find the external boundary
           seg_subdomain(seg_index)=max(region_index(1),region_index(2));
       end
       
   else
       disp('Error');
       
   end
   
   region_index_record(1)=G_matrix(6,ii);
   region_index_record(2)=G_matrix(7,ii);
end

seg_nb=seg_nb(1:seg_index);
seg_subdomain=seg_subdomain(1:seg_index);

h_geo=fopen([geometry_file_name,'.dat'],'wt');

fprintf(h_geo,['! GEOMETRY FILE\n']);
fprintf(h_geo,['! =================================================================================================\n']);

if (stage_marker==1)
    fprintf(h_geo,['! Generation of the line/contour/domain of electrodes and external boundary\n']);
elseif (stage_marker==2)
    fprintf(h_geo,['! Generation the line/contour/domain of the other objectives\n']);
else
    disp('Please check the stage marker!');
end
fprintf(h_geo,['! =================================================================================================\n']);

% Plan define
plan_z=0;
fprintf(h_geo,['.PLAN I 1 Z %g\n'],plan_z);

if(Geo_Parameter.TriSupport==1)
    
    TriFeet_angle=Geo_Parameter.TriFeet_angle;
    TriFeet_radius=Geo_Parameter.TriFeet_radius;
    TriFeet_diameter=Geo_Parameter.TriFeet_diameter;
    TriFeet_piezo_intersection=Geo_Parameter.TriFeet_piezo_intersection;
    TriFeet_location_domain=Geo_Parameter.TriFeet_location_domain;

    Trifeet_1=[TriFeet_radius*cos((TriFeet_angle)/180*pi),TriFeet_radius*sin((TriFeet_angle)/180*pi)];
    Trifeet_2=[TriFeet_radius*cos((TriFeet_angle+120)/180*pi),TriFeet_radius*sin((TriFeet_angle+120)/180*pi)];
    Trifeet_3=[TriFeet_radius*cos((TriFeet_angle+240)/180*pi),TriFeet_radius*sin((TriFeet_angle+240)/180*pi)];
    
    fprintf(h_geo,['.3POINT I 1 X %g Y %g Z %g\n'],Trifeet_1(1),Trifeet_1(2),plan_z);
    fprintf(h_geo,['.3POINT I 2 X %g Y %g Z %g\n'],Trifeet_2(1),Trifeet_2(2),plan_z);
    fprintf(h_geo,['.3POINT I 3 X %g Y %g Z %g\n'],Trifeet_3(1),Trifeet_3(2),plan_z);
    
    if(TriFeet_piezo_intersection==1)
        fprintf(h_geo,['.CON I 1 POINTS 1\n']);
        fprintf(h_geo,['.CON I 2 POINTS 2\n']);
        fprintf(h_geo,['.CON I 3 POINTS 3\n']);
    elseif(TriFeet_piezo_intersection==0)
        fprintf(h_geo,['.CON I 1 POINTS 1 2 3\n']);
    end
    
    Trifeet_marker=zeros(length(seg_nb),1);
    for ii=1:length(TriFeet_location_domain)
        Trifeet_marker(TriFeet_location_domain(ii))=1;
    end
           
end

edge_index=[];
counter_marker_trifeet=1;

for ii=1:length(seg_nb)
    
    part_nb=(sum(seg_nb(1:ii))-seg_nb(ii)+1):1:sum(seg_nb(1:ii));
    
    for jj=1:length(part_nb)
        
        if (G_matrix(1,part_nb(jj))==1) % arc
            
            if(ii==length(seg_nb))
                fprintf(h_geo,['.3ARC I %i CCENTER %g %g %g CPOINTS %g %g %g %g %g %g\n'],part_nb(jj)+line_nb_start,...
                G_matrix(8,part_nb(jj)),G_matrix(9,part_nb(jj)),plan_z,...
                G_matrix(2,part_nb(jj)),G_matrix(4,part_nb(jj)),plan_z,...
                G_matrix(3,part_nb(jj)),G_matrix(5,part_nb(jj)),plan_z);
                OA=[G_matrix(2,part_nb(jj))-G_matrix(8,part_nb(jj)),...
                    G_matrix(4,part_nb(jj))-G_matrix(9,part_nb(jj))];
                OB=[G_matrix(3,part_nb(jj))-G_matrix(8,part_nb(jj)),...
                    G_matrix(5,part_nb(jj))-G_matrix(9,part_nb(jj))];
                theta=acos((OA*OB')/(norm(OA)*norm(OB)));
                dimension_line(part_nb(jj))=G_matrix(10,part_nb(jj))*theta;
                nb_mesh_seed(part_nb(jj))=ceil(G_matrix(10,part_nb(jj))*theta/W_gap*Sparse_mesh_factor_boundary);
                index_line(part_nb(jj))=part_nb(jj)+line_nb_start;
                edge_index=[edge_index,part_nb(jj)+line_nb_start];
            else
                fprintf(h_geo,['.3ARC I %i CCENTER %g %g %g CPOINTS %g %g %g %g %g %g\n'],part_nb(jj)+line_nb_start,...
                G_matrix(8,part_nb(jj)),G_matrix(9,part_nb(jj)),plan_z,...
                G_matrix(2,part_nb(jj)),G_matrix(4,part_nb(jj)),plan_z,...
                G_matrix(3,part_nb(jj)),G_matrix(5,part_nb(jj)),plan_z);
                OA=[G_matrix(2,part_nb(jj))-G_matrix(8,part_nb(jj)),...
                    G_matrix(4,part_nb(jj))-G_matrix(9,part_nb(jj))];
                OB=[G_matrix(3,part_nb(jj))-G_matrix(8,part_nb(jj)),...
                    G_matrix(5,part_nb(jj))-G_matrix(9,part_nb(jj))];
                theta=acos((OA*OB')/(norm(OA)*norm(OB)));
                dimension_line(part_nb(jj))=G_matrix(10,part_nb(jj))*theta;
                nb_mesh_seed(part_nb(jj))=ceil(G_matrix(10,part_nb(jj))*theta/W_gap*Sparse_mesh_factor);
                index_line(part_nb(jj))=part_nb(jj)+line_nb_start;
            end
                
        elseif (G_matrix(1,part_nb(jj))==2) % line
            
            if(ii==length(seg_nb))
                fprintf(h_geo,['.3DROI I %i CPOINTS %g %g %g %g %g %g\n'],part_nb(jj)+line_nb_start,...
                G_matrix(2,part_nb(jj)),G_matrix(4,part_nb(jj)),plan_z,...
                G_matrix(3,part_nb(jj)),G_matrix(5,part_nb(jj)),plan_z);
                OR=[G_matrix(3,part_nb(jj))-G_matrix(2,part_nb(jj)),...
                    G_matrix(5,part_nb(jj))-G_matrix(4,part_nb(jj))];
                dimension_line(part_nb(jj))=norm(OR);
                nb_mesh_seed(part_nb(jj))=ceil(norm(OR)/W_gap*Sparse_mesh_factor_boundary);
                index_line(part_nb(jj))=part_nb(jj)+line_nb_start;
                edge_index=[edge_index,part_nb(jj)+line_nb_start];
            else
                fprintf(h_geo,['.3DROI I %i CPOINTS %g %g %g %g %g %g\n'],part_nb(jj)+line_nb_start,...
                G_matrix(2,part_nb(jj)),G_matrix(4,part_nb(jj)),plan_z,...
                G_matrix(3,part_nb(jj)),G_matrix(5,part_nb(jj)),plan_z);
                OR=[G_matrix(3,part_nb(jj))-G_matrix(2,part_nb(jj)),...
                    G_matrix(5,part_nb(jj))-G_matrix(4,part_nb(jj))];
                dimension_line(part_nb(jj))=norm(OR);
                nb_mesh_seed(part_nb(jj))=ceil(norm(OR)/W_gap*Sparse_mesh_factor);
                index_line(part_nb(jj))=part_nb(jj)+line_nb_start;
            end
                
        else
            disp('Undefined geometry.');
            return;
            
        end
        
        if (ii==length(seg_nb))
            part_nb=(sum(seg_nb(1:ii))-seg_nb(ii)+1):1:sum(seg_nb(1:ii));
            line_nb_end=part_nb(end)+line_nb_start;
        end
            
    end
    
    part_nb_id=[];
    for tt=1:length(part_nb)
        part_nb_id=[part_nb_id,[num2str(part_nb(tt)+line_nb_start),' ']];
    end
    fprintf(h_geo,['.CON I %i LIGNES ',part_nb_id,'FERME\n'],ii+contour_nb_start);
    contour_nb_list(ii)=ii+contour_nb_start;
    
    
    if (ii==length(seg_nb))
        contour_nb_end=ii+contour_nb_start;
    end
    
    if(Geo_Parameter.TriSupport==0)
        if (ii~=length(seg_nb))
            fprintf(h_geo,['.DOM I %i CONTOURS %i SURFACE 1\n'],ii+domain_nb_start,ii+contour_nb_start);
        elseif (ii==length(seg_nb))
            con_nb_id=[int2str(contour_nb_start+1),' A ',int2str(length(seg_nb)-1+contour_nb_start)];
            fprintf(h_geo,['.DOM I %i CONTOURS %i ',con_nb_id,' SURFACE 1\n'],ii+domain_nb_start,ii+contour_nb_start);
            domain_nb_end=ii+domain_nb_start;
        else
            disp('Error');
            return
        end
    elseif(Geo_Parameter.TriSupport==1)
        if(Trifeet_marker(ii)==0)
            if (ii~=length(seg_nb))
                fprintf(h_geo,['.DOM I %i CONTOURS %i SURFACE 1\n'],ii+domain_nb_start,ii+contour_nb_start);
            elseif (ii==length(seg_nb))
                con_nb_id=[int2str(contour_nb_start+1),' A ',int2str(length(seg_nb)-1+contour_nb_start)];
                fprintf(h_geo,['.DOM I %i CONTOURS %i ',con_nb_id,' SURFACE 1\n'],ii+domain_nb_start,ii+contour_nb_start);
                domain_nb_end=ii+domain_nb_start;
            else
                disp('Error');
                return
            end
        elseif(Trifeet_marker(ii)==1)
            
            if(TriFeet_piezo_intersection==1)
                if (ii~=length(seg_nb))
                    fprintf(h_geo,['.DOM I %i CONTOURS %i %i SURFACE 1\n'],ii+domain_nb_start,ii+contour_nb_start,counter_marker_trifeet);
                elseif (ii==length(seg_nb))
                    con_nb_id=[int2str(contour_nb_start+1),' A ',int2str(length(seg_nb)-1+contour_nb_start)];
                    fprintf(h_geo,['.DOM I %i CONTOURS %i ',con_nb_id,' SURFACE 1\n'],ii+domain_nb_start,ii+contour_nb_start);
                    domain_nb_end=ii+domain_nb_start;
                else
                    disp('Error');
                    return
                end
                counter_marker_trifeet=counter_marker_trifeet+1;
            elseif(TriFeet_piezo_intersection==0)
                if (ii~=length(seg_nb))
                    fprintf(h_geo,['.DOM I %i CONTOURS %i SURFACE 1\n'],ii+domain_nb_start,ii+contour_nb_start);
                elseif (ii==length(seg_nb))
                    con_nb_id=[int2str(contour_nb_start+1),' A ',int2str(length(seg_nb)-1+contour_nb_start)];
                    fprintf(h_geo,['.DOM I %i CONTOURS %i 1 ',con_nb_id,' SURFACE 1\n'],ii+domain_nb_start,ii+contour_nb_start);
                    domain_nb_end=ii+domain_nb_start;
                else
                    disp('Error');
                    return
                end
            end
            
        end
    end
    domain_nb_list(ii)=ii+domain_nb_start;
    
end

% if(Geo_Parameter.TriSupport==1)
%     
%     fprintf(h_geo,['.CON I 1 LIGNES 1 2 3 4 FERME\n']);
%     fprintf(h_geo,['.CON I 2 LIGNES 5 6 7 8 FERME\n']);
%     fprintf(h_geo,['.CON I 3 LIGNES 9 10 11 12 FERME\n']);
%     
%     fprintf(h_geo,['.DOM I 1 CONTOURS 1 SURFACE 1\n']);
%     fprintf(h_geo,['.DOM I 2 CONTOURS 2 SURFACE 1\n']);
%     fprintf(h_geo,['.DOM I 3 CONTOURS 3 SURFACE 1\n']);
% 
% 
% end




Pre_mesh_matrix=[G_matrix;zeros(8,size(G_matrix,2))];
Pre_mesh_matrix(14,:)=dimension_line;
Pre_mesh_matrix(15,:)=nb_mesh_seed;
Pre_mesh_matrix(16,:)=index_line;
Pre_mesh_matrix(20,:)=stage_marker*ones(1,size(G_matrix,2));

line_nb_list=index_line';
contour_nb_list=contour_nb_list';
domain_nb_list=domain_nb_list';

Nb_structure.line_nb_list=line_nb_list;
Nb_structure.line_nb_start=line_nb_start;
Nb_structure.line_nb_end=line_nb_end;

Nb_structure.contour_nb_list=contour_nb_list;
Nb_structure.contour_nb_start=contour_nb_start;
Nb_structure.contour_nb_end=contour_nb_end;

Nb_structure.domain_nb_list=domain_nb_list;
Nb_structure.domain_nb_start=domain_nb_start;
Nb_structure.domain_nb_end=domain_nb_end;

Nb_structure.edge_index=edge_index;

% if(Geo_Parameter.Ring==1)
%     Ring_radius=Geo_Parameter.Ring_radius;
%     Outer_radius=Ring_radius(2);
%     fprintf(h_geo,['.3ARC I %i CCENTRE %g %g %g CPOINTS %g %g %g %g %g %g\n'],line_nb_list(end)+1,0,0,0,Outer_radius,0,0,0,Outer_radius,0);
%     fprintf(h_geo,['.3ARC I %i CCENTRE %g %g %g CPOINTS %g %g %g %g %g %g\n'],line_nb_list(end)+2,0,0,0,0,Outer_radius,0,-Outer_radius,0,0);
%     fprintf(h_geo,['.3ARC I %i CCENTRE %g %g %g CPOINTS %g %g %g %g %g %g\n'],line_nb_list(end)+3,0,0,0,-Outer_radius,0,0,0,-Outer_radius,0);
%     fprintf(h_geo,['.3ARC I %i CCENTRE %g %g %g CPOINTS %g %g %g %g %g %g\n'],line_nb_list(end)+4,0,0,0,0,-Outer_radius,0,Outer_radius,0,0);
%     fprintf(h_geo,['.CON I %i LIGNES ',int2str(line_nb_list(end)+1),' ',int2str(line_nb_list(end)+2),' ',int2str(line_nb_list(end)+3),' ',int2str(line_nb_list(end)+4),' ','FERME\n'],contour_nb_list(end)+1);
%     Inner_radius=Ring_radius(1);
%     fprintf(h_geo,['.3ARC I %i CCENTRE %g %g %g CPOINTS %g %g %g %g %g %g\n'],line_nb_list(end)+5,0,0,0,Inner_radius,0,0,0,Inner_radius,0);
%     fprintf(h_geo,['.3ARC I %i CCENTRE %g %g %g CPOINTS %g %g %g %g %g %g\n'],line_nb_list(end)+6,0,0,0,0,Inner_radius,0,-Inner_radius,0,0);
%     fprintf(h_geo,['.3ARC I %i CCENTRE %g %g %g CPOINTS %g %g %g %g %g %g\n'],line_nb_list(end)+7,0,0,0,-Inner_radius,0,0,0,-Inner_radius,0);
%     fprintf(h_geo,['.3ARC I %i CCENTRE %g %g %g CPOINTS %g %g %g %g %g %g\n'],line_nb_list(end)+8,0,0,0,0,-Inner_radius,0,Inner_radius,0,0);
%     fprintf(h_geo,['.CON I %i LIGNES ',int2str(line_nb_list(end)+1),' ',int2str(line_nb_list(end)+2),' ',int2str(line_nb_list(end)+3),' ',int2str(line_nb_list(end)+4),' ','FERME\n'],contour_nb_list(end)+1);
% end

fprintf(h_geo,['RETURN']);
fclose(h_geo);