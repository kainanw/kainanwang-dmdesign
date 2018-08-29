function [G_matrix,nb_electrode]=geometry_matrix_generation_v00(Geo_Parameter)


% Type_index: The type of the pattern layout: 
% 'Keystone': Keystone layout
% 'Hexanogal': Hexanogal layout
% 'Ring_keystone': Ring keystone layout

if(Geo_Parameter.Ring==1)
    
    Ring_dim=Geo_Parameter.Ring_dim;
    Ring_radius_outer=Ring_dim(1,2);
    Ring_radius_inner=Ring_dim(1,1);
    Ring_section_nb=Ring_dim(2,1);
    
end


if (strcmp(Geo_Parameter.layout,'Keystone'))
    
    % Diameter of the wafer [m]
    D=Geo_Parameter.D;

    % Diameter of the circumcircle of the electrode region [m]
    D_cc_region=Geo_Parameter.D_cc_region;

    if (D_cc_region>D)
        disp('The area of the electrode region is larger than the wafer');
        return;
    end

    % Number of the radial order [/]
    Nb_R=Geo_Parameter.Nb_R;

    % Number of electrodes in each crown [/]
    Nb_e_R=Geo_Parameter.Nb_e_R;
    
    nb_electrode=sum(Nb_e_R);

    if (length(Nb_e_R)~=Nb_R)
        disp('The dimension of the crown-electrode array dismatch, please check it again');
        return;
    elseif (Nb_e_R(1)~=1)
        disp('The first entry of the crown-electrode array is not ONE, please check it again');
        return;
    end

    % Angle of the section of the electrode in each crown [/]
    Angle_section_e_R=360./Nb_e_R;

    % Angle of the electrode shift in each crown [/]
    Angle_shift_e_R=Geo_Parameter.Angle_shift_e_R;

    if (length(Angle_shift_e_R)~=Nb_R)
        disp('The dimension of the angle-shift array dismatch, please check it again');
        return;
    elseif (Angle_shift_e_R(1)~=0)
        disp('The first entry of the angle-shift array is not ZERO, please check it again');
        return;
    end

    % Width of the gap [m]
    W_gap=Geo_Parameter.W_gap;

    % Radius of each crown [m]
    R_e_R=Geo_Parameter.R_e_R;

    if (length(R_e_R)~=Nb_R)
        disp('The dimension of the crown-radius array dismatch, please check it again');
        return;
    elseif (R_e_R~=sort(R_e_R))
        disp('The sequence the crown-radius array is not ordered ascendingly, please check it again');
        return;
    elseif (R_e_R(end)~=D_cc_region/2)
        disp('The last entry of the crown-radius array is incorrect, please check it again');
        return;
    end

    % Define the structure e_origin_temp(ii) and e_temp(ii,jj)
    % ii: the id of the radial order
    % jj: the id of the electrode on each crown

    % Central point coordinate of the micro-hexagonal

    % Central electrode
    e_temp(1,1).radius=R_e_R(1)-0.5*W_gap;

    % Crown electrodes
    for ii=2:Nb_R
        e_origin_temp(ii).R(1)=R_e_R(ii-1)+0.5*W_gap;
        e_origin_temp(ii).R(2)=R_e_R(ii)-0.5*W_gap;
        e_origin_temp(ii).R(3)=R_e_R(ii)-0.5*W_gap;
        e_origin_temp(ii).R(4)=R_e_R(ii-1)+0.5*W_gap;
        e_origin_temp(ii).gap_angle(1)=asin(0.5*W_gap/e_origin_temp(ii).R(1))*180/pi;
        e_origin_temp(ii).gap_angle(2)=asin(0.5*W_gap/e_origin_temp(ii).R(2))*180/pi;
        e_origin_temp(ii).gap_angle(3)=Angle_section_e_R(ii)-e_origin_temp(ii).gap_angle(2);
        e_origin_temp(ii).gap_angle(4)=Angle_section_e_R(ii)-e_origin_temp(ii).gap_angle(1);
    end


    for ii=2:Nb_R
        for jj=1:Nb_e_R(ii)
            e_temp(ii,jj).R(1)=e_origin_temp(ii).R(1);
            e_temp(ii,jj).R(2)=e_origin_temp(ii).R(2);
            e_temp(ii,jj).R(3)=e_origin_temp(ii).R(3);
            e_temp(ii,jj).R(4)=e_origin_temp(ii).R(4);
            e_temp(ii,jj).gap_angle(1)=e_origin_temp(ii).gap_angle(1)+(jj-1)*Angle_section_e_R(ii)+Angle_shift_e_R(ii);
            e_temp(ii,jj).gap_angle(2)=e_origin_temp(ii).gap_angle(2)+(jj-1)*Angle_section_e_R(ii)+Angle_shift_e_R(ii);
            e_temp(ii,jj).gap_angle(3)=e_origin_temp(ii).gap_angle(3)+(jj-1)*Angle_section_e_R(ii)+Angle_shift_e_R(ii);
            e_temp(ii,jj).gap_angle(4)=e_origin_temp(ii).gap_angle(4)+(jj-1)*Angle_section_e_R(ii)+Angle_shift_e_R(ii);

            e_temp(ii,jj).point_x(1)=e_temp(ii,jj).R(1)*cos(e_temp(ii,jj).gap_angle(1)*pi/180);
            e_temp(ii,jj).point_x(2)=e_temp(ii,jj).R(2)*cos(e_temp(ii,jj).gap_angle(2)*pi/180);
            e_temp(ii,jj).point_x(3)=e_temp(ii,jj).R(3)*cos(e_temp(ii,jj).gap_angle(3)*pi/180);
            e_temp(ii,jj).point_x(4)=e_temp(ii,jj).R(4)*cos(e_temp(ii,jj).gap_angle(4)*pi/180);

            e_temp(ii,jj).point_y(1)=e_temp(ii,jj).R(1)*sin(e_temp(ii,jj).gap_angle(1)*pi/180);
            e_temp(ii,jj).point_y(2)=e_temp(ii,jj).R(2)*sin(e_temp(ii,jj).gap_angle(2)*pi/180);
            e_temp(ii,jj).point_y(3)=e_temp(ii,jj).R(3)*sin(e_temp(ii,jj).gap_angle(3)*pi/180);
            e_temp(ii,jj).point_y(4)=e_temp(ii,jj).R(4)*sin(e_temp(ii,jj).gap_angle(4)*pi/180);
        end
    end


    % Change the type of the structure
    % Define the structure electrode(ii)
    % ii: the id of the electrode

    electrode(1).id=1;
    electrode(1).R=e_temp(1,1).radius;
    electrode(1).center_x=0;
    electrode(1).center_y=0;

    id=2;
    for ii=2:Nb_R
        for jj=1:Nb_e_R(ii)
            electrode(id).id=id;
            electrode(id).R=e_temp(ii,jj).R;
            electrode(id).gap_angle=e_temp(ii,jj).gap_angle;
            electrode(id).point_x=e_temp(ii,jj).point_x;
            electrode(id).point_y=e_temp(ii,jj).point_y;
            id=id+1;
        end
    end

    clear id e_temp e_origin_temp

    %% Define Geometry matrix for eletrodes

    % Central electrode
    G_matrix_central_electrode=[1,electrode(1).R,0,0,electrode(1).R,1,length(electrode)+1,0,0,electrode(1).R,0,0;
            1,0,-electrode(1).R,electrode(1).R,0,1,length(electrode)+1,0,0,electrode(1).R,0,0;
            1,-electrode(1).R,0,0,-electrode(1).R,1,length(electrode)+1,0,0,electrode(1).R,0,0;
            1,0,electrode(1).R,-electrode(1).R,0,1,length(electrode)+1,0,0,electrode(1).R,0,0];
    G_matrix_central_electrode=G_matrix_central_electrode';

    % Crown electrode
    G_matrix_crown_electrode=[];
    for ii=2:length(electrode)
        G_matrix=[2,electrode(ii).point_x(1),electrode(ii).point_x(2),...
            electrode(ii).point_y(1),electrode(ii).point_y(2),ii,length(electrode)+1,0,0,0,0,0;
            1,electrode(ii).point_x(2),electrode(ii).point_x(3),...
            electrode(ii).point_y(2),electrode(ii).point_y(3),ii,length(electrode)+1,0,0,electrode(ii).R(2),0,0;
            2,electrode(ii).point_x(3),electrode(ii).point_x(4),...
            electrode(ii).point_y(3),electrode(ii).point_y(4),ii,length(electrode)+1,0,0,0,0,0;
            1,electrode(ii).point_x(1),electrode(ii).point_x(4),...
            electrode(ii).point_y(1),electrode(ii).point_y(4),length(electrode)+1,ii,0,0,electrode(ii).R(1),0,0];
        G_matrix=G_matrix';
        G_matrix_crown_electrode=[G_matrix_crown_electrode,G_matrix];
    end

    % Define Geometry matrix for wafer
    G_matrix_wafer=[1,D/2,0,0,D/2,length(electrode)+1,0,0,0,D/2,0,0;
            1,0,-D/2,D/2,0,length(electrode)+1,0,0,0,D/2,0,0;
            1,-D/2,0,0,-D/2,length(electrode)+1,0,0,0,D/2,0,0;
            1,0,D/2,-D/2,0,length(electrode)+1,0,0,0,D/2,0,0];
    G_matrix_wafer=G_matrix_wafer';

    G_matrix=[G_matrix_central_electrode,G_matrix_crown_electrode,G_matrix_wafer];

elseif (strcmp(Geo_Parameter.layout,'Hexanogal'))
    
    % Diameter of the wafer [m]
    D=Geo_Parameter.D;

    % Diameter of the circumcircle of the hexagonal region [m]
    D_Hex_cc_region=Geo_Parameter.D_Hex_cc_region;

    if (D_Hex_cc_region>0.95*D)
        disp('The area of the hexagonal region is large, please check it again');
        return;
    end

    % Number of the radial order [/]
    Nb_R=Geo_Parameter.Nb_R;

    % Width of the gap [m]
    W_gap=Geo_Parameter.W_gap;

    % Diameter of the inscribed circle of the micro-hexagonal [m]
    D_Hex_ic_micro=(D_Hex_cc_region/(2*Nb_R-2)-W_gap)/2

    if (W_gap>D_Hex_ic_micro/3)
        disp('The width of the gap is large, please check it again!');
        return;
    end
    

    D_Hex_cc_micro=D_Hex_ic_micro*2/sqrt(3);

    % Define the structure mi_hex_temp(ii,jj)
    % ii: the id of the radial order
    % jj: the id of the electrode on each crown

    % Central point coordinate of the micro-hexagonal

    % Central electrode
    mi_hex_temp(1,1).center_x=0;
    mi_hex_temp(1,1).center_y=0;

    % Corner electrodes
    for ii=2:Nb_R
        for jj=1:(ii-1):6*ii
            if ii==2
                mi_hex_temp(ii,jj).center_x=(ii-1)*(2*D_Hex_ic_micro+W_gap)...
                    *cos((jj-1)*60*pi/180);
                mi_hex_temp(ii,jj).center_y=(ii-1)*(2*D_Hex_ic_micro+W_gap)...
                    *sin((jj-1)*60*pi/180);
            else
                mi_hex_temp(ii,jj).center_x=(ii-1)*(2*D_Hex_ic_micro+W_gap)...
                    *cos(floor(jj/(ii-1))*60*pi/180);
                mi_hex_temp(ii,jj).center_y=(ii-1)*(2*D_Hex_ic_micro+W_gap)...
                    *sin(floor(jj/(ii-1))*60*pi/180);
            end
        end
    end

    % Middle electrode between corner electrodes  
    for ii=3:Nb_R
        for jj=1:1:6 % section index
            for mm=(ii-1)*(jj-1)+2:(ii-1)*(jj-1)+ii-1
                mm_start=(ii-1)*(jj-1)+1;
                mm_end=(ii-1)*(jj-1)+ii;
                if mm_end>6*ii
                    mm_end=1;
                end
                mi_hex_temp(ii,mm).center_x=mi_hex_temp(ii,mm_start).center_x...
                    +(mm-(ii-1)*(jj-1)-1)/(ii-1)...
                    *(mi_hex_temp(ii,mm_end).center_x-mi_hex_temp(ii,mm_start).center_x);
                mi_hex_temp(ii,mm).center_y=mi_hex_temp(ii,mm_start).center_y...
                    +(mm-(ii-1)*(jj-1)-1)/(ii-1)...
                    *(mi_hex_temp(ii,mm_end).center_y-mi_hex_temp(ii,mm_start).center_y);
            end
        end
    end

    % Change the type of the structure
    % Define the structure mi_hex(ii)
    % ii: the id of the electrode

    mi_hex(1).id=1;
    mi_hex(1).center_x=mi_hex_temp(1,1).center_x;
    mi_hex(1).center_y=mi_hex_temp(1,1).center_y;

    id_micro=2;
    for ii=1:Nb_R
        for jj=1:6*(ii-1)
            mi_hex(id_micro).id=id_micro;
            mi_hex(id_micro).center_x=mi_hex_temp(ii,jj).center_x;
            mi_hex(id_micro).center_y=mi_hex_temp(ii,jj).center_y;
            id_micro=id_micro+1;
        end
    end

    clear id_micro mi_hex_temp

    % Corner point coordinate of the micro-hexagonal
    for ii=1:length(mi_hex)
        for jj=1:6
            mi_hex(ii).point_x(jj)=mi_hex(ii).center_x+D_Hex_cc_micro*cos((60*jj-30)*pi/180);
            mi_hex(ii).point_y(jj)=mi_hex(ii).center_y+D_Hex_cc_micro*sin((60*jj-30)*pi/180);
        end
    end

    % Define Geometry matrix for eletrodes
    G_matrix_electrode=[];
    for ii=1:length(mi_hex)
        G_matrix=[2,mi_hex(ii).point_x(1),mi_hex(ii).point_x(2),...
            mi_hex(ii).point_y(1),mi_hex(ii).point_y(2),ii,length(mi_hex)+1,0,0,0,0,0;
            2,mi_hex(ii).point_x(2),mi_hex(ii).point_x(3),...
            mi_hex(ii).point_y(2),mi_hex(ii).point_y(3),ii,length(mi_hex)+1,0,0,0,0,0;
            2,mi_hex(ii).point_x(3),mi_hex(ii).point_x(4),...
            mi_hex(ii).point_y(3),mi_hex(ii).point_y(4),ii,length(mi_hex)+1,0,0,0,0,0;
            2,mi_hex(ii).point_x(4),mi_hex(ii).point_x(5),...
            mi_hex(ii).point_y(4),mi_hex(ii).point_y(5),ii,length(mi_hex)+1,0,0,0,0,0;
            2,mi_hex(ii).point_x(5),mi_hex(ii).point_x(6),...
            mi_hex(ii).point_y(5),mi_hex(ii).point_y(6),ii,length(mi_hex)+1,0,0,0,0,0;
            2,mi_hex(ii).point_x(6),mi_hex(ii).point_x(1),...
            mi_hex(ii).point_y(6),mi_hex(ii).point_y(1),ii,length(mi_hex)+1,0,0,0,0,0];
        G_matrix=G_matrix';
        nb_electrode=length(mi_hex);
        
        G_matrix_electrode=[G_matrix_electrode,G_matrix];
    end
    
    %%
    if(Geo_Parameter.Ring==1)
        
        % Angle of the section of the electrode in each crown [/]
        Angle_section_ring=360./Ring_section_nb;

        % Radius of each crown [m]
        % Ring radius
        R_e_Ring=[Ring_radius_inner,Ring_radius_outer];
        
        % Central point coordinate of the micro-hexagonal

        % Central electrode
        e_temp(1,1).radius=R_e_Ring(1)-0.5*W_gap;

        % Crown electrodes
        for ii=2
            e_origin_temp(ii).R(1)=R_e_Ring(ii-1)+0.5*W_gap;
            e_origin_temp(ii).R(2)=R_e_Ring(ii)-0.5*W_gap;
            e_origin_temp(ii).R(3)=R_e_Ring(ii)-0.5*W_gap;
            e_origin_temp(ii).R(4)=R_e_Ring(ii-1)+0.5*W_gap;
            e_origin_temp(ii).gap_angle(1)=asin(0.5*W_gap/e_origin_temp(ii).R(1))*180/pi;
            e_origin_temp(ii).gap_angle(2)=asin(0.5*W_gap/e_origin_temp(ii).R(2))*180/pi;
            e_origin_temp(ii).gap_angle(3)=Angle_section_ring-e_origin_temp(ii).gap_angle(2);
            e_origin_temp(ii).gap_angle(4)=Angle_section_ring-e_origin_temp(ii).gap_angle(1);
        end


        for ii=2
            for jj=1:Ring_section_nb
                e_temp(ii,jj).R(1)=e_origin_temp(ii).R(1);
                e_temp(ii,jj).R(2)=e_origin_temp(ii).R(2);
                e_temp(ii,jj).R(3)=e_origin_temp(ii).R(3);
                e_temp(ii,jj).R(4)=e_origin_temp(ii).R(4);
                e_temp(ii,jj).gap_angle(1)=e_origin_temp(ii).gap_angle(1)+(jj-1)*Angle_section_ring;
                e_temp(ii,jj).gap_angle(2)=e_origin_temp(ii).gap_angle(2)+(jj-1)*Angle_section_ring;
                e_temp(ii,jj).gap_angle(3)=e_origin_temp(ii).gap_angle(3)+(jj-1)*Angle_section_ring;
                e_temp(ii,jj).gap_angle(4)=e_origin_temp(ii).gap_angle(4)+(jj-1)*Angle_section_ring;

                e_temp(ii,jj).point_x(1)=e_temp(ii,jj).R(1)*cos(e_temp(ii,jj).gap_angle(1)*pi/180);
                e_temp(ii,jj).point_x(2)=e_temp(ii,jj).R(2)*cos(e_temp(ii,jj).gap_angle(2)*pi/180);
                e_temp(ii,jj).point_x(3)=e_temp(ii,jj).R(3)*cos(e_temp(ii,jj).gap_angle(3)*pi/180);
                e_temp(ii,jj).point_x(4)=e_temp(ii,jj).R(4)*cos(e_temp(ii,jj).gap_angle(4)*pi/180);

                e_temp(ii,jj).point_y(1)=e_temp(ii,jj).R(1)*sin(e_temp(ii,jj).gap_angle(1)*pi/180);
                e_temp(ii,jj).point_y(2)=e_temp(ii,jj).R(2)*sin(e_temp(ii,jj).gap_angle(2)*pi/180);
                e_temp(ii,jj).point_y(3)=e_temp(ii,jj).R(3)*sin(e_temp(ii,jj).gap_angle(3)*pi/180);
                e_temp(ii,jj).point_y(4)=e_temp(ii,jj).R(4)*sin(e_temp(ii,jj).gap_angle(4)*pi/180);
            end
        end


        % Change the type of the structure
        % Define the structure electrode(ii)
        % ii: the id of the electrode

        ring_electrode(1).id=1;
        ring_electrode(1).R=e_temp(1,1).radius;
        ring_electrode(1).center_x=0;
        ring_electrode(1).center_y=0;

        id=2;
        for ii=2
            for jj=1:Ring_section_nb
                ring_electrode(id).id=id;
                ring_electrode(id).R=e_temp(ii,jj).R;
                ring_electrode(id).gap_angle=e_temp(ii,jj).gap_angle;
                ring_electrode(id).point_x=e_temp(ii,jj).point_x;
                ring_electrode(id).point_y=e_temp(ii,jj).point_y;
                id=id+1;
            end
        end

        clear id e_temp e_origin_temp

        %% Define Geometry matrix for eletrodes

        % Crown electrode
        G_matrix_crown_electrode_ring=[];
        for ii=2:Ring_section_nb+1
%             nb_electrode+1:nb_electrode+Ring_section_nb
            G_matrix=[2,ring_electrode(ii).point_x(1),ring_electrode(ii).point_x(2),...
                ring_electrode(ii).point_y(1),ring_electrode(ii).point_y(2),nb_electrode+ii-1,nb_electrode+Ring_section_nb+1,0,0,0,0,0;
                1,ring_electrode(ii).point_x(2),ring_electrode(ii).point_x(3),...
                ring_electrode(ii).point_y(2),ring_electrode(ii).point_y(3),nb_electrode+ii-1,nb_electrode+Ring_section_nb+1,0,0,ring_electrode(ii).R(2),0,0;
                2,ring_electrode(ii).point_x(3),ring_electrode(ii).point_x(4),...
                ring_electrode(ii).point_y(3),ring_electrode(ii).point_y(4),nb_electrode+ii-1,nb_electrode+Ring_section_nb+1,0,0,0,0,0;
                1,ring_electrode(ii).point_x(1),ring_electrode(ii).point_x(4),...
                ring_electrode(ii).point_y(1),ring_electrode(ii).point_y(4),nb_electrode+Ring_section_nb+1,nb_electrode+ii-1,0,0,ring_electrode(ii).R(1),0,0];
            G_matrix=G_matrix';
            G_matrix_crown_electrode_ring=[G_matrix_crown_electrode_ring,G_matrix];
        end
        
        
            % Define Geometry matrix for wafer
        G_matrix_wafer=[1,D/2,0,0,D/2,nb_electrode+Ring_section_nb+1,0,0,0,D/2,0,0;
                1,0,-D/2,D/2,0,nb_electrode+Ring_section_nb+1,0,0,0,D/2,0,0;
                1,-D/2,0,0,-D/2,nb_electrode+Ring_section_nb+1,0,0,0,D/2,0,0;
                1,0,D/2,-D/2,0,nb_electrode+Ring_section_nb+1,0,0,0,D/2,0,0];
        G_matrix_wafer=G_matrix_wafer';

        G_matrix=[G_matrix_electrode,G_matrix_crown_electrode_ring,G_matrix_wafer];
        
        nb_electrode=nb_electrode+Ring_section_nb;
        
        
    else
        
        % Define Geometry matrix for wafer
        G_matrix_wafer=[1,D/2,0,0,D/2,length(mi_hex)+1,0,0,0,D/2,0,0;
                1,0,-D/2,D/2,0,length(mi_hex)+1,0,0,0,D/2,0,0;
                1,-D/2,0,0,-D/2,length(mi_hex)+1,0,0,0,D/2,0,0;
                1,0,D/2,-D/2,0,length(mi_hex)+1,0,0,0,D/2,0,0];
        G_matrix_wafer=G_matrix_wafer';

        G_matrix=[G_matrix_electrode,G_matrix_wafer];
        
    end
   
elseif (strcmp(Geo_Parameter.layout,'Ring_keystone'))
    
    Geo_Parameter.Nb_R=Geo_Parameter.Nb_R+1;
    Geo_Parameter.Nb_e_R=[1,Geo_Parameter.Nb_e_R];
    Geo_Parameter.Angle_shift_e_R=[0,Geo_Parameter.Angle_shift_e_R];
    
    % Diameter of the wafer [m]
    D=Geo_Parameter.D;

    % Diameter of the circumcircle of the electrode region [m]
    D_cc_region=Geo_Parameter.D_cc_region;

    if (D_cc_region>D)
        disp('The area of the electrode region is larger than the wafer');
        return;
    end

    % Number of the radial order [/]
    Nb_R=Geo_Parameter.Nb_R;

    % Number of electrodes in each crown [/]
    Nb_e_R=Geo_Parameter.Nb_e_R;
    

    if (length(Nb_e_R)~=Nb_R)
        disp('The dimension of the crown-electrode array dismatch, please check it again');
        return;
    elseif (Nb_e_R(1)~=1)
        disp('The first entry of the crown-electrode array is not ONE, please check it again');
        return;
    end

    % Angle of the section of the electrode in each crown [/]
    Angle_section_e_R=360./Nb_e_R;

    % Angle of the electrode shift in each crown [/]
    Angle_shift_e_R=Geo_Parameter.Angle_shift_e_R;

    if (length(Angle_shift_e_R)~=Nb_R)
        disp('The dimension of the angle-shift array dismatch, please check it again');
        return;
    elseif (Angle_shift_e_R(1)~=0)
        disp('The first entry of the angle-shift array is not ZERO, please check it again');
        return;
    end

    % Width of the gap [m]
    W_gap=Geo_Parameter.W_gap;

    % Radius of each crown [m]
    R_e_R=Geo_Parameter.R_e_R;

    if (length(R_e_R)~=Nb_R)
        disp('The dimension of the crown-radius array dismatch, please check it again');
        return;
    elseif (R_e_R~=sort(R_e_R))
        disp('The sequence the crown-radius array is not ordered ascendingly, please check it again');
        return;
    elseif (R_e_R(end)~=D_cc_region/2)
        disp('The last entry of the crown-radius array is incorrect, please check it again');
        return;
    end

    % Define the structure e_origin_temp(ii) and e_temp(ii,jj)
    % ii: the id of the radial order
    % jj: the id of the electrode on each crown

    % Central point coordinate of the micro-hexagonal

    % Central electrode
    e_temp(1,1).radius=R_e_R(1)-0.5*W_gap;

    % Crown electrodes
    for ii=2:Nb_R
        e_origin_temp(ii).R(1)=R_e_R(ii-1)+0.5*W_gap;
        e_origin_temp(ii).R(2)=R_e_R(ii)-0.5*W_gap;
        e_origin_temp(ii).R(3)=R_e_R(ii)-0.5*W_gap;
        e_origin_temp(ii).R(4)=R_e_R(ii-1)+0.5*W_gap;
        e_origin_temp(ii).gap_angle(1)=asin(0.5*W_gap/e_origin_temp(ii).R(1))*180/pi;
        e_origin_temp(ii).gap_angle(2)=asin(0.5*W_gap/e_origin_temp(ii).R(2))*180/pi;
        e_origin_temp(ii).gap_angle(3)=Angle_section_e_R(ii)-e_origin_temp(ii).gap_angle(2);
        e_origin_temp(ii).gap_angle(4)=Angle_section_e_R(ii)-e_origin_temp(ii).gap_angle(1);
    end


    for ii=2:Nb_R
        for jj=1:Nb_e_R(ii)
            e_temp(ii,jj).R(1)=e_origin_temp(ii).R(1);
            e_temp(ii,jj).R(2)=e_origin_temp(ii).R(2);
            e_temp(ii,jj).R(3)=e_origin_temp(ii).R(3);
            e_temp(ii,jj).R(4)=e_origin_temp(ii).R(4);
            e_temp(ii,jj).gap_angle(1)=e_origin_temp(ii).gap_angle(1)+(jj-1)*Angle_section_e_R(ii)+Angle_shift_e_R(ii);
            e_temp(ii,jj).gap_angle(2)=e_origin_temp(ii).gap_angle(2)+(jj-1)*Angle_section_e_R(ii)+Angle_shift_e_R(ii);
            e_temp(ii,jj).gap_angle(3)=e_origin_temp(ii).gap_angle(3)+(jj-1)*Angle_section_e_R(ii)+Angle_shift_e_R(ii);
            e_temp(ii,jj).gap_angle(4)=e_origin_temp(ii).gap_angle(4)+(jj-1)*Angle_section_e_R(ii)+Angle_shift_e_R(ii);

            e_temp(ii,jj).point_x(1)=e_temp(ii,jj).R(1)*cos(e_temp(ii,jj).gap_angle(1)*pi/180);
            e_temp(ii,jj).point_x(2)=e_temp(ii,jj).R(2)*cos(e_temp(ii,jj).gap_angle(2)*pi/180);
            e_temp(ii,jj).point_x(3)=e_temp(ii,jj).R(3)*cos(e_temp(ii,jj).gap_angle(3)*pi/180);
            e_temp(ii,jj).point_x(4)=e_temp(ii,jj).R(4)*cos(e_temp(ii,jj).gap_angle(4)*pi/180);

            e_temp(ii,jj).point_y(1)=e_temp(ii,jj).R(1)*sin(e_temp(ii,jj).gap_angle(1)*pi/180);
            e_temp(ii,jj).point_y(2)=e_temp(ii,jj).R(2)*sin(e_temp(ii,jj).gap_angle(2)*pi/180);
            e_temp(ii,jj).point_y(3)=e_temp(ii,jj).R(3)*sin(e_temp(ii,jj).gap_angle(3)*pi/180);
            e_temp(ii,jj).point_y(4)=e_temp(ii,jj).R(4)*sin(e_temp(ii,jj).gap_angle(4)*pi/180);
        end
    end


    % Change the type of the structure
    % Define the structure electrode(ii)
    % ii: the id of the electrode

    electrode(1).id=1;
    electrode(1).R=e_temp(1,1).radius;
    electrode(1).center_x=0;
    electrode(1).center_y=0;

    id=2;
    for ii=2:Nb_R
        for jj=1:Nb_e_R(ii)
            electrode(id).id=id;
            electrode(id).R=e_temp(ii,jj).R;
            electrode(id).gap_angle=e_temp(ii,jj).gap_angle;
            electrode(id).point_x=e_temp(ii,jj).point_x;
            electrode(id).point_y=e_temp(ii,jj).point_y;
            id=id+1;
        end
    end

    clear id e_temp e_origin_temp

    %% Define Geometry matrix for eletrodes

    % Central electrode
    G_matrix_central_electrode=[1,electrode(1).R,0,0,electrode(1).R,1,length(electrode)+1,0,0,electrode(1).R,0,0;
            1,0,-electrode(1).R,electrode(1).R,0,1,length(electrode)+1,0,0,electrode(1).R,0,0;
            1,-electrode(1).R,0,0,-electrode(1).R,1,length(electrode)+1,0,0,electrode(1).R,0,0;
            1,0,electrode(1).R,-electrode(1).R,0,1,length(electrode)+1,0,0,electrode(1).R,0,0];
    G_matrix_central_electrode=G_matrix_central_electrode';

    % Crown electrode
    G_matrix_crown_electrode=[];
    for ii=2:length(electrode)
        G_matrix=[2,electrode(ii).point_x(1),electrode(ii).point_x(2),...
            electrode(ii).point_y(1),electrode(ii).point_y(2),ii,length(electrode)+1,0,0,0,0,0;
            1,electrode(ii).point_x(2),electrode(ii).point_x(3),...
            electrode(ii).point_y(2),electrode(ii).point_y(3),ii,length(electrode)+1,0,0,electrode(ii).R(2),0,0;
            2,electrode(ii).point_x(3),electrode(ii).point_x(4),...
            electrode(ii).point_y(3),electrode(ii).point_y(4),ii,length(electrode)+1,0,0,0,0,0;
            1,electrode(ii).point_x(1),electrode(ii).point_x(4),...
            electrode(ii).point_y(1),electrode(ii).point_y(4),length(electrode)+1,ii,0,0,electrode(ii).R(1),0,0];
        G_matrix=G_matrix';
        G_matrix_crown_electrode=[G_matrix_crown_electrode,G_matrix];
    end
    
    % Define Geometry matrix for wafer
    G_matrix_wafer=[1,D/2,0,0,D/2,length(electrode)+1,0,0,0,D/2,0,0;
            1,0,-D/2,D/2,0,length(electrode)+1,0,0,0,D/2,0,0;
            1,-D/2,0,0,-D/2,length(electrode)+1,0,0,0,D/2,0,0;
            1,0,D/2,-D/2,0,length(electrode)+1,0,0,0,D/2,0,0];
    G_matrix_wafer=G_matrix_wafer';

    G_matrix=[G_matrix_crown_electrode,G_matrix_wafer];

end

% if (Geo_Parameter.Ring==1)
%     
%     Ring_radius=Geo_Parameter.Ring_radius;
%     G_matrix_electrode=G_matrix(:,1:end-4);
%     
%     G_matrix_ring=[1,Ring_radius,0,0,Ring_radius,nb_electrode+1,0,0,0,Ring_radius,0,0;
%             1,0,-Ring_radius,Ring_radius,0,nb_electrode+1,0,0,0,Ring_radius,0,0;
%             1,-Ring_radius,0,0,-Ring_radius,nb_electrode+1,0,0,0,Ring_radius,0,0;
%             1,0,Ring_radius,-Ring_radius,0,nb_electrode+1,0,0,0,Ring_radius,0,0];
%     G_matrix_ring=G_matrix_ring';
% 
%     G_matrix_outer=[1,D/2,0,0,D/2,nb_electrode+2,0,0,0,D/2,0,0;
%             1,0,-D/2,D/2,0,nb_electrode+2,0,0,0,D/2,0,0;
%             1,-D/2,0,0,-D/2,nb_electrode+2,0,0,0,D/2,0,0;
%             1,0,D/2,-D/2,0,nb_electrode+2,0,0,0,D/2,0,0];
%     G_matrix_outer=G_matrix_outer';
%     
%     G_matrix=[G_matrix_electrode,G_matrix_ring,G_matrix_outer];
%     
% end
    
    
