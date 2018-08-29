function active_index=sensor_generation_samcef_v00(sen_file_name,Sensor_Parameter,Geo_Parameter,G_matrix,Nb_electrode,Nb_piezo_index)

h_sen=fopen([sen_file_name,'.dat'],'wt');

% Select the nodes within the pupil to avoid messive nodes in superelement
% generation

fprintf(h_sen,['.SEL\n']);
fprintf(h_sen,['   GROUP NOM "nodes_pupil" NOEUDS\n']);
fprintf(h_sen,['   BOITE CYLINDRIQUE 0 0 RAYON %g\n'],Sensor_Parameter.pupil/2*1.05);

if (strcmp(Sensor_Parameter.type,'Displacement_def'))
    
    active_index=[];
    
    fprintf(h_sen,['RETURN']);
    
elseif (strcmp(Sensor_Parameter.type,'Local_slope'))
    
    Lenslet_width=Sensor_Parameter.pupil/Sensor_Parameter.Nb_lenslet;
    Center_point_coordinate=zeros(Sensor_Parameter.Nb_lenslet^2,4);
    
    % Center point definition
    index=1;
    for ii=1:Sensor_Parameter.Nb_lenslet
        for jj=1:Sensor_Parameter.Nb_lenslet
            x_loc=-Sensor_Parameter.pupil/2+Lenslet_width/2+(ii-1)*Lenslet_width;
            y_loc=-Sensor_Parameter.pupil/2+Lenslet_width/2+(jj-1)*Lenslet_width;
            Center_point_coordinate(index,1)=ii;
            Center_point_coordinate(index,2)=jj;
            Center_point_coordinate(index,3)=x_loc;
            Center_point_coordinate(index,4)=y_loc;
            index=index+1;
            clear x_loc y_loc
        end
    end
    clear ii jj index
    
    % Corner point definition
    Corner_coordinate=zeros(Sensor_Parameter.Nb_lenslet^2,4);
    for ii=1:size(Center_point_coordinate,1)
        x_left=Center_point_coordinate(ii,3)-Lenslet_width/2;
        x_right=Center_point_coordinate(ii,3)+Lenslet_width/2;
        y_bottom=Center_point_coordinate(ii,4)-Lenslet_width/2;
        y_upper=Center_point_coordinate(ii,4)+Lenslet_width/2;
        Corner_coordinate(ii,1)=x_left;
        Corner_coordinate(ii,2)=x_right;
        Corner_coordinate(ii,3)=y_bottom;
        Corner_coordinate(ii,4)=y_upper;
    end
    
    
    figure;
    pdegplot(G_matrix,'EdgeLabels','off','SubdomainLabels','off');
    xlim([-Geo_Parameter.D/2,Geo_Parameter.D/2]);
    ylim([-Geo_Parameter.D/2,Geo_Parameter.D/2]);
    xlabel('X [m]');
    ylabel('Y [m]');
    axis square;
    
    hold on;
    active_index=[];
    for ii=1:size(Center_point_coordinate,1)
        radius_center_point=sqrt(Center_point_coordinate(ii,3)^2+Center_point_coordinate(ii,4)^2);
        if(radius_center_point<=Sensor_Parameter.pupil/2)
            plot(Center_point_coordinate(ii,3),Center_point_coordinate(ii,4),'o','MarkerSize',2,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0]);
            active_index=[active_index,ii];
        else
            plot(Center_point_coordinate(ii,3),Center_point_coordinate(ii,4),'o','MarkerSize',1,'MarkerEdgeColor',[0.8,0.8,0.8],'MarkerFaceColor',[0.8,0.8,0.8]);
        end
    end
    
    for ii=1:size(Center_point_coordinate,1)
        point_1=[Center_point_coordinate(ii,3)-Lenslet_width/2,Center_point_coordinate(ii,4)-Lenslet_width/2];
        point_2=[Center_point_coordinate(ii,3)+Lenslet_width/2,Center_point_coordinate(ii,4)-Lenslet_width/2];
        point_3=[Center_point_coordinate(ii,3)+Lenslet_width/2,Center_point_coordinate(ii,4)+Lenslet_width/2];
        point_4=[Center_point_coordinate(ii,3)-Lenslet_width/2,Center_point_coordinate(ii,4)+Lenslet_width/2];
        plot([point_1(1),point_2(1),point_3(1),point_4(1),point_1(1)],[point_1(2),point_2(2),point_3(2),point_4(2),point_1(2)],'LineWidth',0.5,'Color',[0,0,0]);
    end
    
    draw_circle(Sensor_Parameter.pupil/2);
    
    %% Select the area for S-H sensor lenslet
    fprintf(h_sen,['.SEL\n']);
    for ii=1:length(active_index)
        fprintf(h_sen,['   GROUP NOM "nodes_SH_lenslet_%i" NOEUDS\n'],active_index(ii));
        fprintf(h_sen,['   BOITE XI %g YI %g ZI 0.0 XS %g YS %g ZS 0.0\n'],Corner_coordinate(active_index(ii),1),Corner_coordinate(active_index(ii),3),Corner_coordinate(active_index(ii),2),Corner_coordinate(active_index(ii),4)); 
    end
    
    fprintf(h_sen,['RETURN']);
end
    

fclose(h_sen);