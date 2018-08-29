function []=voltage_map(Voltage,Nb_electrode,v_min,v_max,G_matrix)

%% Check the dimension of the voltage


%% Segmentation of the geometry matrix

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

% figure;
hold on;
colormap(jet);
caxis([v_min,v_max]);
colorbar;
box on;
axis square;

for ii=1:length(seg_nb)-1
    
    domain_nb=seg_subdomain(ii);
	x_vector=[];
    y_vector=[];
    
    part_nb=(sum(seg_nb(1:ii))-seg_nb(ii)+1):1:sum(seg_nb(1:ii));
    
    for jj=1:length(part_nb)
        

        if (G_matrix(1,part_nb(jj))==1) % arc
            
            if(ii==length(seg_nb))

            else
                
                radius_est_1=sqrt((G_matrix(4,part_nb(jj))-G_matrix(9,part_nb(jj)))^2+(G_matrix(2,part_nb(jj))-G_matrix(8,part_nb(jj)))^2);
                radius_est_2=sqrt((G_matrix(5,part_nb(jj))-G_matrix(9,part_nb(jj)))^2+(G_matrix(3,part_nb(jj))-G_matrix(8,part_nb(jj)))^2);
                radius_est=0.5*(radius_est_1+radius_est_2);
                
                
                angle_1=atan2(G_matrix(4,part_nb(jj))-G_matrix(9,part_nb(jj)),G_matrix(2,part_nb(jj))-G_matrix(8,part_nb(jj)));
                angle_2=atan2(G_matrix(5,part_nb(jj))-G_matrix(9,part_nb(jj)),G_matrix(3,part_nb(jj))-G_matrix(8,part_nb(jj)));
                
                
                if (angle_1>angle_2)
                    angle_2=2*pi+angle_2;
                end
                
                angle_seg=linspace(angle_1,angle_2,100);
                x_coord_arc=zeros(1,length(angle_seg));
                y_coord_arc=zeros(1,length(angle_seg));
                
                for tt=1:length(angle_seg)
                    angle_temp=angle_seg(tt);
                    x_coord_arc(tt)=cos(angle_temp)*radius_est+G_matrix(8,part_nb(jj));
                    y_coord_arc(tt)=sin(angle_temp)*radius_est+G_matrix(9,part_nb(jj));
                end
                
                x_vector_temp=[G_matrix(2,part_nb(jj)),x_coord_arc,G_matrix(3,part_nb(jj))];
                y_vector_temp=[G_matrix(4,part_nb(jj)),y_coord_arc,G_matrix(5,part_nb(jj))];
                
                x_vector=[x_vector,x_vector_temp];
                y_vector=[y_vector,y_vector_temp];
                    
            end
                
        elseif (G_matrix(1,part_nb(jj))==2) % line
            
            if(ii==length(seg_nb))

            else

                x_vector_temp=[G_matrix(2,part_nb(jj)),G_matrix(3,part_nb(jj))];
                y_vector_temp=[G_matrix(4,part_nb(jj)),G_matrix(5,part_nb(jj))];
                
                x_vector=[x_vector,x_vector_temp];
                y_vector=[y_vector,y_vector_temp];
                
            end
                
        else
            disp('Undefined geometry.');
            return;
            
        end

    end
    
    patch(x_vector,y_vector,Voltage(domain_nb),'EdgeColor','none');
    
end
