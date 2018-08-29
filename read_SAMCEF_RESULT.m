function [coord_nodes,indices_nodes]=read_SAMCEF_RESULT(filename,module_type,type,nc,group)

if (strcmp(type,'coordinate'))
        
    h_node_coord_x_in=fopen([group,'_x.in'],'w');
    fprintf(h_node_coord_x_in,['$$GET_VALUE  "Code 153" "Component 1" "GROUP [',group,']" " " " "\n']);
    fprintf(h_node_coord_x_in,['$$END\n']);
    evalc(['!samcef sr ',filename,' lcp=',module_type,' request=',group,'_x.in answers=',group,'_x.out']);% call to SAMRES
    
    h_node_coord_x_out=fopen([group,'_x.out'],'r');
    fgetl(h_node_coord_x_out);
    fgetl(h_node_coord_x_out);
    fgetl(h_node_coord_x_out);
    A=fgetl(h_node_coord_x_out);
    nb_nodes=sscanf(A,'%i',[1 1]);
    fgetl(h_node_coord_x_out);

    indices_nodes_x=zeros(nb_nodes,1);
    coord_nodes_x=zeros(nb_nodes,1);
    for tt=1:1:nb_nodes
        A=fgetl(h_node_coord_x_out);
        indices_nodes_x(tt)=sscanf(A,'%i',[1 1]);
    end
    for tt=1:1:nb_nodes
        A=fgetl(h_node_coord_x_out);
        coord_nodes_x(tt)=sscanf(A,'%g',[1 1]);
    end
    fclose(h_node_coord_x_out);
    
    h_node_coord_y_in=fopen([group,'_y.in'],'w');
    fprintf(h_node_coord_y_in,['$$GET_VALUE  "Code 153" "Component 2" "GROUP [',group,']" " " " "\n']);
    fprintf(h_node_coord_y_in,['$$END\n']);
    evalc(['!samcef sr ',filename,' lcp=',module_type,' request=',group,'_y.in answers=',group,'_y.out']);% call to SAMRES
    
    h_node_coord_y_out=fopen([group,'_y.out'],'r');
    fgetl(h_node_coord_y_out);
    fgetl(h_node_coord_y_out);
    fgetl(h_node_coord_y_out);
    A=fgetl(h_node_coord_y_out);
    nb_nodes=sscanf(A,'%i',[1 1]);
    fgetl(h_node_coord_y_out);

    indices_nodes_y=zeros(nb_nodes,1);
    coord_nodes_y=zeros(nb_nodes,1);
    for tt=1:1:nb_nodes
        A=fgetl(h_node_coord_y_out);
        indices_nodes_y(tt)=sscanf(A,'%i',[1 1]);
    end
    for tt=1:1:nb_nodes
        A=fgetl(h_node_coord_y_out);
        coord_nodes_y(tt)=sscanf(A,'%g',[1 1]);
    end
    fclose(h_node_coord_y_out);
    
    [existing_number,ia,ib]=intersect(indices_nodes_x,indices_nodes_y);
    
    indices_nodes_x=indices_nodes_x(ia);
    coord_nodes_x=coord_nodes_x(ia);
    
    indices_nodes_y=indices_nodes_y(ib);
    coord_nodes_y=coord_nodes_y(ib);
    
    
    indices_nodes=indices_nodes_x;
    coord_nodes=[coord_nodes_x,coord_nodes_y];
    fclose('all');
    
%     delete([group,'_x.in']);
%     delete([group,'_y.in']);
%     delete([group,'_x.out']);
%     delete([group,'_y.out']);

elseif (strcmp(type,'displacement'))
    
    h_node_coord_z_in=fopen([group,'_z.in'],'w');
    fprintf(h_node_coord_z_in,['$$GET_VALUE  "Code 163" "Component 3" "GROUP [',group,']" " " "4 %i"\n'],nc);
    fprintf(h_node_coord_z_in,['$$END\n']);
    evalc(['!samcef sr ',filename,' lcp=',module_type,' request=',group,'_z.in answers=',group,'_z.out']);% call to SAMRES
    
    h_node_coord_z_out=fopen([group,'_z.out'],'r');
    fgetl(h_node_coord_z_out);
    fgetl(h_node_coord_z_out);
    fgetl(h_node_coord_z_out);
    A=fgetl(h_node_coord_z_out);
    nb_nodes=sscanf(A,'%i',[1 1]);
    fgetl(h_node_coord_z_out);

    indices_nodes=zeros(nb_nodes,1);
    coord_nodes=zeros(nb_nodes,1);
    for tt=1:1:nb_nodes
        A=fgetl(h_node_coord_z_out);
        indices_nodes(tt)=sscanf(A,'%i',[1 1]);
    end
    for tt=1:1:nb_nodes
        A=fgetl(h_node_coord_z_out);
        coord_nodes(tt)=sscanf(A,'%g',[1 1]);
    end
    fclose(h_node_coord_z_out);
    fclose('all');
    
%     delete([group,'_z.in']);
%     delete([group,'_z.out']);
    
elseif (strcmp(type,'modalshape'))
    
    h_node_coord_z_in=fopen([group,'_z.in'],'w');
    fprintf(h_node_coord_z_in,['$$GET_VALUE  "Code 163" "Component 3" "GROUP [',group,']" " " "6 %i"\n'],nc);
    fprintf(h_node_coord_z_in,['$$END\n']);
    evalc(['!samcef sr ',filename,' lcp=',module_type,' request=',group,'_z.in answers=',group,'_z.out']);% call to SAMRES
    
    h_node_coord_z_out=fopen([group,'_z.out'],'r');
    fgetl(h_node_coord_z_out);
    fgetl(h_node_coord_z_out);
    fgetl(h_node_coord_z_out);
    A=fgetl(h_node_coord_z_out);
    nb_nodes=sscanf(A,'%i',[1 1]);
    fgetl(h_node_coord_z_out);

    indices_nodes=zeros(nb_nodes,1);
    coord_nodes=zeros(nb_nodes,1);
    for tt=1:1:nb_nodes
        A=fgetl(h_node_coord_z_out);
        indices_nodes(tt)=sscanf(A,'%i',[1 1]);
    end
    for tt=1:1:nb_nodes
        A=fgetl(h_node_coord_z_out);
        coord_nodes(tt)=sscanf(A,'%g',[1 1]);
    end
    fclose(h_node_coord_z_out);
    fclose('all');
    
%     delete([group,'_z.in']);
%     delete([group,'_z.out']);
    
end

