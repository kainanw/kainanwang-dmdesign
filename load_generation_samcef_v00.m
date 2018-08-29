function []=load_generation_samcef_v00...
    (load_file_name,Geo_parameter,Load_Parameter,nb_electrode,nb_piezo_index)

BC=Load_Parameter.BC;
Load=Load_Parameter.Load;

h_load=fopen([load_file_name,'.dat'],'wt');

if (strcmp(BC,'Isostatic'))
    fprintf(h_load,['.FRAME I 1 TYPE CYLINDRICAL\n']);
    fprintf(h_load,['       V1 1 0 0 V2 0 1 0\n']);
    fprintf(h_load,['.AXL GROUP "nodes_edge" FRAME 1\n']);
    fprintf(h_load,['.CLM\n']);
    fprintf(h_load,['       FIX GROUP "nodes_edge" COMP 2 3\n']);       
elseif (strcmp(BC,'Pin_supported'))
    fprintf(h_load,['.FRAME I 1 TYPE CYLINDRICAL\n']);
    fprintf(h_load,['       V1 1 0 0 V2 0 1 0\n']);
    fprintf(h_load,['.AXL GROUP "nodes_edge" FRAME 1\n']);
    fprintf(h_load,['.CLM\n']);
    fprintf(h_load,['       FIX GROUP "nodes_edge" COMP 1 2 3 4 6\n']);  
elseif (strcmp(BC,'Clamped'))
    fprintf(h_load,['.CLM\n']);
    fprintf(h_load,['       FIX GROUP "nodes_edge" COMP 1 2 3 4 5 6\n']);
elseif (strcmp(BC,'TriFeet_isostatic'))
    fprintf(h_load,['.FRAME I 1 TYPE CYLINDRICAL\n']);
    fprintf(h_load,['       V1 1 0 0 V2 0 1 0\n']);
    fprintf(h_load,['.AXL GROUP "nodes_trifeet" FRAME 1\n']);
%     fprintf(h_load,['.AXL GROUP "nodes_trifeet_2" FRAME 1\n']);
%     fprintf(h_load,['.AXL GROUP "nodes_trifeet_3" FRAME 1\n']);
    fprintf(h_load,['.CLM\n']);
    fprintf(h_load,['       FIX NOE GROUP "nodes_trifeet" COMP 2 3\n']);
%     fprintf(h_load,['       FIX NOE GROUP "nodes_trifeet_2" COMP 2 3\n']);
%     fprintf(h_load,['       FIX NOE GROUP "nodes_trifeet_3" COMP 2 3\n']);
elseif (strcmp(BC,'TriFeet_clamped'))
    fprintf(h_load,['.CLM\n']);
    fprintf(h_load,['       FIX NOE GROUP "nodes_trifeet" COMP 1 2 3 4 5 6\n']);
%     fprintf(h_load,['       FIX NOE GROUP "nodes_trifeet_2" COMP 1 2 3 4 5 6\n']);
%     fprintf(h_load,['       FIX NOE GROUP "nodes_trifeet_3" COMP 1 2 3 4 5 6\n']);
end

if (strcmp(Load,'IF'))
    fprintf(h_load,['! Electrical load\n']);
    fprintf(h_load,['.CLM\n']);
    load_case_index=1;
    for ii=1:nb_electrode
        V_act=zeros(nb_electrode,1);
        V_act(ii)=1;
        fprintf(h_load,['! IF' int2str(ii) '\n']);
        for ind_act=1:nb_electrode
            fprintf(h_load,['       DEP NOE %i COMP 1 VAL %g NC %i\n'],nb_piezo_index(ind_act),V_act(ind_act),load_case_index);
        end
        load_case_index=load_case_index+1;
    end
    
elseif(strcmp(Load,'Modal_decoupled'))
    
    fprintf(h_load,['! Modal analysis (Decoupled)\n']);
    fprintf(h_load,['.CLM\n']);
    for ind_act=1:nb_electrode
        fprintf(h_load,['       DEP NOE %i COMP 1 VAL 0\n'],nb_piezo_index(ind_act));
    end
    
elseif(strcmp(Load,'Modal_coupled'))
    
    Elec_choice=Load_Parameter.Elec_choice;
    
    fprintf(h_load,['! Modal analysis (Coupled)\n']);
    fprintf(h_load,['.CLM\n']);
    for ind_act=1:nb_electrode
        if(Elec_choice(ind_act)==1)
            fprintf(h_load,['       DEP NOE %i COMP 1 VAL 0\n'],nb_piezo_index(ind_act));
        end
    end
    
elseif(strcmp(Load,'Modal_superelement'))
    
    fprintf(h_load,['.RET\n']);
    fprintf(h_load,['   GROUP "nodes_electric" C 1\n']);
    fprintf(h_load,['   GROUP "nodes_pupil" C 3\n']);
    fprintf(h_load,['.CSUP\n']);
    fprintf(h_load,['   RET DYNAMIQUE COMPOSANTS LANCZOS NVAL %i\n'],Load_Parameter.Modal_reconstruction);
    fprintf(h_load,['	GENERE\n']);
    
elseif(strcmp(Load,'Thermal'))
    
    Temperature=Load_Parameter.Temperature;
    fprintf(h_load,['.CLT\n']);
    fprintf(h_load,['   GROUP "nodes_reflector" TFX val %g\n'],Temperature);
    
elseif (strcmp(Load,'Stroke'))
    
    if(Geo_parameter.Ring==1)
        nb_actuator=nb_electrode-Geo_parameter.Ring_dim(2,1);
    end
        
    fprintf(h_load,['! Electrical load\n']);
    fprintf(h_load,['.CLM\n']);
    for ii=1:nb_actuator
        fprintf(h_load,['       DEP NOE %i COMP 1 VAL %g NC %i\n'],nb_piezo_index(ii),Load_Parameter.Stroke_voltage,1);
    end

end

fprintf(h_load,['RETURN']);
fclose(h_load);
    