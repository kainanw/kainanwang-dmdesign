function []=preparation_samcef_v00(filename,Geo_Parameter,G_matrix,Mesh_Parameter,Pro_Parameter,Load_Parameter)

h_pre=fopen([filename,'.dat'],'wt');

fprintf(h_pre,['.DEL.*\n']);
fprintf(h_pre,['.UNIT SI\n']);
fprintf(h_pre,['MODE PREC 1e-50\n']);

fprintf(h_pre,['INPUT "',filename,'_geo"\n']);
fprintf(h_pre,['INPUT "',filename,'_mesh"\n']);
fprintf(h_pre,['INPUT "',filename,'_pro"\n']);
fprintf(h_pre,['INPUT "',filename,'_sen"\n']);
fprintf(h_pre,['INPUT "',filename,'_load"\n']);

if (strcmp(Load_Parameter.Load,'IF'))
    fprintf(h_pre,['.SAM\n']);
    fprintf(h_pre,['        MF 1\n']);
    fprintf(h_pre,['        NALG 4\n']);
elseif (strcmp(Load_Parameter.Load,'Thermal'))
    fprintf(h_pre,['.SAM\n']);
    fprintf(h_pre,['        MF 1\n']);
    fprintf(h_pre,['        NALG 4\n']);
    fprintf(h_pre,['        NOP6 0\n']);
elseif (strcmp(Load_Parameter.Load,'Modal_decoupled'))
    fprintf(h_pre,['.SAM\n']);
    fprintf(h_pre,['        MF 1\n']);
    fprintf(h_pre,['        NALG 4\n']);
    fprintf(h_pre,['        NVAL %i\n'],Load_Parameter.Modal_order);
elseif (strcmp(Load_Parameter.Load,'Modal_coupled'))
    fprintf(h_pre,['.SAM\n']);
    fprintf(h_pre,['        MF 1\n']);
    fprintf(h_pre,['        NALG 4\n']);
    fprintf(h_pre,['        NVAL %i\n'],Load_Parameter.Modal_order);
elseif (strcmp(Load_Parameter.Load,'Modal_superelement'))
    fprintf(h_pre,['.SAM\n']);
    fprintf(h_pre,['        MF 1\n']);
    fprintf(h_pre,['        NALG 4\n']);
%     fprintf(h_pre,['        NVAL 10\n']);
end
    
fprintf(h_pre,['.FIN 1\n']);
fprintf(h_pre,['RETURN']);
fclose(h_pre);