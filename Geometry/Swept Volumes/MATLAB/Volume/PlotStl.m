% load an ascii STL sample file (STLGETFORMAT and STLREADASCII)

az = [0 45 90 135 180];
el = [90 90 90 90 90];

for ii = 1:length(az)
    
    % VIEW PARAMETERS
    
    param.az = az(ii);
    param.el = el(ii);
    param.alpha = 1.0;
    
    % BUILD FIGURE
    
    [vertices,faces,normals,name] = stlRead('COR_ART.stl');
    name = '  ';
    stlPlot2(vertices,faces,name,param);
    
    % EXPORT FIGURE
    
    figName = ['COR_ART_AZ',num2str(az(ii)),'_EL',num2str(el(ii)),'.png'];
    export_fig(figName,'-transparent','-append')
    
end

