function [Aridx, bridx]=impose_boundary_Stokes(impBound)

Ax               = impBound.Ax;
Ay               = impBound.Ay;
bx               = impBound.bx;
by               = impBound.by;
Bxy              = impBound.Bxy;
Byx              = impBound.Byx;
Px               = impBound.Px;
Py               = impBound.Py;
Qx               = impBound.Qx;
Qy               = impBound.Qy;
numbModesUx      = impBound.numbModesUx;
numbModesUy      = impBound.numbModesUy;
infTagUx         = impBound.infTagUx;
outTagUx         = impBound.outTagUx;
infDataUx        = impBound.infDataUx;
outDataUx        = impBound.outDataUx;
infTagUy         = impBound.infTagUy;
outTagUy         = impBound.outTagUy;
infDataUy        = impBound.infDataUy;
outDataUy        = impBound.outDataUy;
numbControlPtsUx = impBound.numbControlPtsUx;
numbControlPtsUy = impBound.numbControlPtsUx;

%%%%%%%%%%%%%%%%%
% IMPOSE Ax, bx %
%%%%%%%%%%%%%%%%%

bridx = [];

if (strcmp(infTagUx,'dir') && strcmp(outTagUx,'dir'))
    
    for imb = 1 : numbModesUx     % Per ogni frequenza
        
        buffAx = zeros( numbControlPtsUx-2, 1); % Inizializzo il contributo al termine noto del nodo di Dirichlet lungo quanto il numero dei nodi interni
        
        for jmb = 1 : numbModesUx
            
            Aridx( (imb-1)*(numbControlPtsUx-2) + 1 : imb*(numbControlPtsUx-2) , (jmb-1)*(numbControlPtsUx-2) + 1 : jmb*(numbControlPtsUx-2) ) = ...
                Ax( (imb-1)*numbControlPtsUx + 2     : imb*numbControlPtsUx - 1 , (jmb-1)*numbControlPtsUx + 2     : jmb*numbControlPtsUx - 1 );
            
            buffAx = buffAx + Ax( (imb-1)*numbControlPtsUx + 2 : imb*numbControlPtsUx-1, (jmb-1)*numbControlPtsUx + 1 )*infDataUx( jmb ) ...
                + Ax( (imb-1)*numbControlPtsUx + 2 : imb*numbControlPtsUx-1, jmb*numbControlPtsUx         )*outDataUx( jmb );
            
        end
        
        bridx = [ bridx; bx( (imb-1)*numbControlPtsUx + 2 : imb*numbControlPtsUx - 1 ) - buffAx ];
        
    end
    
elseif ( strcmp(infTagUx,'dir') && strcmp(outTagUx,'neu') )
    for imb = 1 : numbModesUx
        
        buffAx = zeros( numbControlPtsUx-1, 1);
        
        for jmb = 1 : numbModesUx
            Aridx( (imb-1)*(numbControlPtsUx-1) + 1 : imb*(numbControlPtsUx-1) , (jmb-1)*(numbControlPtsUx-1) + 1 : jmb*(numbControlPtsUx-1) ) = ...
                Ax( (imb-1)*numbControlPtsUx + 2     : imb*numbControlPtsUx     , (jmb-1)*numbControlPtsUx + 2     : jmb*numbControlPtsUx     );
            
            buffAx = buffAx + Ax( (imb-1)*numbControlPtsUx + 2 : imb*numbControlPtsUx, (jmb-1)*numbControlPtsUx + 1 )*infDataUx(jmb);
            
        end
        
        bridx= [ bridx; bx( (imb-1)*numbControlPtsUx + 2 : imb*numbControlPtsUx ) - buffAx ];
        
        bridx(end) = bridx( end ) + outDataUx ( imb ); % Dato di Neumann
    end
    
elseif ( strcmp(infTagUx,'neu') && strcmp(outTagUx,'dir') )
    
    for imb = 1 : numbModesUx
        
        buffAx = zeros( numbControlPtsUx-1, 1);
        
        for jmb = 1 : numbModesUx
            Aridx( (imb-1)*(numbControlPtsUx-1) + 1 : imb*(numbControlPtsUx-1) , (jmb-1)*(numbControlPtsUx-1) + 1 : jmb*(numbControlPtsUx-1) ) = ...
                Ax( (imb-1)*numbControlPtsUx + 1     : imb*numbControlPtsUx-1   , (jmb-1)*numbControlPtsUx + 1     : jmb*numbControlPtsUx - 1 );
            
            buffAx = buffAx + Ax( (imb-1)*numbControlPtsUx + 1 : imb*numbControlPtsUx - 1 , jmb*numbControlPtsUx )*outDataUx( jmb );
        end
        
        bridx = [ bridx ; bx( (imb-1)*numbControlPtsUx + 1 : imb*numbControlPtsUx - 1 ) - buffAx ];
        
        bridx(end-numbControlPtsUx+2)= bridx(end-numbControlPtsUx+2) + infDataUx( imb );
        
    end
elseif ( strcmp(infTagUx,'rob') && (strcmp(outTagUx,'rob') || (strcmp(outTagUx,'neu') ) ) ) %robrob e robneu condividono lo stesso codice.
    
    Aridx=Ax;
    bridx=bx;
    for imb = 1 : numbModesUx
        bridx((imb-1)*numbControlPtsUx+1) = bridx((imb-1)*numbControlPtsUx+1) + infDataUx ( imb );
        bridx(imb*numbControlPtsUx) = bridx(imb*numbControlPtsUx) + outDataUx ( imb );
    end
elseif ( strcmp(infTagUx,'neu') && (strcmp(outTagUx,'neu') ) )
    Aridx=Ax;
    bridx=bx;
    for imb = 1 : numbModesUx
        bridx((imb-1)*numbControlPtsUx+1) = bridx((imb-1)*numbControlPtsUx+1) + infDataUx ( imb );
        bridx(imb*numbControlPtsUx) = bridx(imb*numbControlPtsUx) + outDataUx ( imb );
    end
elseif ( strcmp(infTagUx,'dir') && strcmp(outTagUx,'rob') )
    
    for imb = 1 : numbModesUx
        
        buffAx = zeros( numbControlPtsUx-1, 1);
        
        for jmb = 1 : numbModesUx
            Aridx( (imb-1)*(numbControlPtsUx-1) + 1 : imb*(numbControlPtsUx-1) , (jmb-1)*(numbControlPtsUx-1) + 1 : jmb*(numbControlPtsUx-1) ) = ...
                Ax( (imb-1)*numbControlPtsUx + 2     : imb*numbControlPtsUx     , (jmb-1)*numbControlPtsUx + 2     : jmb*numbControlPtsUx     );
            buffAx = buffAx + Ax( (imb-1)*numbControlPtsUx + 2 : imb*numbControlPtsUx, (jmb-1)*numbControlPtsUx + 1 )*infDataUx(jmb);
            
        end
        
        bridx= [ bridx; bx( (imb-1)*numbControlPtsUx + 2 : imb*numbControlPtsUx ) - buffAx ];
        bridx(end) = bridx( end ) + outDataUx ( imb ); % Dato di Robin Robin
    end
end

%%%%%%%%%%%%%%%%%
% IMPOSE Ay, by %
%%%%%%%%%%%%%%%%%

bridy = [];

if (strcmp(infTagUy,'dir') && strcmp(outTagUy,'dir'))
    
    for imb = 1 : numbModesUy     % Per ogni frequenza
        
        buffAy = zeros( numbControlPtsUy-2, 1); % Inizializzo il contributo al termine noto del nodo di Dirichlet lungo quanto il numero dei nodi interni
        
        for jmb = 1 : numbModesUy
            
            Aridy( (imb-1)*(numbControlPtsUy-2) + 1 : imb*(numbControlPtsUy-2) , (jmb-1)*(numbControlPtsUy-2) + 1 : jmb*(numbControlPtsUy-2) ) = ...
                Ay( (imb-1)*numbControlPtsUy + 2     : imb*numbControlPtsUy - 1 , (jmb-1)*numbControlPtsUy + 2     : jmb*numbControlPtsUy - 1 );
            
            buffAy = buffAy + Ay( (imb-1)*numbControlPtsUy + 2 : imb*numbControlPtsUy-1, (jmb-1)*numbControlPtsUy + 1 )*infDataUy( jmb ) ...
                + Ay( (imb-1)*numbControlPtsUy + 2 : imb*numbControlPtsUy-1, jmb*numbControlPtsUy         )*outDataUy( jmb );
            
        end
        
        bridy = [ bridy; by( (imb-1)*numbControlPtsUy + 2 : imb*numbControlPtsUy - 1 ) - buffAy ];
        
    end
    
elseif ( strcmp(infTagUy,'dir') && strcmp(outTagUy,'neu') )
    for imb = 1 : numbModesUy
        
        buffAy = zeros( numbControlPtsUy-1, 1);
        
        for jmb = 1 : numbModesUy
            Aridy( (imb-1)*(numbControlPtsUy-1) + 1 : imb*(numbControlPtsUy-1) , (jmb-1)*(numbControlPtsUy-1) + 1 : jmb*(numbControlPtsUy-1) ) = ...
                Ay( (imb-1)*numbControlPtsUy + 2     : imb*numbControlPtsUy     , (jmb-1)*numbControlPtsUy + 2     : jmb*numbControlPtsUy     );
            
            buffAy = buffAy + Ay( (imb-1)*numbControlPtsUy + 2 : imb*numbControlPtsUy, (jmb-1)*numbControlPtsUy + 1 )*infDataUy(jmb);
            
        end
        
        bridy= [ bridy; by( (imb-1)*numbControlPtsUy + 2 : imb*numbControlPtsUy ) - buffAy ];
        
        bridy(end) = bridy( end ) + outDataUy ( imb ); % Dato di Neumann
    end
    
elseif ( strcmp(infTagUy,'neu') && strcmp(outTagUy,'dir') )
    
    for imb = 1 : numbModesUy
        
        buffAy = zeros( numbControlPtsUy-1, 1);
        
        for jmb = 1 : numbModesUy
            Aridy( (imb-1)*(numbControlPtsUy-1) + 1 : imb*(numbControlPtsUy-1) , (jmb-1)*(numbControlPtsUy-1) + 1 : jmb*(numbControlPtsUy-1) ) = ...
                Ay( (imb-1)*numbControlPtsUy + 1     : imb*numbControlPtsUy-1   , (jmb-1)*numbControlPtsUy + 1     : jmb*numbControlPtsUy - 1 );
            
            buffAy = buffAy + Ay( (imb-1)*numbControlPtsUy + 1 : imb*numbControlPtsUy - 1 , jmb*numbControlPtsUy )*outDataUy( jmb );
        end
        
        bridy = [ bridy ; by( (imb-1)*numbControlPtsUy + 1 : imb*numbControlPtsUy - 1 ) - buffAy ];
        
        bridy(end-numbControlPtsUy+2)= bridy(end-numbControlPtsUy+2) + infDataUy( imb );
        
    end
elseif ( strcmp(infTagUy,'rob') && (strcmp(outTagUy,'rob') || (strcmp(outTagUy,'neu') ) ) ) %robrob e robneu condividono lo stesso codice.
    
    Aridy=Ay;
    bridy=by;
    for imb = 1 : numbModesUy
        bridy((imb-1)*numbControlPtsUy+1) = bridy((imb-1)*numbControlPtsUy+1) + infDataUy ( imb );
        bridy(imb*numbControlPtsUy) = bridy(imb*numbControlPtsUy) + outDataUy ( imb );
    end
elseif ( strcmp(infTagUy,'neu') && (strcmp(outTagUy,'neu') ) )
    Aridy=Ay;
    bridy=by;
    for imb = 1 : numbModesUy
        bridy((imb-1)*numbControlPtsUy+1) = bridy((imb-1)*numbControlPtsUy+1) + infDataUy ( imb );
        bridy(imb*numbControlPtsUy) = bridy(imb*numbControlPtsUy) + outDataUy ( imb );
    end
elseif ( strcmp(infTagUy,'dir') && strcmp(outTagUy,'rob') )
    
    for imb = 1 : numbModesUy
        
        buffAy = zeros( numbControlPtsUy-1, 1);
        
        for jmb = 1 : numbModesUy
            Aridy( (imb-1)*(numbControlPtsUy-1) + 1 : imb*(numbControlPtsUy-1) , (jmb-1)*(numbControlPtsUy-1) + 1 : jmb*(numbControlPtsUy-1) ) = ...
                Ay( (imb-1)*numbControlPtsUy + 2     : imb*numbControlPtsUy     , (jmb-1)*numbControlPtsUy + 2     : jmb*numbControlPtsUy     );
            buffAy = buffAy + Ay( (imb-1)*numbControlPtsUy + 2 : imb*numbControlPtsUy, (jmb-1)*numbControlPtsUy + 1 )*infDataUy(jmb);
            
        end
        
        bridy= [ bridy; by( (imb-1)*numbControlPtsUy + 2 : imb*numbControlPtsUy ) - buffAy ];
        bridy(end) = bridy( end ) + outDataUy ( imb ); % Dato di Robin Robin
    end
end
