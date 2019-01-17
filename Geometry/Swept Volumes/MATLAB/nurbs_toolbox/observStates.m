function [Hb,Hm] = observStates(size_mb,cutx,cuty,hx,bc_up,bc_down,Coeff_forma,p,k)

    % SETTING THE NUMBER OF NODES AND INTERVALS

    ne = round((cutx(2)-cutx(1))/hx);   % Number of Intervals

    nnx = ne * k + p + 1 - k;           % Number of Control Point on the X Direction
                                        % for the IGA Mesh

    % SETTING THE NUMBER OF MODALS IN THE Y DIRECTION

    M = 64;

    % CREATION OF THE MESH CONTAINING THE POINTS IN THE Y DIRECTION
    % USED TO EVALUATE THE MODAL BASIS

    evalNodesY = linspace(cuty(1),cuty(2),M);

    y_elements = length(evalNodesY);

    % INITIALIZATION OF THE SOLUTION MATRICES

    sol   = zeros(nnx,y_elements);
    sol_x = zeros(nnx,y_elements);
    sol_y = zeros(nnx,y_elements);

    import Core.BasisHandler

    obj_newModalBasis = BasisHandler();

    obj_newModalBasis.dimModalBasis         = size_mb;
    obj_newModalBasis.evalNodesY            = evalNodesY;
    obj_newModalBasis.labelUpBoundCond      = bc_up{1};
    obj_newModalBasis.labelDownBoundCond    = bc_down{1};
    obj_newModalBasis.coeffForm             = Coeff_forma;

    [coeffModalBase,~] = newModalBasis(obj_newModalBasis);

    %---------------------------------------------------------------------%
    % Note:
    % The coefficients 'coeffModalBase' and 'coeffModalBaseDer' are
    % vectors conrresponding respectively to the coefficients of the
    % modal bases and derivative of the modal basis on the points
    % assigned in the Y direction.
    %---------------------------------------------------------------------%

    % SETTING OF THE ISOGEOMETRIC MESH IN THE X DIRECTION

    %---------------------------------------------------------------------%
    % NOTE:
    % Here we are considering a linear parametrization of the problem
    % and the continuity requirements achieved are given by C^(p-k).
    %---------------------------------------------------------------------% 

    knot      = augknt([cutx(1) cutx(2)],p+1);
    h         = hx;
    internKnot = linspace(0,1,ne+1)';
    internKnot = internKnot(2:end-1);
    ins       = sort(reshape(internKnot*ones(1,k),1,[]));%k*(nel-1)));
    ins       = (cutx(2)-cutx(1))*ins + cutx(1);

    [~,knot] = bspkntins(p,cutx(1):(cutx(2)-cutx(1))*1/p:cutx(2),knot,ins);

    %---------------------------------------------------------------------%
    % Note:
    % The vector 'knot' contains the nodes used in the mesh.
    %---------------------------------------------------------------------%
    
    xx3 = linspace(cutx(1),cutx(2),nnx)';

    Hb = [];
    
    for m = 1:size_mb

        % COMPUTE ISOGEOMETRIC BASIS
        
        nu = numel(xx3);

        N = zeros(p+1,1);                                

        for col=1:nu                                    

            s = findspan(nnx-1, p, xx3(col), knot);
            N = basisfun(s,xx3(col),p,knot);

            tmp1 = s - p + 1;

            for i=0:p                             
                Hb_TEMP(col,tmp1 + i) = N(i + 1);
            end       
        end
        
        Hb_temp{m} = Hb_TEMP;
        
        % CREATION OF THE BASIS MATRIX

        %Hb_temp{m} = basisMatrix(p,u((m-1)*nnx+1:m*nnx)',knot,xx3);
        Hb = blkdiag(Hb,Hb_temp{m});

    end;
    
    Hm = zeros(nnx * M, nnx * size_mb);

    for h = 1:nnx
        for k = 1:M
            for imb = 1:size_mb
                
                % CREATION OF THE MODAL BASIS MATRIX
                
                Hm(h + (k - 1) * nnx , h + (imb - 1) * nnx) = coeffModalBase(k,imb);

            end
        end
    end

    
%     for ii = 1:M
%         solMat(ii,:) = finalSol((ii-1)*nnx + 1 : ii * nnx);
%     end
%     
%     disp(size(Hm));
%     spy(Hm);
%     disp(norm(solMat - sol',inf));

end