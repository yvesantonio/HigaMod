function [Ax,Ay,bx,by] = MyBChandler(Ax,Ay,bx,by,BC,mx,my,fedof,wyhat,modalbx,Jin_out,t)
%       [Ay,by]=DirichletInflow(Ay,by,zeros(my,1),my,fedof); %uy=0 all'inflow
%       [Ay,by]=DirichletOutflow(Ay,by,zeros(my,1),my,fedof); %uy=0 all'outflow
%         
%         %BC per u in inflow
%         if(Data.essential_inflow)
%             [Ax,bx]=DirichletInflow(Ax,bx,Data.inflow_dirichlet,mx,fedof);
%         else
%             bx=NeumannInflow(weights_yhat,bx,Data.Pinf(t),ModalBasisX,fedof,Jin_out);
%         end
%         bx=NeumannOutflow(weights_yhat,bx,Data.Pout(t),ModalBasisX,fedof,Jin_out);
switch BC.INF_x
    case 'dir'
        [Ax,bx] = DirichletInflow(Ax,bx,BC.dataINFx,mx,fedof);
    case 'neu'
        bx = NeumannInflow(wyhat,bx,BC.dataINFx(t),modalbx,fedof,Jin_out');
end

switch BC.INF_y
    case 'dir'
        [Ay,by]=DirichletInflow(Ay,by,BC.dataINFy,my,fedof);
    case 'neu'
        if(BC.dataINFy~=0)
            error(['non homogeneous Neumann conditions for uy not available yet, try to add it\n',...
            'to add the homogeneous one type 0 as the data']);
        end
end

switch BC.OUT_x
    case 'dir'
        error('Dirichlet conditions for ux in outflow not available yet, try to add it');
    case 'neu'
        bx=NeumannOutflow(wyhat,bx,BC.dataOUTx(t),modalbx,fedof,Jin_out');
end


switch BC.OUT_y
    case 'dir'
        [Ay,by]=DirichletOutflow(Ay,by,BC.dataOUTy,my,fedof);
    case 'neu'
        if(BC.dataOUTy~=0)
            error(['non homogeneous Neumann conditions for uy not available yet, try to add it\n',...
            'to add the homogeneous one type 0 as the data']);
        end
end
end