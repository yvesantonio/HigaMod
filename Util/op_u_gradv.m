function varargout = op_u_gradv (spu, spv, msh, coeff)

  shpu = reshape (spu.shape_functions, spu.ncomp, ...
                  msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, ...
                   [],msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradv, 2);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  jacdet_weights = msh.jacdet .* msh.quad_weights .* coeff;
  
  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      shpu_iel = reshape (shpu(:,:,1:spu.nsh(iel),iel), spu.ncomp, msh.nqn, 1, spu.nsh(iel));
%       shpu_iel = repmat (shpu_iel, [1,1,spv.nsh(iel),1]);
      gradv_iel = reshape (gradv(:,:,:,1:spv.nsh(iel),iel), spv.ncomp*ndir, msh.nqn, spv.nsh(iel), 1);
%       gradv_iel = repmat (gradv_iel, [1,1,1,spu.nsh(iel)]);

      jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);
      
      jacdet_shpu = bsxfun (@times, jacdet_iel, shpu_iel);
      tmp1 = sum (bsxfun (@times, jacdet_shpu, gradv_iel), 1);
      values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = reshape (sum (tmp1, 2), spv.nsh(iel), spu.nsh(iel));

      [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
      rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc;
      cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc;
      ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);

    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_shpu_gradv: singular map in element number %d', iel)
    end
  end

  if (nargout == 1 || nargout == 0)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_shpu_gradv: wrong number of output arguments')
  end

end