function fourierCoefficients=project(func,m,bc_up,bc_down,problemCoefficients)
% it projects func on the first m modal basis defined by bc_up,bc_down and
% the coefficients stored in problemCoefficients

% the nodes and correspongind weigths
nbQuadNodes=64;
[~, yq, wyq] = gausslegendre( nbQuadNodes );

%the modal basis must be evaluated on nodes that belongs to [0,1]
modalBasis=new_modal_basis(m, yq, bc_up, bc_down, problemCoefficients);

% we evaluate the input function on the quadrature nodes
funcEvaluated=func(yq);

fourierCoefficients=zeros(m,1);
for jMode=1:m
    fourierCoefficients(jMode,1) = integrate( funcEvaluated'.*modalBasis(:,jMode)' , wyq );
end

% for debugging and visualization porpouses
% we reconstruct the function trough its fourier coefficients
funcReconstructed=zeros(nbQuadNodes,1);
for jMode=1:m
    funcReconstructed=funcReconstructed+modalBasis(:,jMode)*fourierCoefficients(jMode,1);
end

figure
hold on
plot(yq,funcEvaluated)
plot(yq,funcReconstructed,'r');
legend('the function', 'the projection')
title(['Function projected on ', num2str(m), ' modal basis of the type ', bc_up,'-',bc_down])
end