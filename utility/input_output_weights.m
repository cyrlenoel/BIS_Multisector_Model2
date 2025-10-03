% input_output_weights
%
% This file contains the input-output weights as well as the labour and
% capital shares
%
% =========================================================================

% IO requirements matrix
io_requirements = io.io_inputs ;

% Construct intermediate input weights
ii_weights = io_requirements(1 : end-2, :) ; 
ii_weights = ii_weights * diag(sum(ii_weights,1).^(-1)) ;
     
% Labour weights in value added
lva_weights = io_requirements(end-1,:) ./ (io_requirements(end-1,:) + io_requirements(end,:)) ;

% Value added vs intermediate input weights
va_weights = sum(io_requirements(end-1:end,:),1)./sum(io_requirements,1) ;
                




