function C=mul_A_so_2(eqn, opts,opA,B,opB)%#ok<INUSL>

% function C=mul_A_so_2(eqn, opts,opA,B,opB)
%
% The second order system
%
%   M x"(t) + D x'(t) + K x(t) = B u(t)
%                         y(t) = C x(t)
%       
% is transformed to the first order system
%
%   E z'(t) = A z(t) + G u(t)
%      y(t) = L z(t)
%  
% where
%
%      | D  M|
%   E= | M  0|
%   
%      |-K  0|
%   A= | 0  M|
%   
%      | B |
%   G= | 0 |
%   
%   L= [C  0]
%   
%         | x(t)  |
%   z(t)= | x'(t) | .
%   
% Matrices M, D, K are assumed to be quadratic, symmetric and positive definit.
%
%
% This function returns C = A*B, where matrix A given by structure eqn and input matrix B could be transposed. 
% Matrix A is assumed to be quadratic and has a size of 2* size(K).
%
%   Inputs:
%
%   eqn     structure containing  data for matrix A (fields 'K_' and 'M_')
%   opts    structure containing parameters for the algorithm
%   opA     character specifying the shape of A
%           opA = 'N' performs A*opB(B)
%           opA = 'T' performs A'*opB(B)
%   B       m-x-p matrix
%   opB     character specifying the shape of B
%           opB = 'N' performs opA(A)*B
%           opB = 'T' performs opA(A)*B'
%
%   Output:
%
%       |-K  0|
%   C = | 0  M| *opB(B)
%
% This function does not use other so3 functions. 
%
% ATTENTION: opA is not used since matrix A is symmetric

%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others 
%               2009-2020
%


%% check input parameters
if (not(ischar(opA)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opA or opB is not a char');
end

opA = upper(opA); opB = upper(opB);
if(not((opA=='N' || opA=='T')))
    error('MESS:error_arguments','opA is not ''N'' or ''T''');
end

if(not((opB=='N' || opB=='T')))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end

if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if(not(isfield(eqn,'K_')) || not(isnumeric(eqn.K_)) || not(isfield(eqn,'M_')) ...
        || not(isnumeric(eqn.M_)))
    error('MESS:error_arguments',...
        'A consists of K and M, field eqn.K_ or eqn.M_ is not defined');
end


[rowK, colK] = size(eqn.K_);

colA = 2*colK;

%% perform multiplication
switch opB

    % implement operation A*B = C
    case 'N' 
        if(colA ~= size(B,1))
              error('MESS:error_arguments','number of columns of A differs with number of rows of B');
        end
        
        C = [-eqn.K_*B(1:rowK,:) ;...
              eqn.M_*B(rowK+1:end,:)];
            
    % implement operation A*B' = C
    case 'T'
        if(colA ~= size(B,2))
            error('MESS:error_arguments','number of columns of A differs with number of columns of B');
        end
        
        C = [-eqn.K_*B(:,1:colK)';...
              eqn.M_*B(:,colK+1:end)'];
end

end
