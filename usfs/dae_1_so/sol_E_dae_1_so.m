function X = sol_E_dae_1_so(eqn, opts, opE, B, opB)%#ok<INUSL>
%% function sol_E_dae_1_so solves opE(E)*X = opB(B) resp. performs X=opE(E)\opB(B)
%
% Input:
%   eqn     structure contains data for E (M_,K_)
%
%   opts    struct contains parameters for the algorithm
%
%   opE     character specifies the form of opE(E)
%           opE = 'N' solves E *X = opB(B)
%           opE = 'T' sovles E'*X = opB(B)
%
%   B       p-x-q matrix 
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' solves opE(E)*X = B
%           opB = 'T' solves opE(E)*X = B'
%
% Output
%
%   X       matrix fullfills equation opE(E)*X = opB(B)
%
%   uses no other dae_1_so function
%% check input Paramters

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
if (not(ischar(opE)) || not(ischar(opB)))
    error('MESS:error_arguments', 'opE or opB is not a char');
end

opE = upper(opE); opB = upper(opB);
if(not((opE == 'N' || opE == 'T')))
    error('MESS:error_arguments','opE is not ''N'' or ''T''');
end

if(not((opB == 'N' || opB == 'T')))
    error('MESS:error_arguments','opB is not ''N'' or ''T''');
end
if (not(isnumeric(B))) || (not(ismatrix(B)))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
if (not(isfield(eqn,'M_')) || not(isnumeric(eqn.M_)))
    error('MESS:equation_data',...
        'Empty or Corrupted field M detected in equation structure.')
elseif (not(isfield(eqn,'E_')) || not(isnumeric(eqn.E_)))
    error('MESS:equation_data',...
        'Empty or Corrupted field D detected in equation structure.')
end
if (not(isfield(eqn,'isSym')))
    isSym = 0;
else
    isSym = eqn.isSym;
end
if not(isfield(eqn, 'nd'))    || not(isnumeric(eqn.nd))
    error('MESS:nd',...
    'Missing or Corrupted nd field detected in equation structure.');
end

nd = eqn.nd;
if(opB == 'N')
    rowB = size(B, 1);
else
    rowB = size(B, 2);
end

if(2 * nd ~= rowB)
    error('MESS:error_arguments', 'number of rows of B differs from number of cols of E ( 2 * nd)');
end


%% solve
if isSym
    switch opE
        
        case 'N'
            switch opB
                
                
                case 'N'
                    Db = eqn.E_(1 : nd, 1 : nd) * B(1:nd,:);
                    X1 = eqn.M_(1 : nd, 1 : nd) \ (B(nd + 1 : end, :) - Db);
                    X = [X1; B(1 : nd, :)];
                    
                    
                case 'T'
                    Db = eqn.E_ (1 : nd, 1 : nd) * B(:,1:nd)';
                    X1 = eqn.M_(1 : nd, 1 : nd) \ (B( : , nd + 1 : end)' - Db);
                    X = [X1; B( : , 1 : nd)'];
            end
            
        case 'T'
            switch opB
                
                
                case 'N'
                    X2 = eqn.M_(1 : nd, 1 : nd) \ B(1 : nd, :);
                    X1 = eqn.E_(1 : nd, 1 : nd) * X2;
                    X = [B(nd + 1 : end, :) - X1; X2];
                    
                    
                case 'T'
                    X2 = eqn.M_(1 : nd, 1 : nd) \ B( : , 1 : nd)';
                    X1 = eqn.E_(1 : nd, 1 : nd) * X2;
                    X = [B( : , nd + 1 : end)' - X1(1 : nd, :); X2];
            end
            
    end
else
    switch opE
        
        case 'N'
            switch opB
                
                
                case 'N'
                    Db = eqn.E_(1 : nd, 1 : nd) * B(1 : nd , : );
                    X1 = eqn.M_(1 : nd, 1 : nd) \ (B(nd + 1 : end, :) - Db);
                    X = [X1; B(1 : nd, :)];
                    
                    
                case 'T'
                    Db = eqn.E_ * B(:,1:nd)';
                    X1 = eqn.M_(1 : nd, 1 : nd) \ (B( : , nd + 1 : end)' - Db);
                    X = [X1; B( : , 1 : nd)'];
            end
            
        case 'T'
            switch opB
                
                
                case 'N'
                    X2 = eqn.M_(1 : nd, 1 : nd)' \ B(1 : nd, :);
                    X1 = eqn.E_(1 : nd, 1: nd)' * X2;
                    X = [B(nd + 1 : end, :) - X1(1 : nd, :); X2];
                    
                    
                case 'T'
                    X2 = eqn.M_(1 : nd, 1 : nd)' \ B( : , 1 : nd)';
                    X1 = eqn.E_(1 : nd, 1 : nd)' * X2;
                    X = [B( : , nd + 1 : end)' - X1(1 : nd, :); X2];
            end
            
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Performace tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   for opE = 'N' compared three variants:
%   Var 1:
%                 Db = eqn.E_ * [B; zeros(na - nd, cols)];
%                 C1 = B(1 : nd, :) - Db(1 : nd, :);
%   Var 2:
%                 Db = eqn.E_ * [B(1 : nd, :); zeros(na, cols)];
%                 C1 = B(1 : nd, :) - Db(1 : nd, :);
%   Var 3:
%                 C1 = B(1 : nd, :) - eqn.E_(1 : nd, 1 : nd) * B(1 : nd, :);
%   --> Variant 1 fastest (Jul. 2013 MATLAB R2012a 7.14.0.739)
