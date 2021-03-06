function C = mul_E_dae_1(eqn, opts, opE, B, opB)%#ok<INUSL>

%% function mul_A perfoms operation C = opE(E_)*opB(B)
%
% Input:
%   eqn     structure contains field E_
%
%   opts    struct contains parameters for the algorithm
%
%   opE     character specifies the form of opE(E_)
%           opE = 'N' performs E_*opB(B)
%           opE = 'T' performs E_'*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs opE(E_)*B
%           opB = 'T' performs opE(E_)*B'
%
% Output:
% C = opE(E_)*opB(B)
%
%   uses no other dae_1 function

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

%% check input Paramters
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
if(not(isfield(eqn, 'E_'))) || not(isnumeric(eqn.E_))
    error('MESS:error_arguments', ...
        'Missing or Corrupted E_ field detected in equation structure.');
end
if not(isfield(eqn, 'st'))    || not(isnumeric(eqn.st))
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.');
end
st = eqn.st;
%% perfom multiplication
switch opE
    
    case 'N'
        switch opB
            
            %implement operation E_*B
            case 'N'
                if(st ~= size(B, 1))
                    error('MESS:error_arguments','number of cols of E_ differs with rows of B');
                end
                C = eqn.E_(1 : st, 1 : st) * B;
            
            %implement operation E_*B'
            case 'T'
                if(st ~= size(B, 2))
                    error('MESS:error_arguments','number of cols of E_ differs with cols of B');
                end
                C = eqn.E_(1 : st, 1 : st) * B';
        end
        
    case 'T'
        switch opB
            
            %implement operation E_'*B
            case 'N'
                if(st ~= size(B, 1))
                    error('MESS:error_arguments','number of rows of E_ differs with rows of B');
                end
                C = eqn.E_(1 : st, 1 : st)' * B;
                
            %implement operatio E_'*B'
            case 'T'
                if(st ~= size(B, 2))
                    error('MESS:error_arguments','number of rows of E_ differs with cols of B');
                end
                C = eqn.E_(1 : st, 1 : st)' * B';
        end
        
end

end
