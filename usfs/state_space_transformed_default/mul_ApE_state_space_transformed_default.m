function C = mul_ApE_state_space_transformed_default ...
    (eqn, opts, opA, p, opE, B, opB)
%% function C = mul_ApE_state_space_transformed_default ...
%    (eqn, opts, opA, p, opE, B, opB)
%
% This function returns C = EL\(A_ + pE_)/EU*B, where matrices A_ and E_
% given by structure eqn and input matrix B could be transposed.
% Matrices A_ and E_ are assumed to be quadratic.
%
% Inputs
%   eqn             struct contains data for equations
%
%   opts            struct contains parameters for the algorithm
%
%   opA             character specifying the shape of A_
%                   opA = 'N' performs EL\(A_ + p*op(E_))/EU * opB(B)
%                   opA = 'T' performs EL\(A_' + p*op(E_))/EU * opB(B)
%
%   p               scalar value
%
%   opE             character specifying the shape of E_
%                   opE = 'N' solves EL\(op(A_) + p*E_)/EU * opB(B)
%                   opE = 'T' solves EU'\(op(A_) + p*E_')/EL' * opB(B)
%
%   B               m-x-p matrix
%
%   opB             character specifying the shape of B
%                   opB = 'N' performs EL\(op(A_) + p*op(E_))/EU * B
%                   opB = 'T' performs EL\(op(A_) + p*op(E_))/EU * B'
%
% Output
%	C = EL\(op(A_) + p*op(E_))/EU * op(B)
%
% This function uses another default function size_default(eqn, opts) to
% obtain the number of rows of matrix A_ in structure eqn.

%% Check input parameters.
assert(ischar(opA) && ischar(opE) && ischar(opB), ...
    'MESS:error_arguments', ...
    'opA or opB is not a char');

opA = upper(opA);
opE = upper(opE);
opB = upper(opB);

assert((opA == 'N') || (opA == 'T'), ...
    'MESS:error_arguments', ...
    'opA is not ''N'' or ''T''');

assert((opE == 'N') || (opE == 'T'), ...
    'MESS:error_arguments', ...
    'opE is not ''N'' or ''T''');

assert((opB == 'N') || (opB == 'T'), ...
    'MESS:error_arguments', ...
    'opB is not ''N'' or ''T''');

assert(isnumeric(p) && (length(p) == 1), ...
    'MESS:error_arguments', ...
    'p is not a numeric scalar');

assert(isnumeric(B) && ismatrix(B), ...
    'MESS:error_arguments', ...
    'B has to ba a matrix');

%% Check data in eqn structure.
assert(isfield(eqn,'A_'), ...
    'MESS:error_arguments', ...
    'field eqn.A_ is not defined');

assert(isfield(eqn,'E_'), ...
    'MESS:error_arguments', ...
    'field eqn.E_ is not defined');

if isfield(eqn, 'haveE') && eqn.haveE
    assert(isfield(eqn, 'EL'), ...
        'MESS:error_arguments', ...
        'field eqn.EL is not defined');
    assert(isfield(eqn, 'EU'), ...
        'MESS:error_arguments', ...
        'field eqn.EU is not defined');
else
    eqn.haveE = 0;
end

rowA = size_default(eqn, opts);
colA = rowA;

%% Perform multiplication.
if eqn.haveE % Case of non-identity E matrix.
    switch opA
        case 'N'
            switch opE
                case 'N'
                    switch opB
                        case 'N' % Implement EL\(A_ + pE_)/EU*B.
                            assert(colA == size(B, 1), ...
                                'MESS:error_arguments', ...
                                ['number of columns of A_ differs ' ...
                                'with number of rows of B']);
                            C = eqn.EL \ ((eqn.A_ + p*eqn.E_) ...
                                * (eqn.EU \ B));
                        case 'T' % Implement EL\(A_ + pE_)/EU*B'.
                            assert(colA == size(B, 2), ...
                                'MESS:error_arguments', ...
                                ['number of columns of A_ differs ' ...
                                'with number of columns of B']);
                            C = eqn.EL \ ((eqn.A_ + p*eqn.E_) ...
                                * (eqn.EU \ B'));
                    end
                    
                case 'T'
                    switch opB
                        case 'N' % Implement EU'\(A_ + pE_')/EL'*B.
                            assert(colA == size(B, 1), ...
                                'MESS:error_arguments', ...
                                ['number of columns of A_ differs ' ...
                                'with number of rows of B']);
                            C = eqn.EU' \ ((eqn.A_ + p*eqn.E_') ...
                                * (eqn.EL' \ B));
                        case 'T' % Implement EU'\(A_ + pE_')/EL'*B'.
                            assert(colA == size(B, 2), ...
                                'MESS:error_arguments', ...
                                ['number of columns of A_ differs ' ...
                                'with number of columns of B']);
                            C = eqn.EU' \ ((eqn.A_ + p*eqn.E_') ...
                                * (eqn.EL' \ B'));
                    end
                    
            end
        case 'T'
            switch opE
                case 'N'
                    switch opB
                        case 'N' % Implement EL\(A_' + pE_)/EU*B.
                            assert(rowA == size(B, 1), ...
                                'MESS:error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                'number rows of B']);
                            C = eqn.EL \ ((eqn.A_' + p*eqn.E_) ...
                                * (eqn.EU \ B));
                        case 'T' % Implement EL\(A_' + pE_)/EU*B'.
                            assert(rowA == size(B, 2), ...
                                'MESS:error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                'number of columns of B']);
                            C = eqn.EL \ ((eqn.A_' + p*eqn.E_) ...
                                * (eqn.EU \ B'));
                    end
                    
                case 'T'
                    switch opB
                        case 'N' % Implement EU'\(A_' + pE_')/EL'*B.
                            assert(rowA == size(B, 1), ...
                                'MESS:error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                'number rows of B']);
                            C = eqn.EU' \ ((eqn.A_' + p*eqn.E_') ...
                                * (eqn.EL' \ B));
                        case 'T' % Implement EU'\(A_' + pE_')/EL'*X = B'.
                            assert(rowA == size(B, 2), ...
                                'MESS:error_arguments', ...
                                ['number of rows of A_ differs with ' ...
                                'number of columns of B']);
                            C = eqn.EU' \ ((eqn.A_' + p*eqn.E_') ...
                                * (eqn.EL' \ B'));
                    end
                    
            end
    end
else % Case of E_ = I_n, was set by init.
    switch opA
        case 'N'
            switch opB
                case 'N' % Implement (A_ + pI_n)*B.
                    assert(colA == size(B, 1), ...
                        'MESS:error_arguments', ...
                        ['number of columns of A_ differs with ' ...
                        'number of rows of B']);
                    C = (eqn.A_ + p*eqn.E_) * B;
                case 'T' % Implement (A_ + pI_n)*B'.
                    assert(colA == size(B, 2), ...
                        'MESS:error_arguments', ...
                        ['number of columns of A_ differs with ' ...
                        'number of columns of B']);
                    C = (eqn.A_ + p*eqn.E_) * B';
            end
            
        case 'T'
            switch opB
                case 'N' % Implement (A_' + pE_)*B.
                    assert(rowA == size(B, 1), ...
                        'MESS:error_arguments', ...
                        ['number of rows of A_ differs with ' ...
                        'number rows of B']);
                    C = (eqn.A_' + p*eqn.E_) * B;
                case 'T' % Implement (A_' + pE_)*B'.
                    assert(rowA == size(B, 2), ...
                        'MESS:error_arguments', ...
                        ['number of rows of A_ differs with ' ...
                        'number of columns of B']);
                    C = (eqn.A_' + p*eqn.E_) * B';
            end
            
    end
end

end
