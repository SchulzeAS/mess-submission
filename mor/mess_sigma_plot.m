function [out, eqn, opts, oper] = mess_sigma_plot(argone, opts, varargin)
% Computation of simple sigma-magnitude-plots for descriptor systems with
% invertible E and comparison to reduced order models.
%
% Calling sequence:
%
%   [out, eqn, opts, oper] = mess_sigma_plot(eqn, opts, oper, ROM)
% or
%   [out, eqn, opts, oper] = mess_sigma_plot(tr, opts, ROM)
%
% INPUTS:
%   eqn                 struct contains data for equations
%
%   g                   presampled transfer function of the original system
%                       fitting the parameters in opts.sigma
%
%   opts                struct contains parameters for the algorithm
%                       (mandatory with substructure opts.sigma)
%
%   oper                struct contains function handles for operation
%                       with A and E
%
%   ROM                 structure containing reduced order model matrices
%                       either E, A, B, C, D,
%                       or M, E, K, B, Cv, Cp, D
%                       where in the first case E and in both cases D are
%                       optional.
%                       (optional when eqn is present; mandatory otherwise)
%                       
%                       OR a cell array with multiple ROMs as decribed
%                       above. e.g {ROM,ROM} or {ROM;ROM} dimensions either
%                       1xN or Nx1
%
% The relevant substructure of opts is opts.sigma with members:
%
% fmin, fmax   left and right bounds of the frequency range. They will be
%              interpreted as exponents in the logarithmic range if
%              integers are passed.
%              (mandatory for input eqn, ignored for input tr (check w
%              below))
%
% nsample      number of transfer function samples to take in the plot
%              (optional, defaults to 100)
%
% info         verbosity control. (optional, defaults to 2)
%                1   only progress is reported
%                2   also generate the actual sigma plots.
%
%  w           vector of frequency samples fitting
%              (mandatory for input tr, ignored for input eqn)
%
% Outputs:
%
% out.w         vector of sampling frequencies
%
% out.err       vector of sampled maximal singular values of the transfer
%               function of the error system (only when ROM was given)

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

%% Check and assign inputs

narginchk(2,4);
if not(isa(opts,'struct'))
    error('second input must be am options structure');
end
disp(nargin);
if nargin==2
    if not(isa(argone,'numeric'))
        error(['input 1 must be numeric (3d array of doubles)',...
            ' in the 2 input case']);
    else
        g = argone;
        ROM = [];
    end
end

if (nargin > 2)
    if (not(isa(argone, 'struct')) && not(isa(argone, 'numeric')) && not(isa(argone,'cell')))
        error(['input 1 must be numeric (3d array of doubles)',...
            ' or equation structure']);
    end
   
    if not(isa(varargin{1}, 'struct')) || (nargin == 4 && (not(isa(varargin{2}, 'struct'))) && (nargin == 4 && (not(isa(varargin{2}, 'cell')))))
        error(['Either all inputs are structures, or the first input',...
            'is numeric and the rest are structures'])
    end
    if isa(argone, 'numeric')
        g = argone;
        ROM = varargin{1};
        eqn = [];
        oper = [];
    else
        g=[];
        eqn = argone;
        oper = varargin{1};
        if nargin ==4
            ROM = varargin{2};
        else
            ROM = [];
        end
    end
end
%%verify dimensions of @var (ROM) transpose in case of wrong vector orientation
if iscell(ROM)
    dimRom = size(ROM);
    if dimRom(1)~= 1
        ROM = ROM.';
    end
    dimRom = size(ROM);
    if dimRom(1) ~= 1
        error('ROMs not in a 1xN structure');
    end
end
%% check field opts.sigma
if not(isfield(opts,'sigma')) || not(isstruct(opts.sigma))
    error('MESS:control_data',...
        'No sigma plot control data found in options structure.');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check info parameter for output verbosity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isfield(opts.sigma,'info'))
    opts.sigma.info=2;
else
    if not(isnumeric(opts.sigma.info)) && not(islogical(opts.sigma.info))
        error('MESS:info',...
            'opts.sigma.info parameter must be logical or numeric.');
    end
end

if isempty(g) && ( not(isfield(opts.sigma,'w')) && ...
        (not(isfield(opts.sigma,'fmin')) || not(isfield(opts.sigma,'fmax'))))
    error('MESS:control_data',...
        'sigma plot control data does not contain frequency rang bounds.');
end

if not(isempty(g)) && not(isfield(opts.sigma,'w'))
    error('MESS:control_data',...
        ['sigma plot control data must contain frequency vector w', ...
        ' when presampled original is given']);
end

if isempty(g) && not(isfield(opts.sigma,'nsample'))
    opts.sigma.nsample = 100;
end

if isfield(opts.sigma,'w')
    w = opts.sigma.w;
    opts.sigma.nsample = length(w);
else
    if (floor(opts.sigma.fmin)==opts.sigma.fmin) && ...
            (floor(opts.sigma.fmax)==opts.sigma.fmax)
        w=logspace(opts.sigma.fmin,opts.sigma.fmax,opts.sigma.nsample);
    else
        w=logspace(log10(opts.sigma.fmin), ...
            log10(opts.sigma.fmax),opts.sigma.nsample);
    end
end

%% preallocation

szRom = size(ROM);
szRom = szRom(2);


for i = 1:szRom
    out(i).tr1=zeros(1,opts.sigma.nsample);
    if not(isempty(ROM))
        out(i).tr2=out.tr1;
        out(i).err=out.tr1;
        out(i).relerr=out.tr1;
    end
    
    if not(isempty(g))
        [m,p] = size(g(:,:,1));
    else
        m = size(eqn.C,1);
        p = size(eqn.B,2);
    end
    if not(isempty(g))
        out(i).g1 = g;
    else
        out(i).g1 = zeros(m,p,opts.sigma.nsample);
    end
    if not(isempty(ROM(i)))
        out(i).g2 = zeros(m,p,opts.sigma.nsample);
    end
    if opts.sigma.info
        fprintf(['Computing TFMs of original and reduced order systems and ' ...
            'MOR errors\n']);
    end
end
%% preprocess shifted solver if eqn is given
if isempty(g)
    [result, eqn, opts, oper] = oper.init(eqn, opts, oper, 'A','E');
    if not(result)
        error('MESS:control_data', ['system data is not completely',...
            ' defined or corrupted']);
    end
    
    [eqn, opts, oper] = oper.sol_ApE_pre(eqn, opts, oper);
end
%% make sure we have an E in the first-order ROM



if not(isempty(ROM)) && not(iscell(ROM))
    if (not(isfield(ROM,'A')) && ...
            (not(isfield(ROM,'K')) || not(isfield(ROM,'M')))) || ...
            not(isfield(ROM,'B')) || ...
            (not(isfield(ROM,'C')) && ...
            not(isfield(ROM, 'Cv')) && not(isfield(ROM, 'Cp')))
        error('Found incomplete ROM structure!')
    end
    if not(isfield(ROM,'E'))
        if not(isfield(ROM,'A'))
            ROM.E = [];
        else
            ROM.E = eye(size(ROM.A));
        end
    end
    if (isfield(ROM, 'Cp') && not(isfield(ROM, 'Cv'))), ROM.Cv = []; end
    if (isfield(ROM, 'Cv') && not(isfield(ROM, 'Cp'))), ROM.Cp = []; end
end
%added to deal with cells
%basically same thing as before, just dealing with multiple ROMs
if not(isempty(ROM)) && iscell(ROM)
    for c = 1:szRom
        if (not(isfield(ROM{1,c},'A')) && ...
                (not(isfield(ROM{1,c},'K')) || not(isfield(ROM{1,c},'M')))) || ... %
                not(isfield(ROM{1,c},'B')) || ...
                (not(isfield(ROM{1,c},'C')) && ...
                not(isfield(ROM{1,c}, 'Cv')) && not(isfield(ROM{1,c}, 'Cp')))
            error('Found incomplete ROM structure!')
        end
        if not(isfield(ROM{1,c},'E'))
            if not(isfield(ROM{1,c},'A'))
                ROM{1,c}.E = [];
            else
                ROM{1,c}.E = eye(size(ROM{1,c}.A));
            end
        end
        if (isfield(ROM{1,c}, 'Cp') && not(isfield(ROM{1,c}, 'Cv'))), ROM{1,c}.Cv = []; end
        if (isfield(ROM{1,c}, 'Cv') && not(isfield(ROM{1,c}, 'Cp'))), ROM{1,c}.Cp = []; end
    end
end


%% perform the actual sampling
if not(iscell(ROM))
    for k=1:opts.sigma.nsample
        if (opts.sigma.info && not(mod(k,opts.sigma.nsample/10)))
            fprintf('\r Step %3d / %3d',k,opts.sigma.nsample);
        end
        if isempty(g) % sample original model only if it was not given
            if isfield(eqn,'D')&& not(isempty(eqn.D))
                out.g1(:,:,k) = full(eqn.C *...
                    oper.sol_ApE(eqn, opts,'N',-1i*w(k),'N',-eqn.B,'N') + eqn.D);
            else
                out.g1(:,:,k) = full(eqn.C * ...
                    oper.sol_ApE(eqn, opts,'N',-1i*w(k),'N',-eqn.B,'N'));
            end
        end
        if not(isempty(ROM)) % sample reduced model only if it was given
            if isfield(ROM,'D') && not(isempty(ROM.D))
                if isfield(ROM,'A')
                    out.g2(:,:,k) = ROM.C * ((1i*w(k)*ROM.E-ROM.A) \ ROM.B ) ...
                        + ROM.D;
                else
                    out.g2(:,:,k) = (ROM.Cp + 1i*w(k)*ROM.Cv) * ...
                        (( -w(k) * (w(k)*ROM.M-1i*ROM.E) + ROM.K) \ ROM.B )...
                        + ROM.D;
                end
            else
                if isfield(ROM,'A')
                    out.g2(:,:,k) = ROM.C * ( (1i*w(k)*ROM.E - ROM.A) \ ROM.B );
                else
                    out.g2(:,:,k) = (ROM.Cp + 1i*w(k)*ROM.Cv) *...
                        (( -w(k) * ( w(k)*ROM.M - 1i*ROM.E ) + ROM.K) \ ROM.B );
                end
            end
        end
        out.tr1(k) = max(svd(out.g1(:,:,k)));
        if not(isempty(ROM))
            out.err(k) = max(svd(out.g1(:,:,k)-out.g2(:,:,k)));
            out.tr2(k) = max(svd(out.g2(:,:,k)));
            out.relerr(k)=out.err(k)/out.tr1(k);
        end
    end
    
    out.w = w;
end
if iscell(ROM)
    for c = 1:szRom
        for k=1:opts.sigma.nsample
            if (opts.sigma.info && not(mod(k,opts.sigma.nsample/10)))
                fprintf('\r Step %3d / %3d',k,opts.sigma.nsample);
            end
            if isempty(g) % sample original model only if it was not given
                if isfield(eqn,'D')&& not(isempty(eqn.D))
                    out(c).g1(:,:,k) = full(eqn.C *...
                        oper.sol_ApE(eqn, opts,'N',-1i*w(k),'N',-eqn.B,'N') + eqn.D);
                else
                    out(c).g1(:,:,k) = full(eqn.C * ...
                        oper.sol_ApE(eqn, opts,'N',-1i*w(k),'N',-eqn.B,'N'));
                end
            end
            if not(isempty(ROM{1,c})) % sample reduced model only if it was given
                if isfield(ROM{1,c},'D') && not(isempty(ROM{1,c}.D))
                    if isfield(ROM{1,c},'A')
                        out(c).g2(:,:,k) = ROM{1,c}.C * ((1i*w(k)*ROM{1,c}.E-ROM{1,c}.A) \ ROM{1,c}.B ) ...
                            + ROM{1,c}.D;
                    else
                        out(c).g2(:,:,k) = (ROM{1,c}.Cp + 1i*w(k)*ROM{1,c}.Cv) * ...
                            (( -w(k) * (w(k)*ROM{1,c}.M-1i*ROM{1,c}.E) + ROM{1,c}.K) \ ROM{1,c}.B )...
                            + ROM{1,c}.D;
                    end
                else
                    if isfield(ROM{1,c},'A')
                        out(c).g2(:,:,k) = ROM{1,c}.C * ( (1i*w(k)*ROM{1,c}.E - ROM{1,c}.A) \ ROM{1,c}.B );
                    else
                        out(c).g2(:,:,k) = (ROM{1,c}.Cp + 1i*w(k)*ROM{1,c}.Cv) *...
                            (( -w(k) * ( w(k)*ROM{1,c}.M - 1i*ROM{1,c}.E ) + ROM{1,c}.K) \ ROM{1,c}.B );
                    end
                end
            end
            out(c).tr1(k) = max(svd(out(c).g1(:,:,k)));
            if not(isempty(ROM{1,c}))
                out(c).err(k) = max(svd(out(c).g1(:,:,k)-out(c).g2(:,:,k)));
                out(c).tr2(k) = max(svd(out(c).g2(:,:,k)));
                out(c).relerr(k)=out(c).err(k)/out(c).tr1(k);
            end
        end
        out(c).w = w;
    end
    
    
end
if opts.sigma.info
    fprintf('\n\n');
end
%% postprocess shifted solver if eqn is given
if isempty(g)
    [eqn, opts, oper] = oper.sol_ApE_post(eqn, opts, oper);
end

%% Finally,  the plots (if desired)
if (isnumeric(opts.sigma.info) && opts.sigma.info > 1) && not(iscell(ROM))
    if not(isempty(ROM))
        figure;
        subplot(2,1,1);
        loglog(out.w, out.err);
        title('absolute model reduction error');
        xlabel('\omega');
        ylabel('\sigma_{max}(G(j\omega) - G_r(j\omega))');
        axis tight;
        hold on;
        subplot(2,1,2);
        loglog(out.w, out.relerr);
        title('relative model reduction error');
        xlabel('\omega');
        ylabel(['\sigma_{max}(G(j\omega) - G_r(j\omega)) / \' ...
            'sigma_{max}(G(j\omega))']);
        axis tight;
        hold off;
    end
    
    figure;
    loglog(out.w, out.tr1);
    if not(isempty(ROM))
        hold on;
        loglog(out.w, out.tr2, 'r--');
        legend({'original system','reduced system'});
        title('Transfer functions of original and reduced systems');
    else
        legend({'original system'});
        title('Transfer functions of original system');
    end
    xlabel('\omega');
    ylabel('\sigma_{max}(G(j\omega))');
    
    axis tight;
    hold off;
end
%%additional case if ROM iscell -> multiple ROMs
if (isnumeric(opts.sigma.info) && opts.sigma.info > 1) && iscell(ROM)
    
    
    if not(isempty(ROM))
        figure;            
        for i =1:szRom
            hold on;
            subplot(2,1,1);
            loglog(out(i).w, out(i).err);
            title('absolute model reduction error');
            xlabel('\omega');
            ylabel('\sigma_{max}(G(j\omega) - G_r(j\omega))');
            axis tight;
        end
        hold off;
        for i =1:szRom
            hold on;
            subplot(2,1,2);
            loglog(out(i).w, out(i).relerr);
            title('relative model reduction error');
            xlabel('\omega');
            ylabel(['\sigma_{max}(G(j\omega) - G_r(j\omega)) / \' ...
                'sigma_{max}(G(j\omega))']);
            axis tight;
           
        end
        hold off;
    end
        
        figure;
        loglog(out.w, out.tr1);
        if not(isempty(ROM))
            hold on;
            loglog(out.w, out.tr2, 'r--');
            legend({'original system','reduced system'});
            title('Transfer functions of original and reduced systems');
        else
            legend({'original system'});
            title('Transfer functions of original system');
        end
        xlabel('\omega');
        ylabel('\sigma_{max}(G(j\omega))');
        
        axis tight;
        hold off;
    
end

