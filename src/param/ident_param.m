function out = ident_param(varargin)

%% declare persistent variables
persistent p


%% define parameter set to be identified

if isempty(p)

variableNames = {
    % parameters to identify
    {'delta_m';...
    'delta_c_g';...
    'i_0_ref_c';...
    'alpha_c';...
    'k_lambda';...
    'k_D_w';...
    'k_kappa';...
    'k_t'}
    };

set_ident.variableNames=variableNames;

variableDimensions = {
    % parameters to identify
    [1 1;
    1 1;
    1 1;
    1 1;
    1 1;
    1 1;
    1 1;
    1 1]
    };

set_ident.variableDimensions = variableDimensions;

par_lb = { [0.5*0.012134e-3;
            0.5*0.151e-3;
            0.5*1e-5;
            0.5*0.5;
            0.5*1;
            0.5*1;
            0.5*1;
            0.5*1]
    };

set_ident.lowerBounds = par_lb;

par_ub = {
    % adaptable parameters
            [1.5*0.012134e-3;
            1.5*0.151e-3;
            1.5*1e-5;
            1.5*0.5;
            1.5*1;
            1.5*1;
            1.5*1;
            1.5*1]
    };

set_ident.upperBounds = par_ub;

p.set_ident = set_ident;

allLowerBounds = vertcat(par_lb{:});
allUpperBounds = vertcat(par_ub{:});
allVariableNames=cell(1,0);
for i1=1:numel(set_ident.variableNames)
    for i2=1:numel(set_ident.variableNames{i1})
        allVariableNames{end+1} = set_ident.variableNames{i1}{i2};
    end
end
allVariableDimensions = vertcat(set_ident.variableDimensions{:});

p.npar = sum(prod(allVariableDimensions,2));

p.phys2scaled = @(x) rescaling(x, allLowerBounds, allUpperBounds, 0, 1);

p.scaled2phys = @(x) rescaling(x, 0, 1, allLowerBounds, allUpperBounds);

p.par2struct = @(xVec) vec2struct(xVec,variableNames,variableDimensions);

end

%%
out = p;

end
