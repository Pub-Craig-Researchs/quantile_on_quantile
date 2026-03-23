function result = estimate(dataTbl, options)
%ESTIMATE Estimate Quantile-on-Quantile Regression (QQR) model
%
%   result = qqr.estimate(dataTbl) estimates the QQR model using default
%   settings, based on the methodology of Sim & Zhou (2015).
%
%   result = qqr.estimate(dataTbl, Name, Value) specifies additional options
%   using one or more name-value pair arguments.
%
%   Reference:
%       Sim, N. & Zhou, H. (2015). "Oil prices, US stock return, and the
%       dependence between their quantiles." Journal of Banking & Finance,
%       55, 1-8. https://doi.org/10.1016/j.jbankfin.2015.01.013
%
%   Inputs:
%       dataTbl     - Table containing data. The dependent variable should be
%                     in the last column, and independent variables in the
%                     remaining columns.
%
%   Name-Value Pair Arguments:
%       LagIndepVar   - Logical flag to lag the independent variable.
%                       Default: false
%       ThetaGrid     - Row vector of quantile levels for the dependent variable.
%                       Default: linspace(0.05, 0.95, 19)
%       TauGrid       - Column vector of quantile levels for the independent variable.
%                       Default: linspace(0.05, 0.95, 19)'
%       Bandwidth     - Bandwidth for kernel smoothing. Use 'auto' for automatic
%                       selection via ksdensity plug-in method, or a positive scalar.
%                       Default: 'auto'
%       MinBandwidth  - Minimum bandwidth value when using 'auto'.
%                       Default: 0.05
%       LowerBound    - Lower bounds for optimization parameters [beta0, beta1, alpha].
%                       Default: [] (no constraint)
%       UpperBound    - Upper bounds for optimization parameters [beta0, beta1, alpha].
%                       Default: [] (no constraint)
%       InitialPoint  - Initial point for optimization [beta0, beta1, alpha].
%                       Default: [0, 0, 0]
%       Display       - Logical flag to display progress during estimation.
%                       Default: false
%
%   Output:
%       result        - Structure containing estimation results:
%           .intercept   - Intercept estimates (nTheta x nTau matrix)
%           .slope       - Slope estimates (nTheta x nTau matrix)
%           .alpha       - Lagged dependent variable coefficients (nTheta x nTau matrix)
%           .theta       - Quantile levels for dependent variable (1 x nTheta)
%           .tau         - Quantile levels for independent variable (nTau x 1)
%           .bandwidth   - Actual bandwidth used (scalar)
%           .varNames    - Variable names from input table
%           .lagIndepVar - Whether independent variable was lagged
%           .options     - Complete options structure
%
%   Example:
%       % Create sample data
%       rng(42);
%       n = 500;
%       X = randn(n, 1);
%       Y = 0.5 + 0.8 * X + randn(n, 1) * 0.2;
%       dataTbl = table(X, Y);
%
%       % Basic estimation with default settings
%       result = qqr.estimate(dataTbl);
%
%       % Estimation with lagged independent variable
%       result = qqr.estimate(dataTbl, 'LagIndepVar', true);
%
%       % Estimation with custom grids and bounded optimization
%       result = qqr.estimate(dataTbl, ...
%           'ThetaGrid', linspace(0.1, 0.9, 9), ...
%           'TauGrid', linspace(0.1, 0.9, 9)', ...
%           'LowerBound', [-1, -1, -1], ...
%           'UpperBound', [1, 1, 1], ...
%           'Display', true);
%
%       % Plot slope estimates
%       contourf(result.tau, result.theta, result.slope, 20);
%       colorbar;
%
%   See also qqr.loss, qqr.plot, fmincon, ksdensity

% Input validation using arguments block
arguments
    dataTbl table
    options.LagIndepVar (1,1) logical = false
    options.ThetaGrid double = linspace(0.05, 0.95, 19)
    options.TauGrid double = linspace(0.05, 0.95, 19)'
    options.Bandwidth = 'auto'           % 'auto' or positive scalar
    options.MinBandwidth (1,1) double = 0.05
    options.LowerBound double = []       % parameter lower bounds, default unconstrained
    options.UpperBound double = []       % parameter upper bounds, default unconstrained
    options.InitialPoint (1,3) double = [0, 0, 0]
    options.Display (1,1) logical = false
end

% Ensure grids are column vectors for consistency
options.ThetaGrid = options.ThetaGrid(:)';
options.TauGrid = options.TauGrid(:);

% Set random seed for reproducibility
rng(0, 'twister');

% Extract grids
theta = options.ThetaGrid;
tau = options.TauGrid;
nTheta = length(theta);
nTau = length(tau);

% Compute bandwidth
if ischar(options.Bandwidth) || isstring(options.Bandwidth)
    if strcmpi(options.Bandwidth, 'auto')
        % Automatic bandwidth selection using ksdensity plug-in method
        if options.LagIndepVar
            [~, ~, h] = ksdensity(dataTbl{1:end-1, 1}, 'Bandwidth', 'plug-in');
        else
            [~, ~, h] = ksdensity(dataTbl{:, 1}, 'Bandwidth', 'plug-in');
        end
        h = max(h, options.MinBandwidth);
    else
        error('qqr:estimate:InvalidBandwidth', ...
            "Bandwidth must be 'auto' or a positive scalar.");
    end
else
    % Use user-specified bandwidth
    h = options.Bandwidth;
    if ~isscalar(h) || h <= 0
        error('qqr:estimate:InvalidBandwidth', ...
            'Bandwidth must be a positive scalar.');
    end
end

% Pre-allocate coefficient matrices
itcpt = zeros(nTheta, nTau);
slp = zeros(nTheta, nTau);
alp = zeros(nTheta, nTau);

% Set optimization options
fminconOpts = optimoptions('fmincon', 'Display', 'off');

% Initial point and bounds
x0 = options.InitialPoint;
lb = options.LowerBound;
ub = options.UpperBound;

% Display header if requested
if options.Display
    fprintf('QQR Estimation: %d theta values x %d tau values\n', nTheta, nTau);
    fprintf('Bandwidth: %.4f\n', h);
    fprintf('Progress:\n');
end

% Main estimation loop
for i = 1:nTheta
    for j = 1:nTau
        % Define loss function handle
        funcHdl = @(beta) qqr.loss(beta, dataTbl, theta(i), tau(j), h, options.LagIndepVar);
        
        % Optimize using fmincon
        [solution, ~] = fmincon(funcHdl, x0, [], [], [], [], lb, ub, [], fminconOpts);
        
        % Store results
        itcpt(i, j) = solution(1);
        slp(i, j) = solution(2);
        alp(i, j) = solution(3);
    end
    
    % Display progress if requested
    if options.Display
        fprintf('  theta(%d/%d) = %.2f completed\n', i, nTheta, theta(i));
    end
end

% Build result structure
result.intercept = itcpt;
result.slope = slp;
result.alpha = alp;
result.theta = theta;
result.tau = tau;
result.bandwidth = h;
result.varNames = dataTbl.Properties.VariableNames;
result.lagIndepVar = options.LagIndepVar;
result.options = options;

if options.Display
    fprintf('Estimation complete.\n');
end

end
