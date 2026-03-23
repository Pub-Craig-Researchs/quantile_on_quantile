function report = validate(result, dataTbl, options)
%VALIDATE Validate QQR estimates against standard quantile regression
%
%   report = qqr.validate(result, dataTbl) validates the Quantile-on-Quantile
%   Regression estimates by comparing the τ-averaged QQR coefficients with
%   standard Quantile Regression (QR) estimates, as described in Equation 10
%   of Sim & Zhou (2015).
%
%   The validation is based on the theoretical property that QQR estimates
%   averaged over τ should equal standard QR estimates:
%       ĉ₀(θ) = (1/S) Σ_s b̂₀(θ, s)  ≈ QR intercept at θ
%       ĉ₁(θ) = (1/S) Σ_s b̂₁(θ, s)  ≈ QR slope at θ
%
%   Reference:
%       Sim, N. & Zhou, H. (2015). "Oil prices, US stock return, and the
%       dependence between their quantiles." Journal of Banking & Finance,
%       55, 1-8. Section 5, Equation 10.
%       https://doi.org/10.1016/j.jbankfin.2015.01.013
%
%   Inputs:
%       result      - Struct containing QQR estimation results with fields:
%                     intercept (nTheta x nTau), slope (nTheta x nTau),
%                     alpha (nTheta x nTau), theta (1 x nTheta)
%       dataTbl     - Table containing the data (same format as qqr.estimate)
%       options     - Name-value arguments:
%           Display - Logical flag to display results (default: true)
%
%   Output:
%       report      - Struct containing validation results:
%           qqMeanIntercept - τ-averaged QQR intercepts (nTheta x 1)
%           qqMeanSlope     - τ-averaged QQR slopes (nTheta x 1)
%           qqMeanAlpha     - τ-averaged QQR alphas (nTheta x 1)
%           qrIntercept     - Standard QR intercepts (nTheta x 1)
%           qrSlope         - Standard QR slopes (nTheta x 1)
%           qrAlpha         - Standard QR alphas (nTheta x 1)
%           theta           - Quantile levels used
%           interceptMAE    - Mean Absolute Error for intercept
%           slopeMAE        - Mean Absolute Error for slope
%           alphaMAE        - Mean Absolute Error for alpha
%
%   Example:
%       % Generate sample data
%       rng(42);
%       X = randn(200, 1);
%       Y = 0.5 + 0.3 * X + randn(200, 1) * 0.2;
%       dataTbl = table(X, Y);
%
%       % Estimate QQR model
%       result = qqr.estimate(dataTbl);
%
%       % Validate estimates
%       report = qqr.validate(result, dataTbl);
%
%       % Suppress display
%       report = qqr.validate(result, dataTbl, Display=false);
%
%   See also qqr.estimate, qqr.loss, qqr.plot

% Input validation
arguments
    result struct
    dataTbl table
    options.Display (1,1) logical = true
end

%% Step 1: Compute τ-averaged QQR estimates
qqMeanIntercept = mean(result.intercept, 2);  % nTheta x 1
qqMeanSlope = mean(result.slope, 2);          % nTheta x 1
qqMeanAlpha = mean(result.alpha, 2);          % nTheta x 1

%% Step 2: Preprocess data (same as qqr.loss)
% Extract dependent variable (last column) and independent variables
dependentVar = dataTbl{:, end};
IndependentVar = dataTbl{:, 1:end-1};

% Check if lagIndepVar was used (from result struct if available)
if isfield(result, 'lagIndepVar') && result.lagIndepVar
    IndependentVar = [NaN; IndependentVar(1:end-1)];
end

% Create lagged dependent variable
laggedDependentVar = [NaN; dependentVar(1:end-1)];

% Build regression table and remove missing values
regTbl = rmmissing(table(IndependentVar, laggedDependentVar, dependentVar, ...
    VariableNames={'IndependentVar', 'LaggedDependentVar', 'DependentVar'}));

% Extract vectors for optimization
Y = regTbl.DependentVar;
X = regTbl.IndependentVar;
lagY = regTbl.LaggedDependentVar;

%% Step 3: Compute standard QR estimates for each theta
theta = result.theta(:);  % Ensure column vector
nTheta = length(theta);

qrIntercept = zeros(nTheta, 1);
qrSlope = zeros(nTheta, 1);
qrAlpha = zeros(nTheta, 1);

% Optimization options
optOptions = optimoptions('fmincon', 'Display', 'off', ...
    'Algorithm', 'interior-point', 'MaxIterations', 1000);

% Initial values and bounds (same as estimate)
beta0 = [0; 0; 0];
lb = [-Inf; -Inf; -Inf];
ub = [Inf; Inf; Inf];

% Estimate standard QR for each theta
for i = 1:nTheta
    thetaVal = theta(i);
    
    % Define objective function for standard QR
    objFun = @(beta) qrLoss(beta, Y, X, lagY, thetaVal);
    
    % Solve optimization
    [betaOpt, ~, exitflag] = fmincon(objFun, beta0, [], [], [], [], lb, ub, [], optOptions);
    
    if exitflag <= 0
        warning('qqr:validate:OptimizationFailed', ...
            'Optimization did not converge for theta = %.3f', thetaVal);
    end
    
    qrIntercept(i) = betaOpt(1);
    qrSlope(i) = betaOpt(2);
    qrAlpha(i) = betaOpt(3);
end

%% Step 4: Compute consistency metrics
interceptMAE = mean(abs(qqMeanIntercept - qrIntercept));
slopeMAE = mean(abs(qqMeanSlope - qrSlope));
alphaMAE = mean(abs(qqMeanAlpha - qrAlpha));

% Compute correlations (handle constant vectors)
if std(qqMeanIntercept) > 0 && std(qrIntercept) > 0
    interceptCorr = corr(qqMeanIntercept, qrIntercept);
else
    interceptCorr = NaN;
end

if std(qqMeanSlope) > 0 && std(qrSlope) > 0
    slopeCorr = corr(qqMeanSlope, qrSlope);
else
    slopeCorr = NaN;
end

if std(qqMeanAlpha) > 0 && std(qrAlpha) > 0
    alphaCorr = corr(qqMeanAlpha, qrAlpha);
else
    alphaCorr = NaN;
end

%% Step 5: Build report struct
report.qqMeanIntercept = qqMeanIntercept;
report.qqMeanSlope = qqMeanSlope;
report.qqMeanAlpha = qqMeanAlpha;
report.qrIntercept = qrIntercept;
report.qrSlope = qrSlope;
report.qrAlpha = qrAlpha;
report.theta = result.theta;
report.interceptMAE = interceptMAE;
report.slopeMAE = slopeMAE;
report.alphaMAE = alphaMAE;
report.interceptCorr = interceptCorr;
report.slopeCorr = slopeCorr;
report.alphaCorr = alphaCorr;

%% Step 6: Display results (optional)
if options.Display
    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('         QQR Validation: Comparison with Standard QR        \n');
    fprintf('============================================================\n');
    fprintf('Reference: Sim & Zhou (2015), Section 5, Equation 10\n\n');
    
    % Display comparison table
    fprintf('%-8s | %-12s %-12s | %-12s %-12s | %-12s %-12s\n', ...
        'Theta', 'QQ_Int', 'QR_Int', 'QQ_Slope', 'QR_Slope', 'QQ_Alpha', 'QR_Alpha');
    fprintf('%s\n', repmat('-', 1, 86));
    
    for i = 1:nTheta
        fprintf('%8.3f | %12.4f %12.4f | %12.4f %12.4f | %12.4f %12.4f\n', ...
            theta(i), qqMeanIntercept(i), qrIntercept(i), ...
            qqMeanSlope(i), qrSlope(i), qqMeanAlpha(i), qrAlpha(i));
    end
    
    fprintf('%s\n', repmat('-', 1, 86));
    fprintf('\n');
    
    % Display summary metrics
    fprintf('Consistency Metrics:\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('  Intercept MAE:  %.6f    Correlation: %.4f\n', interceptMAE, interceptCorr);
    fprintf('  Slope MAE:      %.6f    Correlation: %.4f\n', slopeMAE, slopeCorr);
    fprintf('  Alpha MAE:      %.6f    Correlation: %.4f\n', alphaMAE, alphaCorr);
    fprintf('------------------------------------------------------------\n');
    fprintf('\nNote: Lower MAE and higher correlation indicate better consistency\n');
    fprintf('      between QQR τ-average and standard QR estimates.\n');
    fprintf('============================================================\n\n');
end

end

%% Local function: Standard Quantile Regression Loss
function l = qrLoss(beta, Y, X, lagY, theta)
%QRLOSS Compute the check function loss for standard quantile regression
%
%   This is the standard QR loss without kernel density weighting:
%       L(beta) = Σ ρ_θ(Y - beta(1) - beta(2)*X - beta(3)*lagY)
%   where ρ_θ(u) = u*(θ - I(u<0)) is the check function.

% Compute fitted values
fit = beta(1) + beta(2) * X + beta(3) * lagY;

% Compute residuals
res = Y - fit;

% Apply check function (asymmetric absolute loss)
res(res > 0) = res(res > 0) * theta;
res(res < 0) = res(res < 0) * (theta - 1);

% Sum of losses
l = sum(res);
end
