function loss = loss(beta, dataTbl, theta, tau, h, lagIndepVar)
%LOSS Compute the loss function for Quantile-on-Quantile Regression (QQR)
%
%   loss = qqr.loss(beta, dataTbl, theta, tau, h, lagIndepVar) computes the
%   weighted check function loss for QQR estimation as defined in Equation 6
%   of Sim & Zhou (2015).
%
%   Reference:
%       Sim, N. & Zhou, H. (2015). "Oil prices, US stock return, and the
%       dependence between their quantiles." Journal of Banking & Finance,
%       55, 1-8. https://doi.org/10.1016/j.jbankfin.2015.01.013
%
%   Inputs:
%       beta        - Coefficient vector [beta0; beta1; alpha] where:
%                     beta0 is the intercept
%                     beta1 is the slope for the independent variable
%                     alpha is the coefficient for the lagged dependent variable
%       dataTbl     - Table containing the data. The dependent variable should
%                     be in the last column, and independent variables in the
%                     remaining columns.
%       theta       - Quantile level for the dependent variable (0 < theta < 1)
%       tau         - Quantile level for the independent variable (0 < tau < 1)
%       h           - Bandwidth parameter for the Gaussian kernel (h > 0)
%       lagIndepVar - Logical flag indicating whether to lag the independent
%                     variable (true) or not (false)
%
%   Output:
%       loss        - Scalar value of the weighted check function loss
%
%   Example:
%       % Create sample data
%       X = randn(100, 1);
%       Y = 0.5 + 0.3 * X + randn(100, 1) * 0.1;
%       dataTbl = table(X, Y);
%
%       % Initial coefficients
%       beta0 = [0; 0; 0];
%
%       % Compute loss at theta=0.5, tau=0.5
%       L = qqr.loss(beta0, dataTbl, 0.5, 0.5, 0.05, false);
%
%   See also qqr.estimate, qqr.plot, fminunc, quantile

% Input validation
arguments
    beta (:,1) {mustBeNumeric, mustBeReal}
    dataTbl table
    theta (1,1) {mustBeNumeric, mustBeReal, mustBeGreaterThan(theta, 0), mustBeLessThan(theta, 1)}
    tau (1,1) {mustBeNumeric, mustBeReal, mustBeGreaterThan(tau, 0), mustBeLessThan(tau, 1)}
    h (1,1) {mustBeNumeric, mustBeReal, mustBePositive}
    lagIndepVar (1,1) logical
end

% Extract dependent variable (last column) and independent variables (remaining columns)
dependentVar = dataTbl{:, end};
IndependentVar = dataTbl{:, 1:end-1};
alpha = beta(end);

% Apply lag to independent variable if specified
if lagIndepVar
    IndependentVar = [NaN; IndependentVar(1:end-1)];
end

% Create lagged dependent variable
laggedDependentVar = [NaN; dependentVar(1:end-1)];

% Build regression table and remove missing values
regTbl = rmmissing(table(IndependentVar, laggedDependentVar, dependentVar, ...
    VariableNames={'IndependentVar', 'LaggedDependentVar', 'DependentVar'}));

% Compute fitted values
fit = beta(1) + beta(2) * (regTbl.IndependentVar - quantile(regTbl.IndependentVar, tau)) + ...
    alpha * regTbl.LaggedDependentVar;

% Compute residuals and apply check function
residual = regTbl.DependentVar - fit;
residual(residual > 0) = residual(residual > 0) * theta;
residual(residual < 0) = residual(residual < 0) * (theta - 1);

% Compute kernel density weights
Idx = tiedrank(regTbl.IndependentVar, 1);
F = (Idx - 1) / height(regTbl);
K = normpdf((F - tau) / h);

% Compute weighted loss
loss = sum(residual .* K);
end
