function fig = plotSurface(result, options)
% PLOTSURFACE Plot QQR estimation results as a contour surface
%
%   fig = qqr.plotSurface(result) plots the slope coefficients from the QQR
%   estimation result as a filled contour plot and returns the figure handle.
%
%   fig = qqr.plotSurface(result, Name, Value) specifies additional options using
%   one or more name-value pair arguments.
%
% Input Arguments:
%   result - Structure returned by qqr.estimate containing:
%       .slope     - Slope coefficient matrix (theta x tau)
%       .intercept - Intercept coefficient matrix (theta x tau)
%       .alpha     - Alpha coefficient matrix (theta x tau)
%       .tau       - Quantile levels for independent variable
%       .theta     - Quantile levels for dependent variable
%       .varNames  - Cell array of variable names (first = indep, last = dep)
%       .lagIndepVar - Logical indicating if lagged independent variable
%
% Name-Value Arguments:
%   'Parameter'       - Which parameter to plot: 'slope' (default), 
%                       'intercept', or 'alpha'
%   'Subplot'         - Subplot position as [rows, cols, index]. Default [1,1,1]
%   'Colormap'        - Colormap name. Default 'gray'
%   'NumContours'     - Number of contour levels. Default 20
%   'ReverseColormap' - Reverse colormap direction. Default true
%
% Example:
%   result = qqr.estimate(dataTbl);
%   figure;
%   qqr.plotSurface(result, 'Parameter', 'slope');
%
%   % Multiple subplots
%   figure;
%   qqr.plotSurface(result, 'Parameter', 'slope', 'Subplot', [1,2,1]);
%   qqr.plotSurface(result, 'Parameter', 'intercept', 'Subplot', [1,2,2]);
%
% See also: qqr.estimate

    arguments
        result struct
        options.Parameter (1,:) char {mustBeMember(options.Parameter, {'slope','intercept','alpha'})} = 'slope'
        options.Subplot (1,3) double = [1, 1, 1]   % [rows, cols, index]
        options.Colormap (1,:) char = 'gray'
        options.NumContours (1,1) double = 20
        options.ReverseColormap (1,1) logical = true
    end

    % Get current figure handle for return value
    fig = gcf;

    % Select data matrix based on parameter choice
    switch options.Parameter
        case 'slope'
            data = result.slope;
        case 'intercept'
            data = result.intercept;
        case 'alpha'
            data = result.alpha;
    end

    % Configure colormap
    c = colormap(options.Colormap);
    if options.ReverseColormap
        c = c(end:-1:1, :);
        colormap(c);
    end

    % Create subplot and contour plot
    subplot(options.Subplot(1), options.Subplot(2), options.Subplot(3))
    contourf(result.tau, result.theta, data, options.NumContours, ...
        'EdgeColor', 'none', 'FaceColor', 'flat')

    % Generate title with corrected variable order
    % Note: varNames{end} = dependent variable (last column)
    %       varNames{1}   = independent variable (first column)
    depVarName = result.varNames{end};    % 因变量 = 最后一列
    indVarName = result.varNames{1};      % 自变量 = 第一列

    % Try to parse simplified names (e.g., "Ret_XXX_future" -> "XXX")
    % Fall back to original names if pattern doesn't match
    try
        depShort = extractBefore(extractAfter(depVarName, "Ret_"), "_future");
        if isempty(depShort)
            depShort = depVarName;
        end
    catch
        depShort = depVarName;
    end
    
    try
        indShort = extractBefore(extractAfter(indVarName, "Ret_"), "_future");
        if isempty(indShort)
            indShort = indVarName;
        end
    catch
        indShort = indVarName;
    end

    % Create title string with correct order: dependent ~ independent
    if result.lagIndepVar
        titleStr = depShort + " ~ lagged " + indShort;
    else
        titleStr = depShort + " ~ " + indShort;
    end
    title(titleStr, 'Interpreter', 'tex', 'FontSize', 10)

    % Axis labels consistent with paper notation
    xlabel('independent variable \tau')
    ylabel('dependent variable \theta')
    colorbar()
    drawnow

end
