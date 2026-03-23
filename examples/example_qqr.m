%% Quantile-on-Quantile Regression (QQR) 完整使用示例
% 基于 Sim & Zhou (2015) 的 QQR 方法
% 本示例展示 +qqr 包的完整工作流

%% 1. 数据准备
% 生成模拟数据（实际使用时替换为真实金融数据）
rng(42, 'twister');
n = 1000;
x = randn(n, 1);         % 自变量（如油价冲击）
y = 0.5 * x + 0.3 * randn(n, 1);  % 因变量（如股票收益率）

dataTbl = table(x, y, 'VariableNames', {'OilShock', 'StockReturn'});

%% 2. QQR 估计
% 使用默认参数（无约束，自适应带宽）
result = qqr.estimate(dataTbl, ...
    'LagIndepVar', false, ...
    'ThetaGrid', linspace(0.1, 0.9, 9), ...
    'TauGrid', linspace(0.1, 0.9, 9)', ...
    'Display', true);

%% 3. 可视化结果
figure('Name', 'QQR Parameter Surfaces', 'Position', [100, 100, 1200, 400]);

% 绘制斜率参数表面
qqr.plotSurface(result, 'Parameter', 'slope', 'Subplot', [1, 3, 1]);
title('Slope \beta_1(\theta, \tau)', 'Interpreter', 'tex');

% 绘制截距参数表面
qqr.plotSurface(result, 'Parameter', 'intercept', 'Subplot', [1, 3, 2]);
title('Intercept \beta_0(\theta, \tau)', 'Interpreter', 'tex');

% 绘制滞后系数表面
qqr.plotSurface(result, 'Parameter', 'alpha', 'Subplot', [1, 3, 3]);
title('Lag coefficient \alpha(\theta)', 'Interpreter', 'tex');

%% 4. QQ 验证（论文 Section 5, Eq. 10）
% 验证 QQ 估计的 τ 平均值是否与标准 QR 一致
report = qqr.validate(result, dataTbl, 'Display', true);

%% 5. 使用滞后自变量
result_lagged = qqr.estimate(dataTbl, ...
    'LagIndepVar', true, ...
    'ThetaGrid', linspace(0.1, 0.9, 9), ...
    'TauGrid', linspace(0.1, 0.9, 9)', ...
    'Display', true);

figure('Name', 'QQR with Lagged Independent Variable');
qqr.plotSurface(result_lagged, 'Parameter', 'slope');
title('Slope with lagged independent variable');

%% 6. 自定义参数约束（可选）
% 如果需要限制参数范围（例如保证平稳性）
result_constrained = qqr.estimate(dataTbl, ...
    'LagIndepVar', true, ...
    'ThetaGrid', linspace(0.1, 0.9, 9), ...
    'TauGrid', linspace(0.1, 0.9, 9)', ...
    'LowerBound', [-Inf, -Inf, -0.99], ...
    'UpperBound', [Inf, Inf, 0.99]);

%% 7. BDS 独立性检验
% 测试残差是否存在非线性依赖结构
% BDS 检验用于检测时间序列的非线性依赖性
residuals = y - mean(y);  % 简单示例：使用去均值后的 y 作为残差
[bdsstat, pval] = qqr.bds(residuals, 5);  % 最大嵌入维度为 5

disp('BDS Independence Test Results:');
disp(table((2:5)', bdsstat', pval', ...
    'VariableNames', {'Dimension', 'BDS_Stat', 'P_Value'}));

disp('=== 示例运行完成 ===');
