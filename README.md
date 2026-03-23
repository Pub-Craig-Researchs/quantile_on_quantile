# Quantile-on-Quantile Regression (QQR) Toolbox for MATLAB

A MATLAB implementation of the Quantile-on-Quantile Regression method for analyzing dependence structures across conditional distributions.

---

## Overview

Quantile-on-Quantile Regression (QQR) extends traditional quantile regression by estimating how the relationship between two variables varies across different quantile levels of both the dependent and independent variables. This approach captures asymmetric and tail-dependent relationships that standard methods may miss.

The method estimates coefficients β₀(θ,τ) and β₁(θ,τ) that vary with:
- **θ**: the quantile of the dependent variable (Y)
- **τ**: the quantile of the independent variable (X)

This toolbox implements the methodology from:

> Sim, N., & Zhou, H. (2015). Oil prices, US stock return, and the dependence between their quantiles. *Journal of Banking & Finance*, 55, 1-8. https://doi.org/10.1016/j.jbankfin.2015.01.013

---

## Requirements

- **MATLAB R2019b** or later (required for `arguments` block validation)
- **Optimization Toolbox** (for `fmincon`)
- **Statistics and Machine Learning Toolbox** (for `ksdensity`, `normpdf`)

---

## Installation

1. Clone or download this repository:
   ```bash
   git clone https://github.com/your-username/quantile_on_quantile.git
   ```

2. Ensure the project root directory is on your MATLAB path:
   ```matlab
   addpath('/path/to/quantile_on_quantile');
   ```


---

## Quick Start

```matlab
% Generate sample data
data = table(randn(500,1), randn(500,1), 'VariableNames', {'X','Y'});

% Estimate QQR model
result = qqr.estimate(data, 'ThetaGrid', linspace(0.05,0.95,19), ...
                            'TauGrid', linspace(0.05,0.95,19));

% Visualize slope surface
qqr.plotSurface(result);

% Validate against standard quantile regression
report = qqr.validate(result, data);
```

---

## API Reference

### qqr.loss

Computes the weighted check-function loss (Equation 6 of the paper).

```matlab
L = qqr.loss(beta, dataTbl, theta, tau, h, lagIndepVar)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `beta` | 3×1 vector | Coefficient vector [β₀; β₁; α] |
| `dataTbl` | table | Data table with independent variable(s) first, dependent variable last |
| `theta` | scalar | Quantile of dependent variable, θ ∈ (0,1) |
| `tau` | scalar | Quantile of independent variable, τ ∈ (0,1) |
| `h` | scalar | Bandwidth, h > 0 |
| `lagIndepVar` | logical | Whether to lag the independent variable |

**Returns:** Scalar loss value.

---

### qqr.estimate

Main estimation function for QQR coefficients.

```matlab
result = qqr.estimate(dataTbl)
result = qqr.estimate(dataTbl, Name, Value, ...)
```

**Parameters:**

| Name-Value Pair | Default | Description |
|-----------------|---------|-------------|
| `LagIndepVar` | `false` | Lag independent variable by one period |
| `ThetaGrid` | `linspace(0.05,0.95,19)` | Grid of θ quantiles for dependent variable |
| `TauGrid` | `linspace(0.05,0.95,19)` | Grid of τ quantiles for independent variable |
| `Bandwidth` | `'auto'` | Bandwidth value or `'auto'` for plug-in selection |
| `MinBandwidth` | `0.05` | Minimum bandwidth when using auto-selection |
| `LowerBound` | `[]` | Lower bounds for [β₀, β₁, α] optimization |
| `UpperBound` | `[]` | Upper bounds for [β₀, β₁, α] optimization |
| `InitialPoint` | `[0, 0, 0]` | Initial coefficient values for optimization |
| `Display` | `false` | Display progress during estimation |

**Returns:** Struct with fields:

| Field | Description |
|-------|-------------|
| `intercept` | nθ × nτ matrix of β₀(θ,τ) estimates |
| `slope` | nθ × nτ matrix of β₁(θ,τ) estimates |
| `alpha` | nθ × nτ matrix of α(θ,τ) estimates |
| `theta` | Vector of θ quantile grid points |
| `tau` | Vector of τ quantile grid points |
| `bandwidth` | Bandwidth used in estimation |
| `varNames` | Variable names from input table |
| `lagIndepVar` | Whether lagging was applied |
| `options` | Struct of all estimation options |

---

### qqr.plotSurface

Creates contour visualization of QQR coefficient surfaces.

```matlab
fig = qqr.plotSurface(result)
fig = qqr.plotSurface(result, Name, Value, ...)
```

**Parameters:**

| Name-Value Pair | Default | Description |
|-----------------|---------|-------------|
| `Parameter` | `'slope'` | Which parameter to plot: `'slope'`, `'intercept'`, or `'alpha'` |
| `Subplot` | `[1,1,1]` | Subplot specification [rows, cols, index] |
| `Colormap` | `'gray'` | Colormap name |
| `NumContours` | `20` | Number of contour levels |
| `ReverseColormap` | `true` | Reverse colormap direction |

**Returns:** Figure handle.

---

### qqr.validate

Validates QQR results by comparing τ-averaged coefficients with standard quantile regression estimates (Section 5, Equation 10 of the paper).

```matlab
report = qqr.validate(result, dataTbl)
report = qqr.validate(result, dataTbl, Name, Value, ...)
```

**Parameters:**

| Name-Value Pair | Default | Description |
|-----------------|---------|-------------|
| `Display` | `true` | Display validation summary |

**Returns:** Struct with fields:

| Field | Description |
|-------|-------------|
| `qqMeanIntercept` | τ-averaged QQR intercepts |
| `qqMeanSlope` | τ-averaged QQR slopes |
| `qqMeanAlpha` | τ-averaged QQR alphas |
| `qrIntercept` | Standard QR intercept estimates |
| `qrSlope` | Standard QR slope estimates |
| `qrAlpha` | Standard QR alpha estimates |
| `theta` | θ grid used |
| `interceptMAE` | Mean absolute error for intercepts |
| `slopeMAE` | Mean absolute error for slopes |
| `alphaMAE` | Mean absolute error for alphas |
| `interceptCorr` | Correlation for intercepts |
| `slopeCorr` | Correlation for slopes |
| `alphaCorr` | Correlation for alphas |

---

### qqr.bds

Brock-Dechert-Scheinkman (BDS) test for independence based on the correlation dimension.

```matlab
[W, SIG, C, C1, K] = qqr.bds(SERIES, MAXDIM, DISTANCE, FLAG, MAXRAM)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `SERIES` | vector | (required) | Time-series vector (typically residuals from a regression) |
| `MAXDIM` | scalar | `2` | Maximum embedding dimension (≥1) |
| `DISTANCE` | scalar | `1.5` | Dimensional distance parameter |
| `FLAG` | 0 or 1 | `0` | Method for determining epsilon: 0=fraction of std, 1=correlation integral target |
| `MAXRAM` | scalar | `150` | Maximum memory (MB) available for computation |

**Returns:**

| Output | Description |
|--------|-------------|
| `W` | BDS statistics for each dimension 2 to MAXDIM |
| `SIG` | Significance levels (two-tailed p-values) |
| `C` | Correlation integral estimates for dimensions 2 to MAXDIM |
| `C1` | First-order correlation integral estimates |
| `K` | Parameter estimate used in variance calculation |

**Example:**
```matlab
% Test residuals for remaining dependence
residuals = y - X * beta_hat;
[W, pval] = qqr.bds(residuals, 5);
fprintf('BDS statistic (dim=2): %.3f, p-value: %.4f\n', W(1), pval(1));
```

---

## Data Format

Input data must be a MATLAB `table` with:
- **Independent variable(s)** in the first column(s)
- **Dependent variable** in the last column

**Example:**
```matlab
% StockReturn is the dependent variable (Y)
% OilPrice is the independent variable (X)
data = table(OilPrice, StockReturn);

% Or with explicit naming:
data = table(oil, stock, 'VariableNames', {'OilPrice', 'StockReturn'});
```

The toolbox automatically uses variable names for plot labels and output.

---

## Examples

### Basic Estimation and Plotting

```matlab
% Load or create data
load financialData.mat  % Assume this contains oilPrice and stockReturn
data = table(oilPrice, stockReturn);

% Estimate with default settings
result = qqr.estimate(data);

% Plot slope surface
figure;
qqr.plotSurface(result, 'Parameter', 'slope');
title('Slope Coefficient Surface');
```

### Multi-Panel Visualization

```matlab
result = qqr.estimate(data);

figure('Position', [100 100 1200 400]);

qqr.plotSurface(result, 'Parameter', 'intercept', 'Subplot', [1,3,1]);
title('Intercept \beta_0(\theta,\tau)');

qqr.plotSurface(result, 'Parameter', 'slope', 'Subplot', [1,3,2]);
title('Slope \beta_1(\theta,\tau)');

qqr.plotSurface(result, 'Parameter', 'alpha', 'Subplot', [1,3,3]);
title('AR Coefficient \alpha(\theta,\tau)');
```

### Constrained Estimation

```matlab
% Impose stationarity constraint on alpha: |α| < 1
result = qqr.estimate(data, ...
    'LowerBound', [-Inf, -Inf, -0.99], ...
    'UpperBound', [Inf, Inf, 0.99], ...
    'Display', true);
```

### Validation Workflow

```matlab
result = qqr.estimate(data, 'ThetaGrid', linspace(0.1, 0.9, 17));
report = qqr.validate(result, data, 'Display', true);

% Check validation metrics
fprintf('Slope correlation: %.4f\n', report.slopeCorr);
fprintf('Slope MAE: %.4f\n', report.slopeMAE);
```

### Complete Workflow

See `examples/example_qqr.m` for a full demonstration including:
- Data preparation
- Estimation with custom grids
- Visualization of all coefficient surfaces
- Validation against standard quantile regression
- Interpretation guidelines

---

## Mathematical Background

### Check Function

The asymmetric check (loss) function for quantile θ:

```
ρ_θ(u) = u · (θ − I(u < 0))
```

where I(·) is the indicator function.

### QQR Model

The quantile-on-quantile regression model:

```
Q_θ(Y_t | X_t) = β₀(θ,τ) + β₁(θ,τ) · (X_t − Q_τ(X)) + α(θ) · Y_{t-1}
```

where:
- Q_θ(Y_t | X_t) is the θ-th conditional quantile of Y given X
- Q_τ(X) is the τ-th unconditional quantile of X
- β₀(θ,τ) captures the intercept surface
- β₁(θ,τ) captures the slope surface (dependence structure)
- α(θ) is the autoregressive coefficient

### Kernel Weighting

Coefficients are estimated using Gaussian kernel weights based on the empirical CDF:

```
K_h(F_n(X_t) − τ) = (1/h) · φ((F_n(X_t) − τ)/h)
```

where φ(·) is the standard normal density and h is the bandwidth.

### Bandwidth Selection

Default bandwidth is computed using the `ksdensity` plug-in method applied to the empirical CDF values, with a minimum floor of 0.05 to ensure numerical stability.

---

## Project Structure

```
quantile_on_quantile/
├── +qqr/                      # Main package (MATLAB namespace)
│   ├── bds.m                  # BDS independence test
│   ├── estimate.m             # QQR estimation
│   ├── loss.m                 # Check-function loss
│   ├── plotSurface.m          # Contour visualization
│   └── validate.m             # Validation against QR
├── examples/
│   └── example_qqr.m          # Full workflow example
└── README.md
```

---

## Citation

If you use this toolbox in your research, please cite the original paper:

```bibtex
@article{sim2015oil,
  title={Oil prices, {US} stock return, and the dependence between their quantiles},
  author={Sim, Nicholas and Zhou, Hongtao},
  journal={Journal of Banking \& Finance},
  volume={55},
  pages={1--8},
  year={2015},
  publisher={Elsevier},
  doi={10.1016/j.jbankfin.2015.01.013}
}
```

---

## License

This project is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html) (GPL-3.0).

See [LICENSE](LICENSE) for the full license text.
