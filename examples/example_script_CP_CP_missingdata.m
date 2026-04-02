%% Example script: coupled CP-CP with missing data in one block
%
% Two rank-3 CP tensors share an exact coupling in their first mode.
%
%   Tensor 1 (complete):  size [20 x 30 x 40], modes [1 2 3]
%   Tensor 2 (incomplete): size [20 x 50 x 60], modes [4 5 6]
%                          ~20% of entries missing at random
%
% Both tensors are fitted with squared Frobenius loss. The shared factor
% (mode 1 / mode 4) is recovered jointly; missing entries in tensor 2 are
% imputed at each outer iteration via the EM loop.
%
% The example shows:
%   - How to build Z.miss for the incomplete block
%   - Convergence of the tensor fit alongside the EM imputation metric
%   - How to evaluate fit on observed entries and imputation quality

close all
clear all
run ../setup.m
rng(123, 'twister');

%% Synthetic data specification
sz     = {20, 30, 40, 20, 50, 60};   % sizes of all 6 modes
P      = 2;                           % number of tensors
R      = 3;                           % number of components (same in both tensors)
modes  = {[1 2 3], [4 5 6]};         % modes belonging to each tensor
lambdas_data = {ones(1,R), ones(1,R)};
noise  = 0.05;                        % small Gaussian noise
normalize_columns = 0;
distr_data = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
              @(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};

%% Tensor model
model{1} = 'CP';
model{2} = 'CP';

%% Coupling: tensors share mode 1 and mode 4 exactly (coupling type 0)
coupling.lin_coupled_modes   = [1 0 0 1 0 0];
coupling.coupling_type       = [0];
coupling.coupl_trafo_matrices = cell(6, 1);

%% Loss functions
loss_function{1}      = 'Frobenius';
loss_function{2}      = 'Frobenius';
loss_function_param{1} = [];
loss_function_param{2} = [];
cmtf.utils.check_data_input(sz, modes, lambdas_data, coupling, loss_function, model);

%% No constraints on any mode (focusing on the missing-data feature)
constrained_modes = [0 0 0 0 0 0];
constraints       = cell(6, 1);

%% Equal weights
weights = [1/2, 1/2];

%% Assemble Z
Z.loss_function       = loss_function;
Z.loss_function_param = loss_function_param;
Z.model               = model;
Z.modes               = modes;
Z.size                = sz;
Z.coupling            = coupling;
Z.constrained_modes   = constrained_modes;
Z.constraints         = constraints;
Z.weights             = weights;

%% Create synthetic ground-truth data
[X, Atrue] = cmtf.utils.create_coupled_data('model', model, 'size', sz, ...
    'modes', modes, 'lambdas', lambdas_data, 'noise', noise, ...
    'coupling', coupling, 'normalize_columns', normalize_columns, ...
    'distr_data', distr_data, 'loss_function', Z.loss_function);

%% Normalize and store data in Z
normZ = cell(P, 1);
for p = 1:P
    Z.object{p} = X{p};
    normZ{p}    = norm(Z.object{p});
    Z.object{p} = Z.object{p} / normZ{p};
end

%% Build missing-data mask for tensor 2 only (~20% missing)
sz2       = [sz{4}, sz{5}, sz{6}];    % [20 50 60]
n_entries = prod(sz2);
miss_frac = 0.20;
miss_mask = true(sz2);
miss_mask(randperm(n_entries, round(miss_frac * n_entries))) = false;

% Store the ground-truth values at missing positions (for later evaluation)
X2_true_norm = double(Z.object{2});   % noiseless model not available directly;
                                       % we use the noisy but observed data as proxy

% Provide the mask only for tensor 2; tensor 1 is fully observed (no Z.miss{1})
Z.miss{1} = [];                        % no missing data in tensor 1
Z.miss{2} = miss_mask;                 % logical mask: true = observed

fprintf('Tensor 2: %d total entries, %d observed (%.0f%% missing)\n', ...
    n_entries, sum(miss_mask(:)), 100 * miss_frac);

%% Initialization
init_options.lambdas_init = lambdas_data;
init_options.nvecs        = 0;
init_options.distr        = distr_data;
init_options.normalize    = 1;
init_fac = cmtf.utils.init_coupled_AOADMM_CMTF(Z, 'init_options', init_options);

%% Solver options
options.Display               = 'iter';
options.DisplayIters          = 25;
options.MaxOuterIters         = 2000;
options.MaxInnerIters         = 5;
options.AbsFuncTol            = 1e-7;
options.OuterRelTol           = 1e-8;
options.innerRelPrTol_coupl   = 1e-4;
options.innerRelDualTol_coupl = 1e-4;
options.innerRelPrTol_constr  = 1e-4;
options.innerRelDualTol_constr = 1e-4;
options.bsum                  = 0;
options.eps_log               = 1e-10;

%% Run AO-ADMM with EM imputation
fprintf('\nRunning AO-ADMM with EM missing-data imputation...\n');
tic
[Zhat, Fac, ~, out] = cmtf.aoadmm.cmtf_AOADMM(Z, 'alg_options', options, ...
    'init', init_fac, 'init_options', init_options);
t_elapsed = toc;
fprintf('Done in %.1f s (%d outer iterations)\n', t_elapsed, out.OuterIterations);
fprintf('EM convergence (f_rel_missing_final): %.2e\n', out.f_rel_missing_final);

%% Fit on tensor 1 (fully observed)
Fit1 = 100 * (1 - norm(Z.object{1} - full(Zhat{1}))^2 / norm(Z.object{1})^2);

%% Fit on tensor 2 — observed entries only
X2_rec  = double(full(Zhat{2}));
X2_data = double(Z.object{2});   % contains imputed values at missing positions after EM
obs_res = X2_data(miss_mask) - X2_rec(miss_mask);
Fit2_obs = 100 * (1 - norm(obs_res)^2 / norm(X2_data(miss_mask))^2);

%% Factor match scores
true_ktensor{1} = ktensor(lambdas_data{1}' ./ normZ{1}, Atrue(modes{1}));
true_ktensor{2} = ktensor(lambdas_data{2}' ./ normZ{2}, Atrue(modes{2}));
FMS1 = score(Zhat{1}, true_ktensor{1}, 'lambda_penalty', false);
FMS2 = score(Zhat{2}, true_ktensor{2}, 'lambda_penalty', false);

fprintf('\nResults:\n');
fprintf('  Tensor 1  — Fit: %.2f%%   FMS: %.4f\n', Fit1, FMS1);
fprintf('  Tensor 2  — Fit (observed): %.2f%%   FMS: %.4f\n', Fit2_obs, FMS2);

%% Imputation quality: compare recovered values at missing positions to
%  the noiseless ground-truth model at those positions
M2_true = double(full(true_ktensor{2}));  % noiseless rank-R model
imputed_vals     = X2_rec(~miss_mask);
groundtruth_vals = M2_true(~miss_mask);
err_imputation   = norm(imputed_vals - groundtruth_vals) / norm(groundtruth_vals);
fprintf('  Imputation error (vs noiseless model): %.4f\n', err_imputation);

%% Plots
iters = 0:out.OuterIterations;

figure('Name', 'Convergence — coupled CP-CP with missing data', 'NumberTitle', 'off');

subplot(2, 2, 1)
semilogy(iters, out.func_val_conv, 'b-',  'LineWidth', 1.5)
hold on
semilogy(iters, out.func_coupl_conv, 'r--', 'LineWidth', 1.5)
xlabel('Outer iteration')
ylabel('Value')
legend('Tensor fit', 'Coupling residual', 'Location', 'southwest')
title('Objective convergence')
grid on

subplot(2, 2, 2)
semilogy(out.time_at_it, out.func_val_conv, 'b-', 'LineWidth', 1.5)
hold on
semilogy(out.time_at_it, out.func_coupl_conv, 'r--', 'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('Value')
legend('Tensor fit', 'Coupling residual', 'Location', 'southwest')
title('Convergence vs time')
grid on

subplot(2, 2, 3)
% Compare true and recovered factor columns for the shared mode
A_shared_true = Atrue{1};               % true shared factor (mode 1 = mode 4)
A_shared_rec  = Fac.fac{1};            % recovered (in normalised space)
[~, perm]     = max(abs(A_shared_true' * A_shared_rec), [], 2);
plot(A_shared_true(:, 1:R), '-',  'LineWidth', 1.5)
hold on
plot(A_shared_rec(:, perm), '--', 'LineWidth', 1.5)
xlabel('Index')
ylabel('Factor loading')
title(sprintf('Shared factor (FMS=%.3f)', FMS1))
legend('True 1','True 2','True 3','Est 1','Est 2','Est 3','Location','best')
grid on

subplot(2, 2, 4);
% Compare hidden missing entries to imputed entries
miss_mask = ~Z.miss{2};
tmp = double(Z.object{2});
old_vals = tmp(miss_mask);
Xhat = double(Zhat{2});
new_vals = Xhat(miss_mask);
scatter(old_vals, new_vals);
xlabel("Hidden values");
ylabel("Imputed values");
title("Imputation of hidden entries")
grid on
sgtitle(sprintf('Coupled CP-CP, tensor 2 with %.0f%% missing (EM imputation)', 100*miss_frac))
