%  example script 6 AOADMM for CMTF 
% In this example, we create a synthetic dataset consisting of one CP tensor
% of order 3 (modes 1,2,3) and size [50,60,40] with 3 components and 
% two matrices (modes 4,5 and modes 6,7) of size [50,70] and [60,80], both with 3 components. 
% The tensor and the first matrix are directly coupled in the first mode (coupling
% 1 for modes 1 and 4) and the tensor is coupled with the second matrix in
% teh second mode (coupling 2 for modes 2 and 6). Modes
% 1,2,4,5,6 and 7 are also generated to be non-negative.
% Some Gaussian noise is added to the synthetic datasets.
% The loss functions is set to be squared Frobenius norm.
% Non-negativity constraints are enforced for modes 1,2,4,5,6 and 7.

function tests = test_matrix_matrix_CP_nonneg
    tests = functiontests(localfunctions);
end

function test_workflow(testCase)

    % specify synthetic data
    sz     = {50,60,40,50,70,60,80}; %size of each mode
    P      = 3; %number of tensors
    lambdas_data= {[1 1 1], [1 1 1], [1 1 1]}; % norms of components in each data set (length of each array specifies the number of components in each dataset)
    modes  = {[1 2 3], [4 5],[6 7]}; % which modes belong to which dataset: every mode should have its unique number d, sz(d) corresponds to size of that mode
    noise = 0.2; %level of noise, for gaussian noise only!
    distr_data = {@(x,y) rand(x,y), @(x,y) rand(x,y),@(x,y) randn(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y)}; % function handle of distribution of data within each factor matrix /or Delta if linearly coupled, x,y are the size inputs %coupled modes need to have same distribution! If not, just the first one will be considered
    normalize_columns = 0; %wether or not to normalize columns of the created factor matrices, this might destroy the distribution!
    % specify tensor model
    model{1} = 'CP';
    model{2} = 'CP';
    model{3} = 'CP';
    % specify couplings
    coupling.lin_coupled_modes = [1 2 0 1 0 2 0]; % which modes are coupled, coupled modes get the same number (0: uncoupled)
    coupling.coupling_type = [0 0]; % for each coupling number in the array lin_coupled_modes, set the coupling type: 0 exact coupling, 1: HC=Delta, 2: CH=Delta, 3: C=HDelta, 4: C=DeltaH
    coupling.coupl_trafo_matrices = cell(7,1); % cell array with coupling transformation matrices for each mode (if any, otherwise keep empty)
    
    % set the fitting function for each dataset: 'Frobenius' for squared
    % Frobenius norm, 'KL' for KL divergence, IS for Itakura-Saito, 'beta' for other beta divergences (give beta in loss_function_param),...more todo
    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';
    loss_function{3} = 'Frobenius';
    loss_function_param{1} = [];
    loss_function_param{2} = [];
    loss_function_param{3} = [];
    % check model
    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
    
    % set initialization options
    init_options.lambdas_init = {[1 1 1], [1 1 1], [1 1 1]}; %norms of components in each data set for initialization
    init_options.nvecs = 0; % wether or not to use cmtf_nvecs.m funcion for initialization of factor matrices Ci (if true, distr_data and normalize are ignored for Ci, not for Zi)
    init_options.distr = distr_data; % distribution of the initial factor matrices and their auxiliary variables
    init_options.normalize = 1; % wether or not to normalize the columns of the initial factor matrices (might destroy the distribution)
    
    % set constraints
    constrained_modes = [1 0 0 1 1 1 1]; % 1 if the mode is constrained in some way, 0 otherwise, put the same for coupled modes!
    
    constraints = cell(length(constrained_modes),1); % cell array of length number of modes containing the type of constraint or regularization for each mode, empty if no constraint
    %specify constraints-regularizations for each mode, find the options in the file "List of constraints and regularizations.txt"
    constraints{1} = {'non-negativity'};
    constraints{2} = {'non-negativity'};
    constraints{4} = {'non-negativity'};
    constraints{5} = {'non-negativity'};
    constraints{6} = {'non-negativity'};
    constraints{7} = {'non-negativity'};

    % set weights
    weights = [1/3 1/3 1/3]; %weight w_i for each data set
    
    % build model
    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;
    
    % create data
    [X, Atrue, Deltatrue,sigmatrue] = cmtf.utils.create_coupled_data('model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, 'noise', noise,'coupling',coupling,'normalize_columns',normalize_columns,'distr_data',distr_data,'loss_function',Z.loss_function); %create data
    % create Z.object and normalize
    normZ=cell(P,1);
    for p=1:P
        Z.object{p} = X{p};
        if strcmp(model{p},'CP')
            normZ{p} = norm(Z.object{p});
            Z.object{p} = Z.object{p}/normZ{p};
        elseif strcmp(model{p},'PAR2')
            normZ{p} = 0;
            for k=1:length(Z.object{p})
                normZ{p} = normZ{p} + norm(Z.object{p}{k},'fro')^2;
            end
            normZ{p} = sqrt(normZ{p});
            for k=1:length(Z.object{p})
                Z.object{p}{k} = Z.object{p}{k}/normZ{p};
            end
        end
    end
    
    % Create random initialization
    init_fac = cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);
    
    % set options 
    options.Display ='no'; %  set to 'iter' or 'final' or 'no'
    options.DisplayIters = 10;
    options.MaxOuterIters = 4000;
    options.MaxInnerIters = 5;
    options.AbsFuncTol   = 1e-4;
    options.OuterRelTol = 1e-8;
    options.innerRelPrTol_coupl = 1e-3;
    options.innerRelPrTol_constr = 1e-3;
    options.innerRelDualTol_coupl = 1e-3;
    options.innerRelDualTol_constr = 1e-3;
    options.bsum = 0; % wether or not to use AO with BSUM regularization
    %options.bsum_weight = 1e-3; %set the penalty parameter (mu) for BSUM regularization
    options.eps_log = 1e-10; % for KL divergence log(x+eps) for numerical stability
    %options.lbfgsb_options = lbfgsb_options;
    
    % run algorithm
    [Zhat,Fac,FacInit,out] = cmtf.aoadmm.cmtf_AOADMM(Z,'alg_options',options,'init',init_fac,'init_options',init_options); 

    % FIT
    Fit1 = 100*(1-norm(Z.object{1}-full(Zhat{1}))^2/norm(Z.object{1})^2);
    Fit2 = 100*(1-norm(Z.object{2}-full(Zhat{2}))^2/norm(Z.object{2})^2);
    Fit3 = 100*(1-norm(Z.object{3}-full(Zhat{3}))^2/norm(Z.object{3})^2);

    % FMS 
    true_ktensor{1} =(ktensor(lambdas_data{1}'./normZ{1},Atrue(modes{1})));
    FMS1 = score(Zhat{1},true_ktensor{1});
    true_ktensor{2} =(ktensor(lambdas_data{2}'./normZ{2},Atrue(modes{2})));
    FMS2 = score(Zhat{2},true_ktensor{2});
    true_ktensor{3} =(ktensor(lambdas_data{3}'./normZ{3},Atrue(modes{3})));
    FMS3 = score(Zhat{3},true_ktensor{3});    
    
    % Report output for logs
    Fit1
    Fit2
    Fit3
    FMS1
    FMS2
    FMS3

    % Test expected output (see ~/examples/expectedOutput)
    testCase.verifyGreaterThanOrEqual(Fit1, 95);
    testCase.verifyGreaterThanOrEqual(FMS1, 0.99);

    testCase.verifyGreaterThanOrEqual(Fit2, 95);
    testCase.verifyGreaterThanOrEqual(FMS2, 0.95);

    testCase.verifyGreaterThanOrEqual(Fit3, 95);
    testCase.verifyGreaterThanOrEqual(FMS3, 0.95);
end

