%  example script 4 AOADMM for CMTF 
% In this example, we create a synthetic ragged PARAFAC2 tensor
% of order 3 (modes 1,2,3) and size [40,[61:1:120],60] with 3 components  We add
% Gaussian noise with level 0.2 and use Frobenius norm loss.
%%

function tests = test_irregularPAR2
    tests = functiontests(localfunctions);
end

function test_workflow(testCase)

    % Control randomness
    rng(123, 'twister');
    
    % specify synthetic data
    sz     = {40,[61:1:120],60}; %size of each mode
    P      = 1; %number of tensors
    lambdas_data= {[1 1 1]}; % norms of components in each data set (length of each array specifies the number of components in each dataset)
    modes  = {[1 2 3]}; % which modes belong to which dataset: every mode should have its unique number d, sz(d) corresponds to size of that mode
    noise = 0.2; %level of noise, for gaussian noise only!
    distr_data = {@(x,y) randn(x,y), @(x,y) randn(x,y),@(x,y) rand(x,y)+0.1}; % function handle of distribution of data within each factor matrix /or Delta if linearly coupled, x,y are the size inputs %coupled modes need to have same distribution! If not, just the first one will be considered
    normalize_columns = 0; %wether or not to normalize columns of the created factor matrices, this might destroy the distribution!
    % specify tensor model
    model{1} = 'PAR2';
    % specify couplings
    coupling.lin_coupled_modes = [0 0 0]; % which modes are coupled, coupled modes get the same number (0: uncoupled)
    coupling.coupling_type = []; % for each coupling number in the array lin_coupled_modes, set the coupling type: 0 exact coupling, 1: HC=Delta, 2: CH=Delta, 3: C=HDelta, 4: C=DeltaH
    coupling.coupl_trafo_matrices = cell(3,1); % cell array with coupling transformation matrices for each mode (if any, otherwise keep empty)
    
    
    % set the fitting function for each dataset: 'Frobenius' for squared
    % Frobenius norm, 'KL' for KL divergence, IS for Itakura-Saito, 'beta' for other beta divergences (give beta in loss_function_param),...more todo
    loss_function{1} = 'Frobenius';
    loss_function_param{1} = [];
    % check model
    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
    
    % set initialization options
    init_options.lambdas_init = {[1 1 1]}; %norms of components in each data set for initialization
    init_options.nvecs = 0; % wether or not to use cmtf_nvecs.m funcion for initialization of factor matrices Ci (if true, distr_data and normalize are ignored for Ci, not for Zi)
    init_options.distr =distr_data; % distribution of the initial factor matrices and their auxiliary variables
    init_options.normalize = 1; % wether or not to normalize the columns of the initial factor matrices (might destroy the distribution)
    
    % set constraints
    constrained_modes = [0 0 1]; % 1 if the mode is constrained in some way, 0 otherwise, put the same for coupled modes!
    
    constraints = cell(length(constrained_modes),1); % cell array of length number of modes containing the type of constraint or regularization for each mode, empty if no constraint
    %specify constraints-regularizations for each mode, find the options in the file "List of constraints and regularizations.txt"
    constraints{3} = {'non-negativity'};
    
    % set weights
    weights = [1]; %weight w_i for each data set
    
    
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
    [X, Atrue, Deltatrue,sigmatrue] = cmtf.utils.create_irregularPARAFAC2_coupled_data('model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, 'noise', noise,'coupling',coupling,'normalize_columns',normalize_columns,'distr_data',distr_data,'loss_function',Z.loss_function); %create data
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
    options.AbsFuncTol   = 1e-7;
    options.OuterRelTol = 1e-8;
    options.innerRelPrTol_coupl = 1e-5;
    options.innerRelPrTol_constr = 1e-5;
    options.innerRelDualTol_coupl = 1e-5;
    options.innerRelDualTol_constr = 1e-5;
    options.bsum = 0; % wether or not to use AO with BSUM regularization
    %options.bsum_weight = 1e-3; %set the penalty parameter (mu) for BSUM regularization
    options.eps_log = 1e-10; % for KL divergence log(x+eps) for numerical stability
    %options.lbfgsb_options = lbfgsb_options;
    
    % run algorithm
    [Zhat,Fac,FacInit,out] = cmtf.aoadmm.cmtf_AOADMM(Z,'alg_options',options,'init',init_fac,'init_options',init_options); 
      
    %% FIT
    Fit1 = 0;
    Fitx = 0;
    for k=1:length(sz{2})
        Fit1 = Fit1 + norm(Z.object{1}{k}-Zhat{1}.A*diag(Zhat{1}.C(k,:))*Zhat{1}.Bk{k}','fro')^2;
        Fitx    = Fitx    + norm(Z.object{1}{k},'fro')^2;
    end
    Fit1 = 100*(1-Fit1/Fitx);

    % FMS 
    FMS_A = score(ktensor(ones(3,1),Zhat{1}.A),ktensor(ones(3,1),Atrue{1}),'lambda_penalty',false);
    FMS_C = score(ktensor(ones(3,1),Zhat{1}.C),ktensor(ones(3,1),Atrue{3}),'lambda_penalty',false);
    SollargeB = [];
    largeB = [];
    for k=1:length(sz{2})
        SollargeB = [SollargeB;Zhat{1}.Bk{k}];
        largeB = [largeB;Atrue{2}{k}];
    end
    FMS_B = score(ktensor(ones(3,1),SollargeB),ktensor(ones(3,1),largeB),'lambda_penalty',false);

    % Test expected output (see ~/examples/expectedOutput)
    testCase.verifyEqual(Fit1, 96.442, "AbsTol", 1e-2); % larger tolerance due to flakiness between platforms
    testCase.verifyEqual(FMS_A, 1, "AbsTol", 1e-3);
    testCase.verifyEqual(FMS_B, 0.9955, "AbsTol", 1e-3);
    testCase.verifyEqual(FMS_C, 1, "AbsTol", 1e-3);
end
    

