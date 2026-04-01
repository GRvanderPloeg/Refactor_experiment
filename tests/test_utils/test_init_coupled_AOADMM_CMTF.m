function tests = test_init_coupled_AOADMM_CMTF
    tests = functiontests(localfunctions);
end

function test_matrix(testCase)
    rng(0,'twister');

    model{1} = 'CP';
    sz     = {60,50};
    modes  = {[1 2]};
    lambdas_data= {[1 1 1]};

    coupling.lin_coupled_modes = [0 0];
    coupling.coupling_type = [];
    coupling.coupl_trafo_matrices = cell(2,1);

    distr_data = {@(x,y) randn(x,y), @(x,y) randn(x,y)};

    loss_function{1} = 'Frobenius';
    loss_function_param{1} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;
    Z.object{1} = tensor(rand(60, 50));

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 1;

    cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);
end

function test_CP(testCase)
    rng(0,'twister');

    model{1} = 'CP';
    sz     = {60,50,70};
    modes  = {[1 2 3]};
    lambdas_data= {[1 1 1]};
    
    coupling.lin_coupled_modes = [0 0 0];
    coupling.coupling_type = [];
    coupling.coupl_trafo_matrices = cell(3,1);
    
    distr_data = {@(x,y) randn(x,y), @(x,y) randn(x,y),@(x,y) randn(x,y)};

    loss_function{1} = 'Frobenius';    
    loss_function_param{1} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 1;
    Z.object{1} = tensor(rand(60, 50, 70));

    cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);
end

function test_regular_PAR2(testCase)
    rng(0,'twister');

    P = 1;
    model{1} = 'PAR2';
    sz     = {40,120*ones(1,60),60};
    modes  = {[1 2 3]};
    lambdas_data= {[1 1 1]};
    noise = 0.2;

    coupling.lin_coupled_modes = [0 0 0];
    coupling.coupling_type = [];
    coupling.coupl_trafo_matrices = cell(3,1);
    
    normalize_columns = 0;
    distr_data = {@(x,y) rand(x,y), @(x,y) rand(x,y),@(x,y) rand(x,y)+0.1};

    loss_function{1} = 'Frobenius';
    loss_function_param{1} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 1;
    
    [X, ~, ~,~] = cmtf.utils.create_coupled_data('model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, 'noise', noise,'coupling',coupling,'normalize_columns',normalize_columns,'distr_data',distr_data,'loss_function',Z.loss_function); %create data

    for p=1:P
        Z.object{p} = X{p};
    end

    cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);
end

function test_irregular_PAR2(testCase)
    rng(0,'twister');

    P = 1;
    model{1} = 'PAR2';
    sz     = {40,[61:1:120],60};
    modes  = {[1 2 3]};
    lambdas_data= {[1 1 1]};
    noise = 0.2;

    coupling.lin_coupled_modes = [0 0 0];
    coupling.coupling_type = [];
    coupling.coupl_trafo_matrices = cell(3,1);
    
    normalize_columns = 0;
    distr_data = {@(x,y) randn(x,y), @(x,y) randn(x,y),@(x,y) rand(x,y)+0.1};

    loss_function{1} = 'Frobenius';
    loss_function_param{1} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 1;

    [X, ~, ~,~] = cmtf.utils.create_coupled_data('model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, 'noise', noise,'coupling',coupling,'normalize_columns',normalize_columns,'distr_data',distr_data,'loss_function',Z.loss_function); %create data

    for p=1:P
        Z.object{p} = X{p};
    end

    cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);
end

function test_matrix_CP(testCase)
    rng(0,'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {20,30, 20,30,20};
    modes  = {[1 2], [3 4 5]};
    lambdas_data= {[1 1 1], [1 1 1]};

    coupling.lin_coupled_modes = [1 0 1 0 0];
    coupling.coupling_type = 0;
    coupling.coupl_trafo_matrices = cell(5,1);

    distr_data = {@(x,y) rand(x,y),@(x,y) randn(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y)};

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';
    loss_function_param{1} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 1;
    Z.object{1} = tensor(rand(20,30));
    Z.object{2} = tensor(rand(20,30,20));

    cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);
end

function test_CP_CP(testCase)
    rng(0,'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {20,30,40, 20,30,20};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data= {[1 1 1], [1 1 1]};

    coupling.lin_coupled_modes = [1 0 0 1 0 0];
    coupling.coupling_type = 0;
    coupling.coupl_trafo_matrices = cell(6,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y),@(x,y) randn(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y)};

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';
    loss_function_param{1} = [];
    loss_function_param{2} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 1;
    Z.object{1} = tensor(rand(20,30,40));
    Z.object{2} = tensor(rand(20,30,20));

    cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);
end

function test_matrix_PAR2(testCase)
    rng(0,'twister');

    P = 2;
    model{1} = 'CP';
    model{2} = 'PAR2';
    sz     = {40,60, 40,120*ones(1,50),50};
    modes  = {[1 2], [3 4 5]};
    lambdas_data= {[1 1 1], [1 1 1]};
    noise = 0.5;

    coupling.lin_coupled_modes = [1 0 1 0 0];
    coupling.coupling_type = 0;
    coupling.coupl_trafo_matrices = cell(5,1);
    
    normalize_columns = 0;
    distr_data = {@(x,y) rand(x,y), @(x,y) rand(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y)+0.1}; % function handle of distribution of data within each factor matrix /or Delta if linearly coupled, x,y are the size inputs %coupled modes need to have same distribution! If not, just the first one will be considered

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';
    loss_function_param{1} = [];
    loss_function_param{2} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 1;
    
    [X, ~, ~,~] = cmtf.utils.create_coupled_data('model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, 'noise', noise,'coupling',coupling,'normalize_columns',normalize_columns,'distr_data',distr_data,'loss_function',Z.loss_function); %create data

    for p=1:P
        Z.object{p} = X{p};
    end

    cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);
end

function test_CP_PAR2(testCase)
    rng(0,'twister');

    P = 2;
    model{1} = 'CP';
    model{2} = 'PAR2';
    sz     = {20,30,40, 20,30*ones(1,20),20};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data= {[1 1 1], [1 1 1]};
    noise = 0;

    coupling.lin_coupled_modes = [1 0 0 1 0 0];
    coupling.coupling_type = 0;
    coupling.coupl_trafo_matrices = cell(6,1);
    
    normalize_columns = 0;
    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y),@(x,y) randn(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y)+0.1};

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';
    loss_function_param{1} = [];
    loss_function_param{2} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 1;
    
    [X, ~, ~,~] = cmtf.utils.create_coupled_data('model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, 'noise', noise,'coupling',coupling,'normalize_columns',normalize_columns,'distr_data',distr_data,'loss_function',Z.loss_function); %create data

    for p=1:P
        Z.object{p} = X{p};
    end

    cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);
end

function test_matrix_CP_partialcoupling(testCase)
    rng(0,'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {50,30,40,50,70};
    modes  = {[1 2 3], [4 5]};
    lambdas_data= {[1 1 1 1], [1 1 1]};

    coupling.lin_coupled_modes = [1 0 0 1 0];
    coupling.coupling_type = 4;
    coupling.coupl_trafo_matrices = cell(5,1);
    coupling.coupl_trafo_matrices{1} = eye(4,4);
    coupling.coupl_trafo_matrices{4} = [eye(3,3);0 0 0];

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y),@(x,y) randn(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y)};

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';
    loss_function_param{1} = [];
    loss_function_param{2} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 1;
    Z.object{1} = tensor(rand(50,30,40));
    Z.object{1} = tensor(rand(50,70));

    cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);
end

function test_CP_CP_doublesamplingrate(testCase)
    rng(0,'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {50,30,40,100,70,80};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data= {[1 1 1 1], [1 1 1 1]};

    coupling.lin_coupled_modes = [1 0 0 1 0 0];
    coupling.coupling_type = 1;
    coupling.coupl_trafo_matrices = cell(6,1);
    coupling.coupl_trafo_matrices{1} = eye(50,50);
    coupling.coupl_trafo_matrices{4} = zeros(50,100);
    for i=1:50
      coupling.coupl_trafo_matrices{4}(i,i+(i-1)) = 1; %take every second entry
    end

    distr_data = {@(x,y) randn(x,y),@(x,y) randn(x,y), @(x,y) randn(x,y),@(x,y) randn(x,y),@(x,y) randn(x,y),@(x,y) rand(x,y)};

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';
    loss_function_param{1} = [];
    loss_function_param{2} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 1;
    Z.object{1} = tensor(rand(50,30,40));
    Z.object{2} = tensor(rand(100,70,80));

    cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);
end

function test_CP_CP_coupling_type2(testCase)
% Delta initialization for coupling type 2 (component-space coupling).
    rng(0,'twister');

    R = 3;
    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {20,30,40, 20,70,80};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,R), ones(1,R)};

    coupling.lin_coupled_modes = [1 0 0 1 0 0];
    coupling.coupling_type = 2;
    coupling.coupl_trafo_matrices = cell(6,1);
    d = 2;
    coupling.coupl_trafo_matrices{1} = rand(R, d);
    coupling.coupl_trafo_matrices{4} = rand(R, d);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';
    loss_function_param{1} = [];
    loss_function_param{2} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 0;
    Z.object{1} = tensor(rand(20,30,40));
    Z.object{2} = tensor(rand(20,70,80));

    A = cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);

    % coupling_fac{1} should have sz{1} rows and d columns
    verifyEqual(testCase, size(A.coupling_fac{1}), [sz{1}, d]);
end

function test_CP_CP_coupling_type3(testCase)
% Delta initialization for coupling type 3 (latent-space coupling).
    rng(0,'twister');

    R = 3;
    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {20,30,40, 25,70,80};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,R), ones(1,R)};

    coupling.lin_coupled_modes = [1 0 0 1 0 0];
    coupling.coupling_type = 3;
    coupling.coupl_trafo_matrices = cell(6,1);
    d = 5;
    coupling.coupl_trafo_matrices{1} = rand(sz{1}, d);
    coupling.coupl_trafo_matrices{4} = rand(sz{4}, d);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';
    loss_function_param{1} = [];
    loss_function_param{2} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 0;
    Z.object{1} = tensor(rand(20,30,40));
    Z.object{2} = tensor(rand(25,70,80));

    A = cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);

    % coupling_fac{1} should have d rows and R columns
    verifyEqual(testCase, size(A.coupling_fac{1}), [d, R]);
end

function test_matrix_matrix_CP(testCase)
    rng(0,'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    model{3} = 'CP';
    sz     = {50,60,40,50,70,60,80};
    modes  = {[1 2 3], [4 5],[6 7]};
    lambdas_data= {[1 1 1], [1 1 1], [1 1 1]};

    coupling.lin_coupled_modes = [1 2 0 1 0 2 0];
    coupling.coupling_type = [0 0]; 
    coupling.coupl_trafo_matrices = cell(7,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) rand(x,y),@(x,y) randn(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y)};

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';
    loss_function{3} = 'Frobenius';
    loss_function_param{1} = [];
    loss_function_param{2} = [];
    loss_function_param{3} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 1;
    Z.object{1} = tensor(rand(50,60,40));
    Z.object{2} = tensor(rand(50,70));
    Z.object{3} = tensor(rand(60,80));

    cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);
end

function test_CP_nvecs(testCase)
% nvecs=1: use SVD-based initialization for a CP tensor.
    rng(0,'twister');

    model{1} = 'CP';
    sz     = {20,30,40};
    modes  = {[1 2 3]};
    lambdas_data = {[1 1 1]};

    coupling.lin_coupled_modes = [0 0 0];
    coupling.coupling_type = [];
    coupling.coupl_trafo_matrices = cell(3,1);

    distr_data = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};

    loss_function{1} = 'Frobenius';
    loss_function_param{1} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;
    Z.object{1} = tensor(rand(20,30,40));

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 1;
    init_options.distr = distr_data;
    init_options.normalize = 0;

    A = cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);

    verifyEqual(testCase, size(A.fac{1}), [sz{1}, length(lambdas_data{1})]);
end

function test_PAR2_nvecs_regular(testCase)
% nvecs=1: use SVD-based initialization for regular PAR2. Covers PAR2 SVD
% branches including A.DeltaB, A.P, A.mu_DeltaB initialization.
    rng(0,'twister');

    model{1} = 'PAR2';
    sz     = {40, 120*ones(1,10), 10};
    modes  = {[1 2 3]};
    lambdas_data = {[1 1 1]};
    noise = 0.2;

    coupling.lin_coupled_modes = [0 0 0];
    coupling.coupling_type = [];
    coupling.coupl_trafo_matrices = cell(3,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) rand(x,y), @(x,y) rand(x,y)+0.1};

    loss_function{1} = 'Frobenius';
    loss_function_param{1} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    normalize_columns = 0;
    [X, ~, ~, ~] = cmtf.utils.create_coupled_data( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'normalize_columns', normalize_columns, ...
        'distr_data', distr_data, 'loss_function', loss_function);
    Z.object{1} = X{1};

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 1;
    init_options.distr = distr_data;
    init_options.normalize = 0;

    A = cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);

    verifyTrue(testCase, iscell(A.fac{2}));
    verifyEqual(testCase, length(A.fac{2}), length(sz{2}));
    verifyTrue(testCase, isfield(A, 'DeltaB'));
    verifyTrue(testCase, isfield(A, 'P'));
    verifyTrue(testCase, isfield(A, 'mu_DeltaB'));
end

function test_constrained_PAR2_mode2(testCase)
% Constrained PAR2 varying mode (mode 2) with non-negativity.
% Covers the constrained PAR2 branch (lines 102-112 of init function).
    rng(0,'twister');

    model{1} = 'PAR2';
    sz     = {40, 120*ones(1,10), 10};
    modes  = {[1 2 3]};
    lambdas_data = {[1 1 1]};
    noise = 0.2;

    coupling.lin_coupled_modes = [0 0 0];
    coupling.coupling_type = [];
    coupling.coupl_trafo_matrices = cell(3,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) rand(x,y), @(x,y) rand(x,y)+0.1};

    loss_function{1} = 'Frobenius';
    loss_function_param{1} = [];

    % Constrain mode 2 (PAR2 varying mode) with non-negativity
    constrained_modes = [0 1 0];
    constraints = cell(length(constrained_modes),1);
    constraints{2} = {'non-negativity'};
    weights = [1];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    normalize_columns = 0;
    [X, ~, ~, ~] = cmtf.utils.create_coupled_data( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'normalize_columns', normalize_columns, ...
        'distr_data', distr_data, 'loss_function', loss_function);
    Z.object{1} = X{1};

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 0;

    A = cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);

    % constraint_fac for the varying mode should be a cell of non-negative matrices
    verifyTrue(testCase, iscell(A.constraint_fac{2}));
    for k = 1:length(sz{2})
        verifyGreaterThanOrEqual(testCase, A.constraint_fac{2}{k}, 0);
    end
end

function test_matrix_CP_Poisson(testCase)
    rng(0,'twister');

    P = 2;
    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {50,30,40,50,70};
    modes  = {[1 2 3], [4 5]};
    lambdas_data= {[1 1 1], [1 1 1]};
    noise = 0.2;
    shape = 1;
    scale = 1;

    coupling.lin_coupled_modes = [1 0 0 1 0];
    coupling.coupling_type = 0;
    coupling.coupl_trafo_matrices = cell(5,1);
    
    normalize_columns = 0;
    distr_data = {@(x,y) gamrnd(shape,scale,x,y),@(x,y) gamrnd(shape,scale,x,y),@(x,y) gamrnd(shape,scale,x,y),@(x,y) gamrnd(shape,scale,x,y),@(x,y) gamrnd(shape,scale,x,y)};

    loss_function{1} = 'KL';
    loss_function{2} = 'KL';
    loss_function_param{1} = [];
    loss_function_param{2} = [];

    constrained_modes = [0 0 0];
    constraints = cell(length(constrained_modes),1);
    weights = [1/2 1/2];

    Z.loss_function = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model = model;
    Z.modes = modes;
    Z.size  = sz;
    Z.coupling = coupling;
    Z.constrained_modes = constrained_modes;
    Z.constraints = constraints;
    Z.weights = weights;

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs = 0;
    init_options.distr = distr_data;
    init_options.normalize = 1;
    
    [X, ~, ~,~] = cmtf.utils.create_coupled_data('model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, 'noise', noise,'coupling',coupling,'normalize_columns',normalize_columns,'distr_data',distr_data,'loss_function',Z.loss_function); %create data

    for p=1:P
        Z.object{p} = X{p};
    end

    cmtf.utils.init_coupled_AOADMM_CMTF(Z,'init_options', init_options);
end