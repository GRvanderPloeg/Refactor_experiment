function tests = test_check_data_input
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

    loss_function{1} = 'Frobenius';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
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

    loss_function{1} = 'Frobenius';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
end

function test_regular_PAR2(testCase)
    rng(0,'twister');

    model{1} = 'PAR2';
    sz     = {40,120*ones(1,60),60};
    modes  = {[1 2 3]};
    lambdas_data= {[1 1 1]};

    coupling.lin_coupled_modes = [0 0 0];
    coupling.coupling_type = [];
    coupling.coupl_trafo_matrices = cell(3,1);

    loss_function{1} = 'Frobenius';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
end

function test_irregular_PAR2(testCase)
    rng(0,'twister');

    model{1} = 'PAR2';
    sz     = {40,[61:1:120],60};
    modes  = {[1 2 3]};
    lambdas_data= {[1 1 1]};

    coupling.lin_coupled_modes = [0 0 0];
    coupling.coupling_type = [];
    coupling.coupl_trafo_matrices = cell(3,1);

    loss_function{1} = 'Frobenius';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
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

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
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

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
end

function test_matrix_PAR2(testCase)
    rng(0,'twister');

    model{1} = 'CP';
    model{2} = 'PAR2';
    sz     = {40,60, 40,120*ones(1,50),50};
    modes  = {[1 2], [3 4 5]};
    lambdas_data= {[1 1 1], [1 1 1]};

    coupling.lin_coupled_modes = [1 0 1 0 0];
    coupling.coupling_type = 0;
    coupling.coupl_trafo_matrices = cell(5,1);

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
end

function test_CP_PAR2(testCase)
    rng(0,'twister');

    model{1} = 'CP';
    model{2} = 'PAR2';
    sz     = {20,30,40, 20,30*ones(1,20),20};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data= {[1 1 1], [1 1 1]};

    coupling.lin_coupled_modes = [1 0 0 1 0 0];
    coupling.coupling_type = 0;
    coupling.coupl_trafo_matrices = cell(6,1);

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
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

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
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

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
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

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';
    loss_function{3} = 'Frobenius';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
end

function test_matrix_CP_Poisson(testCase)
    rng(0,'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {50,30,40,50,70};
    modes  = {[1 2 3], [4 5]};
    lambdas_data= {[1 1 1], [1 1 1]};
    shape = 1;
    scale = 1;

    coupling.lin_coupled_modes = [1 0 0 1 0];
    coupling.coupling_type = 0;
    coupling.coupl_trafo_matrices = cell(5,1);

    loss_function{1} = 'KL';
    loss_function{2} = 'KL';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
end

function test_coupling_type2(testCase)
% Two CP tensors coupled via component-space coupling (type 2).
% H_m has R rows (= number of lambdas) and d columns (shared space dimension).
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
    d = 2; % d <= R
    coupling.coupl_trafo_matrices{1} = rand(R, d);
    coupling.coupl_trafo_matrices{4} = rand(R, d);

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
end

function test_coupling_type3(testCase)
% Two CP tensors coupled via latent-space coupling (type 3).
% H_m has sz{m} rows and d columns (shared latent dimension).
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
    d = 5; % shared latent dimension
    coupling.coupl_trafo_matrices{1} = rand(sz{1}, d);
    coupling.coupl_trafo_matrices{4} = rand(sz{4}, d);

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';

    cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model);
end

function test_warning_uncoupled_mode_trafo(testCase)
% A non-empty coupl_trafo_matrix for an uncoupled mode (lin_coupled_modes=0)
% should trigger a warning.
    rng(0,'twister');

    model{1} = 'CP';
    sz     = {60,50,70};
    modes  = {[1 2 3]};
    lambdas_data = {[1 1 1]};

    coupling.lin_coupled_modes = [0 0 0];
    coupling.coupling_type = [];
    coupling.coupl_trafo_matrices = cell(3,1);
    coupling.coupl_trafo_matrices{2} = eye(50); % spurious matrix for uncoupled mode 2

    loss_function{1} = 'Frobenius';

    verifyWarning(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        '');
end

function test_warning_exact_coupling_trafo(testCase)
% A non-empty coupl_trafo_matrix for a mode with coupling_type=0 (exact coupling)
% should trigger a warning.
    rng(0,'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {20,30,40, 20,70};
    modes  = {[1 2 3], [4 5]};
    lambdas_data = {[1 1 1], [1 1 1]};

    coupling.lin_coupled_modes = [1 0 0 1 0];
    coupling.coupling_type = 0;
    coupling.coupl_trafo_matrices = cell(5,1);
    coupling.coupl_trafo_matrices{1} = eye(20); % spurious matrix for exactly coupled mode

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';

    verifyWarning(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        '');
end