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

% =========================================================================
% Error-branch tests
% =========================================================================

function test_error_sz_modes_mismatch(testCase)
% max mode index exceeds length(sz)
    model = {'CP'};
    sz    = {60, 50};          % only 2 mode sizes
    modes = {[1 2 3]};         % but mode 3 referenced
    lambdas_data = {[1 1 1]};
    coupling.lin_coupled_modes   = [0 0 0];
    coupling.coupling_type       = [];
    coupling.coupl_trafo_matrices = cell(3,1);
    loss_function = {'Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_coupling_count_mismatch(testCase)
% max(lin_coupled_modes)=1 but coupling_type has 2 entries
    model = {'CP','CP'};
    sz    = {20,30,40, 20,70,80};
    modes = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,3), ones(1,3)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = [0 0];  % should be length 1
    coupling.coupl_trafo_matrices = cell(6,1);
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_PAR2_size_mismatch(testCase)
% PAR2: sz{modes{p}(3)} must equal length(sz{modes{p}(2)})
    model = {'PAR2'};
    sz    = {40, 120*ones(1,60), 50};  % sz{3}=50 but length(sz{2})=60
    modes = {[1 2 3]};
    lambdas_data = {[1 1 1]};
    coupling.lin_coupled_modes   = [0 0 0];
    coupling.coupling_type       = [];
    coupling.coupl_trafo_matrices = cell(3,1);
    loss_function = {'Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_PAR2_non_frobenius(testCase)
% PAR2 model only supports Frobenius loss
    model = {'PAR2'};
    sz    = {40, 120*ones(1,60), 60};
    modes = {[1 2 3]};
    lambdas_data = {[1 1 1]};
    coupling.lin_coupled_modes   = [0 0 0];
    coupling.coupling_type       = [];
    coupling.coupl_trafo_matrices = cell(3,1);
    loss_function = {'KL'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_PAR2_mode2_coupled(testCase)
% Coupling in mode 2 (varying mode) of PAR2 is not supported
    model = {'CP','PAR2'};
    sz    = {40,60, 40,120*ones(1,50),50};
    modes = {[1 2], [3 4 5]};
    lambdas_data = {[1 1 1], [1 1 1]};
    coupling.lin_coupled_modes   = [1 0 1 1 0];  % mode 4 (mode-2 of PAR2) is coupled
    coupling.coupling_type       = [0 0];
    coupling.coupl_trafo_matrices = cell(5,1);
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type0_sz_mismatch(testCase)
% Type 0 (exact coupling): coupled modes must have same sz
    model = {'CP','CP'};
    sz    = {20,30,40, 25,70,80};  % sz{1}=20, sz{4}=25 — mismatch
    modes = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,3), ones(1,3)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 0;
    coupling.coupl_trafo_matrices = cell(6,1);
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type0_lambda_mismatch(testCase)
% Type 0: coupled modes must have same number of components
    model = {'CP','CP'};
    sz    = {20,30,40, 20,70,80};
    modes = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,3), ones(1,4)};  % 3 vs 4 components
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 0;
    coupling.coupl_trafo_matrices = cell(6,1);
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type1_missing_trafo(testCase)
% Type 1: coupling matrix must be provided
    model = {'CP','CP'};
    sz    = {20,30,40, 40,70,80};
    modes = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,3), ones(1,3)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 1;
    coupling.coupl_trafo_matrices = cell(6,1);  % both empty
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type1_not_right_invertible(testCase)
% Type 1: coupling matrix must be right-invertible (full row rank)
    model = {'CP','CP'};
    sz    = {20,30,40, 40,70,80};
    modes = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,3), ones(1,3)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 1;
    coupling.coupl_trafo_matrices = cell(6,1);
    % Tall matrix: more rows than columns -> rank < rows -> not right-invertible
    coupling.coupl_trafo_matrices{1} = [eye(20); zeros(1,20)];  % 21x20, rank=20 < 21
    coupling.coupl_trafo_matrices{4} = [eye(20); zeros(1,20)];
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type1_sz_cols_mismatch(testCase)
% Type 1: number of columns of trafo must equal sz{m}
    model = {'CP','CP'};
    sz    = {20,30,40, 40,70,80};
    modes = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,3), ones(1,3)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 1;
    coupling.coupl_trafo_matrices = cell(6,1);
    coupling.coupl_trafo_matrices{1} = eye(20,20);   % ok: 20 cols == sz{1}=20
    coupling.coupl_trafo_matrices{4} = eye(20,30);   % bad: 30 cols ~= sz{4}=40
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type1_trafo_rows_mismatch(testCase)
% Type 1: both coupled trafos must have same number of rows
    model = {'CP','CP'};
    sz    = {20,30,40, 40,70,80};
    modes = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,3), ones(1,3)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 1;
    coupling.coupl_trafo_matrices = cell(6,1);
    coupling.coupl_trafo_matrices{1} = eye(20,20);   % 20 rows
    coupling.coupl_trafo_matrices{4} = zeros(15,40); % 15 rows — mismatch
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type1_lambda_mismatch(testCase)
% Type 1: coupled modes must have same number of components
    model = {'CP','CP'};
    sz    = {20,30,40, 40,70,80};
    modes = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,3), ones(1,4)};  % 3 vs 4
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 1;
    coupling.coupl_trafo_matrices = cell(6,1);
    coupling.coupl_trafo_matrices{1} = eye(20,20);
    coupling.coupl_trafo_matrices{4} = eye(20,40);
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type2_missing_trafo(testCase)
% Type 2: coupling matrix must be provided
    model = {'CP','CP'};
    sz    = {20,30,40, 20,70,80};
    modes = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,3), ones(1,3)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 2;
    coupling.coupl_trafo_matrices = cell(6,1);  % both empty
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type2_lambda_rows_mismatch(testCase)
% Type 2: rows of trafo must equal length(lambdas{p})
    model = {'CP','CP'};
    sz    = {20,30,40, 20,70,80};
    modes = {[1 2 3], [4 5 6]};
    R = 3; d = 2;
    lambdas_data = {ones(1,R), ones(1,R)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 2;
    coupling.coupl_trafo_matrices = cell(6,1);
    coupling.coupl_trafo_matrices{1} = rand(R+1, d);  % R+1 rows ~= R lambdas
    coupling.coupl_trafo_matrices{4} = rand(R+1, d);
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type2_sz_mismatch(testCase)
% Type 2: coupled modes must have same sz
    model = {'CP','CP'};
    sz    = {20,30,40, 25,70,80};  % sz{1}=20, sz{4}=25 — mismatch
    modes = {[1 2 3], [4 5 6]};
    R = 3; d = 2;
    lambdas_data = {ones(1,R), ones(1,R)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 2;
    coupling.coupl_trafo_matrices = cell(6,1);
    coupling.coupl_trafo_matrices{1} = rand(R, d);
    coupling.coupl_trafo_matrices{4} = rand(R, d);
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type2_trafo_cols_mismatch(testCase)
% Type 2: both coupled trafos must have same number of columns
    model = {'CP','CP'};
    sz    = {20,30,40, 20,70,80};
    modes = {[1 2 3], [4 5 6]};
    R = 3;
    lambdas_data = {ones(1,R), ones(1,R)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 2;
    coupling.coupl_trafo_matrices = cell(6,1);
    coupling.coupl_trafo_matrices{1} = rand(R, 2);  % 2 cols
    coupling.coupl_trafo_matrices{4} = rand(R, 3);  % 3 cols — mismatch
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type2_more_cols_than_rows(testCase)
% Type 2: trafo cannot have more columns than rows (d <= R)
    model = {'CP','CP'};
    sz    = {20,30,40, 20,70,80};
    modes = {[1 2 3], [4 5 6]};
    R = 3;
    lambdas_data = {ones(1,R), ones(1,R)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 2;
    coupling.coupl_trafo_matrices = cell(6,1);
    d = R + 2;  % d > R — invalid
    coupling.coupl_trafo_matrices{1} = rand(R, d);
    coupling.coupl_trafo_matrices{4} = rand(R, d);
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type3_missing_trafo(testCase)
% Type 3: coupling matrix must be provided
    model = {'CP','CP'};
    sz    = {20,30,40, 25,70,80};
    modes = {[1 2 3], [4 5 6]};
    R = 3;
    lambdas_data = {ones(1,R), ones(1,R)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 3;
    coupling.coupl_trafo_matrices = cell(6,1);  % both empty
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type3_sz_rows_mismatch(testCase)
% Type 3: rows of trafo must equal sz{m}
    model = {'CP','CP'};
    sz    = {20,30,40, 25,70,80};
    modes = {[1 2 3], [4 5 6]};
    R = 3; d = 5;
    lambdas_data = {ones(1,R), ones(1,R)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 3;
    coupling.coupl_trafo_matrices = cell(6,1);
    coupling.coupl_trafo_matrices{1} = rand(sz{1}+1, d);  % rows ~= sz{1}
    coupling.coupl_trafo_matrices{4} = rand(sz{4}, d);
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type3_trafo_cols_mismatch(testCase)
% Type 3: both coupled trafos must have same number of columns
    model = {'CP','CP'};
    sz    = {20,30,40, 25,70,80};
    modes = {[1 2 3], [4 5 6]};
    R = 3;
    lambdas_data = {ones(1,R), ones(1,R)};
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 3;
    coupling.coupl_trafo_matrices = cell(6,1);
    coupling.coupl_trafo_matrices{1} = rand(sz{1}, 5);  % 5 cols
    coupling.coupl_trafo_matrices{4} = rand(sz{4}, 4);  % 4 cols — mismatch
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type3_lambda_mismatch(testCase)
% Type 3: coupled modes must have same number of components
    model = {'CP','CP'};
    sz    = {20,30,40, 25,70,80};
    modes = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,3), ones(1,4)};  % 3 vs 4
    coupling.lin_coupled_modes   = [1 0 0 1 0 0];
    coupling.coupling_type       = 3;
    coupling.coupl_trafo_matrices = cell(6,1);
    coupling.coupl_trafo_matrices{1} = rand(sz{1}, 5);
    coupling.coupl_trafo_matrices{4} = rand(sz{4}, 5);
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type4_missing_trafo(testCase)
% Type 4: coupling matrix must be provided
    model = {'CP','CP'};
    sz    = {50,30,40, 50,70};
    modes = {[1 2 3], [4 5]};
    lambdas_data = {ones(1,4), ones(1,3)};
    coupling.lin_coupled_modes   = [1 0 0 1 0];
    coupling.coupling_type       = 4;
    coupling.coupl_trafo_matrices = cell(5,1);  % both empty
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type4_lambda_cols_mismatch(testCase)
% Type 4: columns of trafo must equal length(lambdas{p})
    model = {'CP','CP'};
    sz    = {50,30,40, 50,70};
    modes = {[1 2 3], [4 5]};
    lambdas_data = {ones(1,4), ones(1,3)};
    coupling.lin_coupled_modes   = [1 0 0 1 0];
    coupling.coupling_type       = 4;
    coupling.coupl_trafo_matrices = cell(5,1);
    coupling.coupl_trafo_matrices{1} = eye(5,5);          % 5 cols ~= 4 lambdas
    coupling.coupl_trafo_matrices{4} = [eye(3,3);0 0 0];
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type4_sz_mismatch(testCase)
% Type 4: coupled modes must have same sz
    model = {'CP','CP'};
    sz    = {50,30,40, 55,70};  % sz{1}=50, sz{4}=55 — mismatch
    modes = {[1 2 3], [4 5]};
    lambdas_data = {ones(1,4), ones(1,3)};
    coupling.lin_coupled_modes   = [1 0 0 1 0];
    coupling.coupling_type       = 4;
    coupling.coupl_trafo_matrices = cell(5,1);
    coupling.coupl_trafo_matrices{1} = eye(4,4);
    coupling.coupl_trafo_matrices{4} = [eye(3,3);0 0 0];
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end

function test_error_type4_trafo_rows_mismatch(testCase)
% Type 4: both coupled trafos must have same number of rows
    model = {'CP','CP'};
    sz    = {50,30,40, 50,70};
    modes = {[1 2 3], [4 5]};
    lambdas_data = {ones(1,4), ones(1,3)};
    coupling.lin_coupled_modes   = [1 0 0 1 0];
    coupling.coupling_type       = 4;
    coupling.coupl_trafo_matrices = cell(5,1);
    coupling.coupl_trafo_matrices{1} = eye(4,4);    % 4 rows
    coupling.coupl_trafo_matrices{4} = eye(3,3);    % 3 rows — mismatch
    loss_function = {'Frobenius','Frobenius'};
    verifyError(testCase, ...
        @() cmtf.utils.check_data_input(sz,modes,lambdas_data,coupling,loss_function,model), ...
        ?MException);
end