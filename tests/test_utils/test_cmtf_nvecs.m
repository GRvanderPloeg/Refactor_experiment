function tests = test_cmtf_nvecs
    tests = functiontests(localfunctions);
end

function test_tensor_mode1(testCase)
% Dense tensor: singular vectors along mode 1.
    rng(0, 'twister');
    sz1 = 15; sz2 = 10; sz3 = 12;
    r = 2;
    Z.object{1} = tensor(rand(sz1, sz2, sz3));
    Z.modes{1}  = [1 2 3];

    U = cmtf.utils.cmtf_nvecs(Z, 1, r);

    verifyEqual(testCase, size(U), [sz1, r]);
end

function test_tensor_mode2(testCase)
% Dense tensor: singular vectors along mode 2.
    rng(1, 'twister');
    sz1 = 15; sz2 = 10; sz3 = 12;
    r = 3;
    Z.object{1} = tensor(rand(sz1, sz2, sz3));
    Z.modes{1}  = [1 2 3];

    U = cmtf.utils.cmtf_nvecs(Z, 2, r);

    verifyEqual(testCase, size(U), [sz2, r]);
end

function test_matrix_mode1(testCase)
% Plain matrix as object, requesting mode-1 vectors (uses double(obj) path).
    rng(2, 'twister');
    sz1 = 20; sz2 = 15;
    r = 2;
    Z.object{1} = rand(sz1, sz2); % plain matrix, not a tensor object
    Z.modes{1}  = [1 2];

    U = cmtf.utils.cmtf_nvecs(Z, 1, r);

    verifyEqual(testCase, size(U), [sz1, r]);
end

function test_matrix_mode2(testCase)
% Plain matrix as object, requesting mode-2 vectors (uses double(obj)' path).
    rng(3, 'twister');
    sz1 = 20; sz2 = 15;
    r = 2;
    Z.object{1} = rand(sz1, sz2);
    Z.modes{1}  = [1 2];

    U = cmtf.utils.cmtf_nvecs(Z, 2, r);

    verifyEqual(testCase, size(U), [sz2, r]);
end

function test_coupled_tensors(testCase)
% Two tensors sharing mode 1: concatenated unfoldings should be used.
    rng(4, 'twister');
    sz1 = 12; sz2 = 8; sz3 = 6; sz4 = 9; sz5 = 7;
    r = 2;
    Z.object{1} = tensor(rand(sz1, sz2, sz3));
    Z.object{2} = tensor(rand(sz1, sz4, sz5));
    Z.modes{1}  = [1 2 3];
    Z.modes{2}  = [1 4 5];

    U = cmtf.utils.cmtf_nvecs(Z, 1, r);

    verifyEqual(testCase, size(U), [sz1, r]);
end
