function [exit_flag] = make_exit_flag(iter,f_tensors,f_couplings,f_constraints,f_PAR2_couplings,options,illconditioned)
% Determine the stopping condition that terminated the AOADMM outer loop.
%
% Returns 'maxIterations' if the iteration limit was reached, a string if
% an ill-conditioned linear system was encountered, or a struct with one
% field per objective component indicating whether it stopped due to
% 'AbsFuncTol' or 'RelFuncTol'.
%
% Syntax:
%   exit_flag = cmtf.utils.make_exit_flag(iter, f_tensors, f_couplings, ...
%                   f_constraints, f_PAR2_couplings, options, illconditioned)
%
% Inputs:
%   iter             - Final iteration count (one past the last completed iteration)
%   f_tensors        - Final tensor fit objective value
%   f_couplings      - Final coupling objective value
%   f_constraints    - Final constraint objective value
%   f_PAR2_couplings - Final PARAFAC2 coupling objective value
%   options          - Struct with fields MaxOuterIters and AbsFuncTol
%   illconditioned   - Logical flag; true if a near-singular system was detected
%
% Output:
%   exit_flag - String 'maxIterations', string 'illconditioned lin system',
%               or struct with fields f_tensors, f_couplings, f_constraints,
%               f_PAR2_couplings each set to 'AbsFuncTol' or 'RelFuncTol'

if iter>options.MaxOuterIters
    exit_flag = 'maxIterations';
elseif illconditioned 
    exit_flag = "illconditioned lin system";
else
    if f_tensors < options.AbsFuncTol 
        exit_flag.f_tensors = "AbsFuncTol";
    else
        exit_flag.f_tensors = "RelFuncTol";
    end
    if f_couplings < options.AbsFuncTol 
        exit_flag.f_couplings = "AbsFuncTol";
    else
        exit_flag.f_couplings = "RelFuncTol";
    end
    if f_constraints < options.AbsFuncTol 
        exit_flag.f_constraints = "AbsFuncTol";
    else
        exit_flag.f_constraints = "RelFuncTol";
    end
    if f_PAR2_couplings < options.AbsFuncTol
        exit_flag.f_PAR2_couplings = "AbsFuncTol";
    else
        exit_flag.f_PAR2_couplings = "RelFuncTol";
    end
end
   
end

