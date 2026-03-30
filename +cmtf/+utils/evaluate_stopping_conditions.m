function [stop] = evaluate_stopping_conditions(f_tensors,f_couplings,f_constraints,f_PAR2_couplings,f_tensors_old,f_couplings_old,f_constraints_old,f_PAR2_couplings_old,options)
% Evaluate convergence stopping conditions for the AOADMM outer loop.
%
% Returns true if ALL four objective components (tensor fit, coupling,
% constraints, PARAFAC2 coupling) have individually converged, where
% convergence means either the value is below options.AbsFuncTol or its
% relative change from the previous iteration is below options.OuterRelTol.
%
% Syntax:
%   stop = cmtf.utils.evaluate_stopping_conditions( ...
%              f_tensors, f_couplings, f_constraints, f_PAR2_couplings, ...
%              f_tensors_old, f_couplings_old, f_constraints_old, f_PAR2_couplings_old, ...
%              options)
%
% Inputs:
%   f_tensors, f_couplings, f_constraints, f_PAR2_couplings         - Current objective components
%   f_tensors_old, f_couplings_old, f_constraints_old, f_PAR2_couplings_old - Previous iteration values
%   options - Struct with fields AbsFuncTol and OuterRelTol
%
% Output:
%   stop - Logical scalar; true if all components have converged
stop_tensors = false;
stop_couplings = false;
stop_constraints = false;
stop_PAR2_couplings = false;

if f_tensors_old>0
    f_tensors_rel_change = abs(f_tensors_old-f_tensors)/f_tensors_old;
else
    f_tensors_rel_change = abs(f_tensors_old-f_tensors);
end
if f_tensors < options.AbsFuncTol || f_tensors_rel_change < options.OuterRelTol
    stop_tensors = true;
end

if f_couplings_old>0
    f_couplings_rel_change = abs(f_couplings_old-f_couplings)/f_couplings_old;
else
    f_couplings_rel_change = abs(f_couplings_old-f_couplings);
end
if f_couplings < options.AbsFuncTol || f_couplings_rel_change < options.OuterRelTol
    stop_couplings = true;
end

if f_constraints_old>0
    f_constraints_rel_change = abs(f_constraints_old-f_constraints)/f_constraints_old;
else
    f_constraints_rel_change = abs(f_constraints_old-f_constraints);
end
if f_constraints < options.AbsFuncTol || f_constraints_rel_change < options.OuterRelTol
    stop_constraints = true;
end

if f_PAR2_couplings_old>0
    f_PAR2_couplings_rel_change = abs(f_PAR2_couplings_old-f_PAR2_couplings)/f_PAR2_couplings_old;
else
    f_PAR2_couplings_rel_change = abs(f_PAR2_couplings_old-f_PAR2_couplings);
end
if f_PAR2_couplings < options.AbsFuncTol || f_PAR2_couplings_rel_change < options.OuterRelTol
    stop_PAR2_couplings = true;
end

stop = stop_tensors & stop_couplings & stop_constraints & stop_PAR2_couplings;


end

