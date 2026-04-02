function [rel_primal_res_coupling,rel_dual_res_coupling] = eval_res_ADMM_coupl(Z,G,modes,coupl_id,ctype,oldDelta)
% Computes relative primal and dual residuals of ADMM for all coupling types.
%
% Primal residual: distance from the coupling constraint (e.g. H*C = Delta).
% Dual residual:   change in Delta (scaled by H for types 3 and 4).

    rel_primal_res_coupling = 0;
    rel_dual_res_coupling   = 0;
    Delta    = G.coupling_fac{coupl_id};
    dDelta   = Delta - oldDelta;

    for mm = modes
        C  = G.fac{mm};
        mu = G.coupling_dual_fac{mm};

        switch ctype
            case 0
                primal_vec       = C - Delta;
                primal_norm_ref  = norm(C, 'fro');
                dual_vec         = dDelta;
            case 1
                HC               = Z.coupling.coupl_trafo_matrices{mm}*C;
                primal_vec       = HC - Delta;
                primal_norm_ref  = norm(HC, 'fro');
                dual_vec         = dDelta;
            case 2
                CH               = C*Z.coupling.coupl_trafo_matrices{mm};
                primal_vec       = CH - Delta;
                primal_norm_ref  = norm(CH, 'fro');
                dual_vec         = dDelta;
            case 3
                primal_vec       = C - Z.coupling.coupl_trafo_matrices{mm}*Delta;
                primal_norm_ref  = norm(C, 'fro');
                dual_vec         = Z.coupling.coupl_trafo_matrices{mm}*dDelta;
            case 4
                primal_vec       = C - Delta*Z.coupling.coupl_trafo_matrices{mm};
                primal_norm_ref  = norm(C, 'fro');
                dual_vec         = dDelta*Z.coupling.coupl_trafo_matrices{mm};
        end

        rel_primal_res_coupling = rel_primal_res_coupling + norm(primal_vec,'fro') / primal_norm_ref;
        scaling = norm(mu,'fro');
        if scaling > 0
            rel_dual_res_coupling = rel_dual_res_coupling + norm(dual_vec,'fro') / scaling;
        else
            rel_dual_res_coupling = rel_dual_res_coupling + norm(dual_vec,'fro');
        end
    end

    rel_primal_res_coupling = rel_primal_res_coupling / length(modes);
    rel_dual_res_coupling   = rel_dual_res_coupling   / length(modes);
end
