function [rel_primal_res_coupling,rel_dual_res_coupling] = eval_res_ADMM_coupl_case0(G,modes,coupl_id,oldDelta)
% computes relative primal and dual residuals of ADMM iteration for factor
% matrices in in coupled modes connected (coupling case 0 only!)
    rel_primal_res_coupling = 0;
    rel_dual_res_coupling = 0;
    for mm=modes
        rel_primal_res_coupling = rel_primal_res_coupling + norm(G.fac{mm}-G.coupling_fac{coupl_id},'fro')/norm(G.fac{mm},'fro');
        scaling = norm(G.coupling_dual_fac{mm},'fro');
        if scaling>0
            rel_dual_res_coupling = rel_dual_res_coupling + norm(G.coupling_fac{coupl_id}-oldDelta,'fro')/scaling;
        else
            rel_dual_res_coupling = rel_dual_res_coupling + norm(G.coupling_fac{coupl_id}-oldDelta,'fro');
        end
    end
    rel_primal_res_coupling = rel_primal_res_coupling/length(modes);
    rel_dual_res_coupling = rel_dual_res_coupling/length(modes);
end