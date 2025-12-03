function [lbfgsb_iterations,G] = lbfgsb_update(Z,G,lscalar,uscalar,fh,gh,p,m,constrained,coupling_type,rho,options)
    % updates G.fac{m} using lbfgsb
    lbfgsb_options = options.lbfgsb_options;
    ll = lscalar{p}*ones(numel(G.fac{m}),1);
    uu = uscalar{p}*ones(numel(G.fac{m}),1);
    lbfgsb_loss_func_inner = @(x) compute_gen_f_g(Z,G,x,p,m,fh{p},gh{p},constrained,coupling_type,rho,options);
    lbfgsb_options.x0 = G.fac{m}(:); %vectorize starting point
    [xfacm_,~,lbfgsb_info] = lbfgsb(lbfgsb_loss_func_inner,ll,uu,lbfgsb_options);
    lbfgsb_iterations = lbfgsb_info.iterations;
    G.fac{m} = reshape(xfacm_,size(G.fac{m}));
end