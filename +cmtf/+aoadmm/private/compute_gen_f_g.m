function [function_value,gradient_vector] = compute_gen_f_g(Z,G,x,p,m,ff,gg,constrained,coupling_type,rho,options)
    % constrained: boolean
    % coupling_type: 0,...,4 or something else for not coupled
    MM = update(ktensor(G.fac(Z.modes{p})),find(Z.modes{p}==m),x);
    Mfull = full(MM);
    function_value = Z.weights(p)*collapse(ff(Z.object{p},Mfull));
    Y = gg(Z.object{p},Mfull);
    Gmatrix = mttkrp(Y,MM.U,find(Z.modes{p}==m));
    gradient_vector = Z.weights(p).*Gmatrix(:); %vectorize
    if constrained
        function_value = function_value + rho/2*sum((x-G.constraint_fac{m}(:)+G.constraint_dual_fac{m}(:)).^2);
        gradient_vector = gradient_vector + rho*(x-G.constraint_fac{m}(:)+G.constraint_dual_fac{m}(:));
    end
    n = Z.coupling.lin_coupled_modes(m);
    switch coupling_type
        case 0
            function_value = function_value + rho/2*sum((x-G.coupling_fac{n}(:)+G.coupling_dual_fac{m}(:)).^2);
            gradient_vector = gradient_vector + rho*(x-G.coupling_fac{n}(:)+G.coupling_dual_fac{m}(:));
        case 1
            function_value = function_value + rho/2*sum(sum(Z.coupling.coupl_trafo_matrices{m}*reshape(x,size(G.fac{m}))-G.coupling_fac{n}+G.coupling_dual_fac{m}).^2); 
            gradient_vector = gradient_vector + rho*reshape(Z.coupling.coupl_trafo_matrices{m}'*(Z.coupling.coupl_trafo_matrices{m}*reshape(x,size(G.fac{m}))-G.coupling_fac{n}+G.coupling_dual_fac{m}),[],1);
        case 2
            function_value = function_value + rho/2*sum(sum(reshape(x,size(G.fac{m}))*Z.coupling.coupl_trafo_matrices{m}-G.coupling_fac{n}+G.coupling_dual_fac{m}).^2);
            gradient_vector = gradient_vector + rho*reshape((reshape(x,size(G.fac{m}))*Z.coupling.coupl_trafo_matrices{m}-G.coupling_fac{n}+G.coupling_dual_fac{m})*Z.coupling.coupl_trafo_matrices{m}',[],1);
        case 3
            function_value = function_value + rho/2*sum((x-reshape(Z.coupling.coupl_trafo_matrices{m}*G.coupling_fac{n},[],1)+G.coupling_dual_fac{m}(:)).^2);
            gradient_vector = gradient_vector + rho*(x-reshape(Z.coupling.coupl_trafo_matrices{m}*G.coupling_fac{n},[],1)+G.coupling_dual_fac{m}(:));
        case 4
            function_value = function_value + rho/2*sum((x-reshape(G.coupling_fac{n}*Z.coupling.coupl_trafo_matrices{m},[],1)+G.coupling_dual_fac{m}(:)).^2);
            gradient_vector = gradient_vector + rho*(x-reshape(G.coupling_fac{n}*Z.coupling.coupl_trafo_matrices{m},[],1)+G.coupling_dual_fac{m}(:));
    end
    if isfield(Z,'ridge') % If f=lambda*sum(x.^2) then surely df/fx=2*lambda*x and NOT lambda/2*x??? So the gradient should be 2*Z.ridge(m)*x
        function_value = function_value + Z.ridge(m)*sum(x.^2);
        gradient_vector = gradient_vector + Z.ridge(m)/2*x;
    end
    if options.bsum
        function_value = function_value + options.bsum_weight/2*sum((x-G.fac{m}(:)).^2);
        gradient_vector = gradient_vector + options.bsum_weight*(x-G.fac{m}(:));
    end
end