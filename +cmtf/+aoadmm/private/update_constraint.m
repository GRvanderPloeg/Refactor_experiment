function [oldZ,G] = update_constraint(Z,G,m,rho)
    %updates G.constraint_fac{mm} and G.constraint_dual_fac{mm}
    oldZ = G.constraint_fac{m};
    if length(rho)>1 % if it is the third mode (C,D_k) of Parafac2 model (2ns mode(B_k) is handled in own function)
        G.constraint_fac{m} = feval(Z.prox_operators{m},(G.fac{m} + G.constraint_dual_fac{m}),max(rho));
    else
        G.constraint_fac{m} = feval(Z.prox_operators{m},(G.fac{m} + G.constraint_dual_fac{m}),rho);
    end
    G.constraint_dual_fac{m} = G.constraint_dual_fac{m} + G.fac{m} - G.constraint_fac{m};
end