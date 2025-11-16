function S = prox_ment(St,went)

limit = 5e2;
 
T1 = -went.*log(went) + St - went;
L1 = -log(went) + St./went - 1;
L2 = log(L1);
T2 = went.*L2;
S2 = T1 - T2;

C = St./went -1 -log(went);
  
S2(C<=limit) = went.*lambertw(exp(C(C<=limit)));
S = S2;
 