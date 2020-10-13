function x = solve_gpu(A,B)
%in preparation for future implementation of some solver that does not 
%have a problem with B being NxM, where both N;M are != 1

x = A\B;