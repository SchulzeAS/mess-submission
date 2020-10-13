function x = mul_gpu(A,B)

C = gpuArray(A);
D = gpuArray(B);

x = C*D;