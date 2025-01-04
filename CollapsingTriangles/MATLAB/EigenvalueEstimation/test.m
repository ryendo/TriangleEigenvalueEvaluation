syms t s
K = [s^2 + 1/t,   s*t;
     t + 1/s,     2*t^3];
M = [  s*t,       3*s + 1/t^2;
      4*s^3,      s*t^2];

[fctr, Kscaled, Mscaled] = normalizeMatrices(K, M)
