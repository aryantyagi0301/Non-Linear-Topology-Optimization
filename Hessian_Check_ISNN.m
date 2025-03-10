F = eye(3);

F(1, 1) = F(1, 1) + 0.2;
F(1, 2) = F(1, 2) - 0.1;
F(1, 3) = F(1, 3) + 0.11;

F(2, 1) = F(2, 1) - 0.05;
F(2, 2) = F(2, 2) + 0.1;
F(2, 3) = F(2, 3) + 0.21;

F(3, 1) = F(3, 1) - 0.24;
F(3, 2) = F(3, 2) + 0.14;
F(3, 3) = F(3, 3) - 0.1;

C = F' * F;
analytical_modulus_withC(C)