function modulus = analytical_modulus_withC(C)
    c1 = 2.0;
    c2 = 3.0;
    beta = c2;

    
    C_inv = inv(C);
    I3 = det(C);
    J = sqrt(det(C));

    
    dyadic_product = zeros(3, 3, 3, 3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    dyadic_product(i, j, k, l) = C_inv(i, k) * C_inv(j, l);
                end
            end
        end
    end
    
    dyadic_product = permute(dyadic_product, [1, 3, 2, 4]);
    
    dyadic_dot = zeros(3, 3, 3, 3);
    for A = 1:3
        for B = 1:3
            for C = 1:3
                for D = 1:3
                    dyadic_dot(A, B, C, D) = 0.5 * (C_inv(A, C) * C_inv(B, D) + C_inv(A, D) * C_inv(B, C));
                end
            end
        end
    end


    modulus = ((2.0 * c1 - 2.0 * c2 * J * (J - 1.0)) * dyadic_dot) + ...
              ((2.0 * c2 * J^2 - c2 * J) * dyadic_product);


    modulus = toVoigt(modulus);
end
