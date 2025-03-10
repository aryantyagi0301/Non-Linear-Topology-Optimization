function out = toVoigt(inp)

    out = zeros(6, 6);


    out(1, 1) = inp(1, 1, 1, 1); 
    out(2, 2) = inp(2, 2, 2, 2); 
    out(3, 3) = inp(3, 3, 3, 3); 
    out(6, 6) = inp(1, 2, 1, 2); 
    out(5, 5) = inp(1, 3, 1, 3); 
    out(4, 4) = inp(2, 3, 2, 3); 

    out(2, 1) = inp(2, 2, 1, 1); 
    out(3, 1) = inp(3, 3, 1, 1); 
    out(6, 1) = inp(1, 2, 1, 1); 
    out(5, 1) = inp(1, 3, 1, 1); 
    out(4, 1) = inp(2, 3, 1, 1); 

    out(1, 2) = inp(1, 1, 2, 2); 
    out(3, 2) = inp(3, 3, 2, 2); 
    out(6, 2) = inp(1, 2, 2, 2); 
    out(5, 2) = inp(1, 3, 2, 2); 
    out(4, 2) = inp(2, 3, 2, 2); 

    out(1, 3) = inp(1, 1, 3, 3); 
    out(2, 3) = inp(2, 2, 3, 3); 
    out(6, 3) = inp(1, 2, 3, 3); 
    out(5, 3) = inp(1, 3, 3, 3); 
    out(4, 3) = inp(2, 3, 3, 3); 

    out(1, 6) = inp(1, 1, 1, 2); 
    out(2, 6) = inp(2, 2, 1, 2); 
    out(3, 6) = inp(3, 3, 1, 2); 
    out(5, 6) = inp(1, 3, 1, 2); 
    out(4, 6) = inp(2, 3, 1, 2); 

    out(1, 5) = inp(1, 1, 1, 3); 
    out(2, 5) = inp(2, 2, 1, 3); 
    out(3, 5) = inp(3, 3, 1, 3); 
    out(6, 5) = inp(1, 2, 1, 3); 
    out(4, 5) = inp(2, 3, 1, 3); 

    out(1, 4) = inp(1, 1, 2, 3); 
    out(2, 4) = inp(2, 2, 2, 3); 
    out(3, 4) = inp(3, 3, 2, 3); 
    out(6, 4) = inp(1, 2, 2, 3); 
    out(5, 4) = inp(1, 3, 2, 3); 
end
