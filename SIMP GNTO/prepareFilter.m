% Prepare Filter
function [H,Hs] =prepareFilter(Nele,Nelx,Nely,Nelz,El,rmin)
rmin=rmin/El;
iH = ones(Nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:Nelz
    for i1 = 1:Nelx
        for j1 = 1:Nely
            e1 = (k1-1)*Nelx*Nely + (i1-1)*Nely+j1;  
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),Nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),Nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),Nely)
                        e2 = (k2-1)*Nelx*Nely + (i2-1)*Nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2)); 
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
end