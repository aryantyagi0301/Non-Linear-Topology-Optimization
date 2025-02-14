% To remove bars with particularly small thicknesses
function [xy00,xold,xold1,xmin,xmax,low,upp,N,nn,xmma]=deleteBar(xy00,xold,xold1,xmin,xmax,low,upp,N,Var_num,xmma)
bb=0;    %Calculate the number of bars to be deleted
t00=0.1; % Minimum value for the retained bar thickness
nf=[];   % Number of design variables 
for k=1 : N
    aa=xy00(Var_num*k);
    if aa<t00 % Control the minimum inner radius of the bars
        bb=bb+1;  fff=Var_num*k-Var_num+1:Var_num*k; nf=[fff nf];
    end
end
xy00(nf)=[]; xold(nf)=[]; xold1(nf)=[]; xmin(nf)=[];  xmax(nf)=[]; low(nf)=[]; upp(nf)=[];xmma(nf)=[];
N=N-bb; nn=Var_num*N; % The total number of design variables
end