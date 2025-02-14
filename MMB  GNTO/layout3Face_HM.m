% Initialize the layout of the design domain by MMBs
function [variable,N]=layout3Face_HM(DL,DW,DH,xn,yn,zn,rk)
% Design domain of length ,width and height
L=DL-2*rk;  x_int=L/(xn);
W=DW-2*rk;  y_int=W/(yn);
H=DH-2*rk;  z_int=H/(zn);
% Bars layout in x-direction
xh1=rk:x_int:(L-x_int+rk);             
xh1=kron(xh1,ones(1,(zn+1)*(yn+1)));   
xh2=(rk+x_int):x_int:L+rk;    
xh2=kron(xh2,ones(1,(zn+1)*(yn+1)));   
yh=DW/2;   
yh=repmat(kron(yh,ones(1,zn+1)),1,xn);    
zh=rk:z_int:H+rk;
zh=repmat(repmat(zh,1,yn+1),1,xn);  
% Bars layout in z-direction
xv=rk:x_int:L+rk;
xv=repmat(kron(xv,ones(1,zn)),1,yn+1);      
yv=DW/2:y_int:(W+rk); 
yv=kron(yv,ones(1,zn*(xn+1)));
zv1=rk:z_int:H-z_int+rk;
zv1=repmat(repmat(zv1,1,xn+1),1,yn+1);       
zv2=z_int+rk:z_int:H+rk;
zv2=repmat(repmat(zv2,1,xn+1),1,yn+1);        
% Bars layout in y-direction
xy=rk:x_int:L+rk;
xy=kron(xy,ones(1,yn*(zn+1)));          
yy1=rk:y_int:W-y_int+rk;
yy1=repmat(repmat(yy1,1,zn+1),1,xn+1);        
yy2=y_int+rk:y_int:W+rk;
yy2=repmat(repmat(yy2,1,zn+1),1,xn+1);       
zy=rk:z_int:H+rk;
zy=repmat(kron(zy,ones(1,yn)),1,xn+1);
N=length(xh1)+length(xv)+length(xy);       
Rk=repmat(rk,1,N);            
variable=[xh1,xv,xy; yh,yv,yy1; zh,zv1,zy; xh2,xv,xy; yh,yv,yy2; zh,zv2,zy; Rk];
end