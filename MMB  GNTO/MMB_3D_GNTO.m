% Geometric Nonlinear Topology Optimization base on MMB
function MMB_3D_GNTO(DL,DW,DH,nelx,nely,nelz,xn,yn,zn,volfrac,rk)
penal=3;
scale=1.0;  % Objective function amplification factor for MMA
E=3e9;    % Young¡¯s modulus
nu=0.4; % Poisson ratio
xthreshold=0.01; % Low density threshold
force=200 ;  % Magnitude of external force
IterMax=16;  % Maximum number of Newton-Raphson iterations
pmin=1e-4;
level = 0.6; % Density isosurface values
parent_dir_name ='GNTO results'; % File name for saving the figures
if  exist(parent_dir_name)
rmdir(parent_dir_name, 's') % Delete the existed file
end
mkdir(parent_dir_name);
% Calculate the coordinates of the center points for all elements
EL=DL/nelx; EW=DW/nely; EH=DH/nelz;
[ey,ex,ez]=meshgrid(EW * [0.5 : nely],EL * [0.5 : nelx],EH * [0.5 : nelz]);
EC.dim = 3;EC.xs = cell(EC.dim, 1);
EC.xs{1}=ex; EC.xs{2}=ey; EC.xs{3}=ez;
% Calculate all nodal coordinates
[ey1,ex1,ez1] = meshgrid(0:nely, 0:nelx, 0:nelz);
Nodes=[ex1(:) ey1(:) ez1(:)];
% Initial layout of the design domain by MMBs
[variable,N]=layout3Face_HM(DL,DW,DH,xn,yn,zn,rk);
%  Maximum and minimum values for the design variables
xmin=[ 0.0 ; 0.0 ;  0.0 ; 0.0 ;0.0 ; 0.0 ; -5];
[Var_num,~]=size(xmin);
xmin=repmat(xmin,N,1);
xmax=[DL ;DW;DH ;DL ;DW; DH;  DH/4];
xmax=repmat(xmax,N,1);
% MMA parameters
m=1;   nn=Var_num*N;
[low,upp]=deal(ones(nn,1));
c=5000*ones(m,1);
a0=1;
[d,b]=deal(zeros(m,1));
xy00=variable(:);
[xold1,xold2]=deal(xy00);
% Nodal DOFs for the loading force
loadnid=(nelx+1)*(nely+1)*nelz+(nelx/2+1):(nelx+1):(nelx+1)*(nely+1)*(nelz+1);
loaddof = 3*loadnid;
% Nodal DOFs at the fixed end
fixednidL =(nelx+1)*(nely+1)*nelz/2+1:(nelx+1):(nelx+1)*(nely+1)*(nelz/2+1);
fixednidR=fixednidL+nelx;
fixednid =[fixednidL fixednidR];
Fixeddofs = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2];
% Prepare for FEA
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
Alldofs=1:ndof;
Freedofs=setdiff(Alldofs,Fixeddofs);
ForceVector = sparse(loaddof,1,force,ndof,1);
nodegrd = reshape(1:(nelx+1)*(nely+1),nelx+1,nely+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:);
edofMat=repmat(edofVec,1,24)+repmat([-2 -1 0 1 2 3 3*nelx+[4 5 6 1 2 3] ...
    3*(nely+1)*(nelx+1)+[-2 -1 0 1 2 3 [3*nelx+[4 5 6 1 2 3]]]] ,nele,1);
Eles= repmat(nodeids(:),1,8)+repmat([0 1 nelx+[2 1] ...
    (nely+1)*(nelx+1)+[0 1 nelx+[2 1]]],nele,1);
[Nodes2Ele]=NodeSurroundingElement(Nodes,Eles); % Find the adjacent element for each node
Loop=1; change=1; maxloop=101; Lamada=zeros(ndof,1);
while change>0.00001 && Loop<maxloop
    % Calculate the element density by MMBs
    bata=min(4+Loop*(4-2)/30,6);
    alctr=min(0.993+Loop*0.005/50,0.9999); %MMA parameter
    pa=1.;
    dek=cell(N,1);
    dekdx=cell(N,1);
    for k=1:N
        [dek{k},dekdx{k}]=DEK(xy00(Var_num*k-Var_num+1:Var_num*k),ex(:),ey(:),ez(:));
        mm=-bata*(dek{k}-xy00(Var_num*k));
        pa=pa*1./(1.+exp(mm));
    end
    PE=1.-pa;
    Eledensity=reshape(PE,nelx,nely,nelz);
    % Show the figure for the optimized result
    clf;   h = visualizeLevelSet(EC, Eledensity, 'surface', level);
    FileName=[parent_dir_name,'\Fig_',int2str(Loop),'.png'];
    saveas(h,FileName);
    % Find nodes surrounded by low-density elements
    [uselessNode]=NodeSurroundWithVoidElement(Nodes2Ele, PE,xthreshold);
    % Geometrically nonlinear FEA + Residual force sensitivity with the element density
    [U,GKF,dRdpe]=GNFEA(ndof,nele,Freedofs,IterMax,ForceVector,Nodes,Eles,edofMat,E,nu,PE,penal,uselessNode,pmin);
    Lamada(Freedofs,:)=GKF(Freedofs,Freedofs)\ForceVector(Freedofs,:);
    dfdpe=full(-Lamada'*dRdpe)';
    % Objective function and sensitivity analysis
    Comp=ForceVector'*U; % End-compliance
    df0dx=zeros(Var_num*N,1); % The sensitivity of the objective function
    dfdx=zeros(Var_num*N,1); % The sensitivity of the constraint function
    for k=1:N
        penal=3;
        mid=exp(-bata*(dek{k}-xy00(Var_num*k)));
        dpedek=-(1.-PE).*bata.* mid./(1.+ mid);
        dpedrk=-1*dpedek;
        drkdx=[0,0,0,0,0,0,1];
        df0dx(Var_num*k-Var_num+1:Var_num*k,1)=(dfdpe.*dpedek)'*dekdx{k}+dfdpe'*dpedrk*drkdx;
        dfdx(Var_num*k-Var_num+1:Var_num*k,1)=(dpedek'*dekdx{k})+sum(dpedrk*drkdx);
    end
    %MMA solve
    xval=xy00;
    xold=xy00;
    f0val =Comp*scale;
    df0dx=df0dx*scale;
    df0dx2=0*df0dx;
    fval=sum(PE)/(nelx*nely*nelz*volfrac)-1;
    dfdx=dfdx/(nelx*nely*nelz*volfrac);
    dfdx2=0*dfdx;
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,ss,low,upp] = mmasub(m,nn,Loop,xval,xmin,xmax,xold1,xold2,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,b,c,d,alctr);
     % Delete the minimal bars to accelerate the process
    [xy00,xold,xold1,xmin,xmax,low,upp,N,nn,xmma]=deleteBar(xy00,xold,xold1,xmin,xmax,low,upp,N,Var_num,xmma);
    xold2 = xold1;
    xold1 = xy00;
       xy00=xmma;
    change=max(max(abs( xy00-xold)));
    disp([' It.: ' sprintf('%4i\t',Loop) ' Obj.: ' sprintf('%6.3f\t',Comp) ' Vol.: ' ...
        sprintf('%6.4f\t',mean(PE)) 'ch.:' sprintf('%6.4f\t',change)]);
    Loop = Loop + 1;
end
end
%===Geometrically nonlinear FEA + Residual force sensitivity with the element density===
function [U,GKF,dRdpe]=GNFEA(ndof,nele,Freedofs,IterMax,ForceVector,Nodes,Elements,edofMat,E,nu,x,p,uselessNode,pmin)
U=zeros(ndof,1); loop = 0; % Iteration number
du=zeros(ndof,1);  % Displacement increment
ResidualForceMax=max(abs(ForceVector(Freedofs)));
D0=LinearElasticD(E,nu);  % Given constitutive tensor
TOL=1;% Newtonian iteration tolerance
while loop<IterMax  && ResidualForceMax>TOL
    Dimension=3;
    Fint=zeros(ndof,1);% Initialize the internal force vector
    GaussCoordinate=[-0.57735026918963D0, 0.57735026918963D0];
    GaussWeight=[1.00000000000000D0, 1.00000000000000D0];
    iK = zeros(24*24,nele);
    jK = zeros(24*24,nele);
    sK = zeros(24*24,nele);
    for i=1:nele
        D=(x(i)^(p)*(1-pmin)+pmin)*D0;
        ElementNodeCoordinate=Nodes(Elements(i,:),:); % Element node coordinate matrix
        edof=edofMat(i,:);
        ElementDisplacement0=U(edof);
        ElementDisplacement=reshape(ElementDisplacement0,Dimension,8);
        % Calculate all Gaussian points in each element
        tmp=zeros(24);
        for LX=1:2, for LY=1:2, for LZ=1:2
                    E1=GaussCoordinate(LX); E2=GaussCoordinate(LY); E3=GaussCoordinate(LZ);
                    [dNdx, JacobiDET] = ShapeFunction([E1 E2 E3], ElementNodeCoordinate);
                    FAC=GaussWeight(LX)*GaussWeight(LY)*GaussWeight(LZ)*JacobiDET;
                    F=ElementDisplacement*dNdx' + eye(3);  % Calculate the deformation gradient
                    StrainMatrix=0.5*(F'*F-eye(3));   % Calculate Lagrangian strain
                    Strain=[StrainMatrix(1,1) StrainMatrix(2,2) StrainMatrix(3,3) 2*StrainMatrix(1,2) 2*StrainMatrix(2,3) 2*StrainMatrix(1,3)]';
                    GaussPointStress=D*Strain;
                    BN=zeros(6,24);
                    BG=zeros(9,24);
                    for I=1:8
                        COL=(I-1)*3+1:(I-1)*3+3;
                        BN(:,COL)=[dNdx(1,I)*F(1,1) dNdx(1,I)*F(2,1) dNdx(1,I)*F(3,1);
                            dNdx(2,I)*F(1,2) dNdx(2,I)*F(2,2) dNdx(2,I)*F(3,2);
                            dNdx(3,I)*F(1,3) dNdx(3,I)*F(2,3) dNdx(3,I)*F(3,3);
                            dNdx(1,I)*F(1,2)+dNdx(2,I)*F(1,1) dNdx(1,I)*F(2,2)+dNdx(2,I)*F(2,1) dNdx(1,I)*F(3,2)+dNdx(2,I)*F(3,1);
                            dNdx(2,I)*F(1,3)+dNdx(3,I)*F(1,2) dNdx(2,I)*F(2,3)+dNdx(3,I)*F(2,2) dNdx(2,I)*F(3,3)+dNdx(3,I)*F(3,2);
                            dNdx(1,I)*F(1,3)+dNdx(3,I)*F(1,1) dNdx(1,I)*F(2,3)+dNdx(3,I)*F(2,1) dNdx(1,I)*F(3,3)+dNdx(3,I)*F(3,1)];
                        BG(:,COL)=[dNdx(1,I) 0         0;
                            dNdx(2,I) 0         0;
                            dNdx(3,I) 0         0;
                            0         dNdx(1,I) 0;
                            0         dNdx(2,I) 0;
                            0         dNdx(3,I) 0;
                            0         0         dNdx(1,I);
                            0         0         dNdx(2,I);
                            0         0         dNdx(3,I)];
                    end
                    Fint(edof) = Fint(edof) + FAC*BN'*GaussPointStress; % Internal force
                    SIG=[GaussPointStress(1) GaussPointStress(4) GaussPointStress(6);
                        GaussPointStress(4) GaussPointStress(2) GaussPointStress(5);
                        GaussPointStress(6) GaussPointStress(5) GaussPointStress(3)];
                    SHEAD=zeros(9);
                    SHEAD(1:3,1:3)=SIG;  SHEAD(4:6,4:6)=SIG;  SHEAD(7:9,7:9)=SIG;
                    EKF = BN'*D*BN + BG'*SHEAD*BG; % Calculate element stiffness matrix
                    tmp=tmp+FAC*EKF;
        end; end; end
        [dof1, dof2] = meshgrid(edof);
        iK(:,i) = dof1(:);
        jK(:,i) = dof2(:);
        sK(:,i) = tmp(:);
    end
    GKF=sparse(iK,jK,sK,ndof,ndof);% Assemble global stiffness matrix
    ResidualForce = ForceVector-Fint;
    du(Freedofs,:) = GKF(Freedofs,Freedofs)\ResidualForce(Freedofs,:); % Calculate the increment of the displacement
    U = U + du;
    if  loop>0 && ~isempty(uselessNode)
        uselessDOF=[3*uselessNode-1;3*uselessNode];% Nodal DOFs belongs the low density elements
        Freedofs=setdiff(Freedofs,uselessDOF);  % Update the freedoms DOFs for the convergence
        ResidualForceMax=max(abs(ResidualForce(Freedofs)));
        disp([' N-R process.' sprintf('%4i\t:',loop) 'Residual: ' sprintf('%6.3f\t',full(ResidualForceMax)) ]);
    end
    loop = loop + 1;
end
iK=zeros(24,size(Elements,1));
jK=zeros(24,size(Elements,1));
sK=zeros(24,size(Elements,1));
for i=1:nele
    dDdx=(p*x(i)^(p-1)*(1-pmin))*D0;
    ElementNodeCoordinate=Nodes(Elements(i,:),:);
    edof=edofMat(i,:);
    ElementDisplacement0=U(edof);
    ElementDisplacement=reshape(ElementDisplacement0,Dimension,8);
    tmp=zeros(24,1);
    for LX=1:2, for LY=1:2, for LZ=1:2
                E1=GaussCoordinate(LX); E2=GaussCoordinate(LY); E3=GaussCoordinate(LZ);
                [dNdx, JacobiDET] = ShapeFunction([E1 E2 E3], ElementNodeCoordinate); 
                FAC=GaussWeight(LX)*GaussWeight(LY)*GaussWeight(LZ)*JacobiDET;
                F=ElementDisplacement*dNdx' + eye(3); 
                StrainMatrix=0.5*(F'*F-eye(3));
                Strain=[StrainMatrix(1,1) StrainMatrix(2,2) StrainMatrix(3,3) 2*StrainMatrix(1,2) 2*StrainMatrix(2,3) 2*StrainMatrix(1,3)]';
                BN=zeros(6,24);  
                for I=1:8  
                    COL=(I-1)*3+1:(I-1)*3+3;
                    BN(:,COL)=[dNdx(1,I)*F(1,1) dNdx(1,I)*F(2,1) dNdx(1,I)*F(3,1);
                        dNdx(2,I)*F(1,2) dNdx(2,I)*F(2,2) dNdx(2,I)*F(3,2);
                        dNdx(3,I)*F(1,3) dNdx(3,I)*F(2,3) dNdx(3,I)*F(3,3);
                        dNdx(1,I)*F(1,2)+dNdx(2,I)*F(1,1) dNdx(1,I)*F(2,2)+dNdx(2,I)*F(2,1) dNdx(1,I)*F(3,2)+dNdx(2,I)*F(3,1);
                        dNdx(2,I)*F(1,3)+dNdx(3,I)*F(1,2) dNdx(2,I)*F(2,3)+dNdx(3,I)*F(2,2) dNdx(2,I)*F(3,3)+dNdx(3,I)*F(3,2);
                        dNdx(1,I)*F(1,3)+dNdx(3,I)*F(1,1) dNdx(1,I)*F(2,3)+dNdx(3,I)*F(2,1) dNdx(1,I)*F(3,3)+dNdx(3,I)*F(3,1)];
                end
                tmp=tmp+ FAC*BN'*dDdx*Strain;
%                 dRdpe(edof,i) = dRdpe(edof,i) + FAC*BN'*dDdx*Strain;
    end; end; end
    iK(:,i)=edof;
    jK(:,i)=i;
    sK(:,i)=tmp;
end
dRdpe=sparse(iK,jK,sK,ndof,size(Elements,1));
end
%=======Obtain elements and nodes duringdeformation large displacement=======
function [Nodes3Ele]=NodeSurroundingElement(Nodes,Eles)
Nodes3Ele=sparse(size(Nodes,1),size(Eles,1));
for i=1:size(Nodes,1)
    [hang,~]=find(Eles==i);
    Nodes3Ele(i,hang)=ones(1,size(hang,1));
end
end
function [uselessNode]=NodeSurroundWithVoidElement(Nodes2Ele, xPhy,xmin)
xPhyMatrix=repmat(xPhy',size(Nodes2Ele,1),1);
a=Nodes2Ele.*xPhyMatrix;
[max_a,~]=max(a,[],2);
[uselessNode]=find(max_a<xmin);% Find node numbers surrounded by low densities
end
% Calculate the derivation of shape function and Jacobi matrix
function [dNdx, JacobiDET] = ShapeFunction(GaussPoint, ElementNode)
ParentNodes=[-1  1  1 -1 -1  1  1 -1;
    -1 -1  1  1 -1 -1  1  1;
    -1 -1 -1 -1  1  1  1  1]; %The coordinates for the reference element
ParentNDerivative=zeros(3,8);
for I=1:8
    XPoint = ParentNodes(1,I);
    YPoint = ParentNodes(2,I);
    ZPoint = ParentNodes(3,I);
    ShapePart = [1+GaussPoint(1)*XPoint 1+GaussPoint(2)*YPoint 1+GaussPoint(3)*ZPoint];  
    ParentNDerivative(1,I) = 0.125*XPoint*ShapePart(2)*ShapePart(3);
    ParentNDerivative(2,I) = 0.125*YPoint*ShapePart(1)*ShapePart(3);
    ParentNDerivative(3,I) = 0.125*ZPoint*ShapePart(1)*ShapePart(2);
end
Jacobi = ParentNDerivative*ElementNode;% Calculate Jacobi matrix
JacobiDET = det(Jacobi);
JacobiINV=inv(Jacobi);
dNdx=JacobiINV*ParentNDerivative;% Calculate the derivation of shape function
end
%==== Forming modified dek and dekdx for each MMB ====
function [dek,dekdx]=DEK(xk,xe,ye,ze)
x1=xk(1); y1=xk(2);z1=xk(3);x2=xk(4); y2=xk(5);  z2=xk(6);
Smin=1E-9;
lk=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
a1=xe-x1;  a2=ye-y1;  a3=ze-z1;
b1=xe-x2;  b2=ye-y2;  b3=ze-z2;
n=length(xe); %n is the number of the elements
dek=zeros(n,1); dekdx=zeros(n,7);   num_all=[1:n]';
num1=find((a1.*(x2-x1)+a2.*(y2-y1)+a3.*(z2-z1))<=0);
dek(num1)=max(Smin,sqrt(a1(num1).^2+ a2(num1).^2+ a3(num1).^2));
dekdx(num1,1)=-a1(num1)./dek(num1);
dekdx(num1,2)=-a2(num1)./dek(num1);
dekdx(num1,3)=-a3(num1)./dek(num1);
num2=find((b1.*(x1-x2)+b2.*(y1-y2)+b3.*(z1-z2))<=0);
dek(num2)=max(Smin,sqrt(b1(num2).^2+ b2(num2).^2+ b3(num2).^2));
dekdx(num2,4)=-b1(num2)./dek(num2);
dekdx(num2,5)=-b2(num2)./dek(num2);
dekdx(num2,6)=-b3(num2)./dek(num2);
num3=setdiff(num_all,[num1;num2]);
A=b2.*a3-a2.*b3; B=a1.*b3-b1.*a3; C=b1.*a2-a1.*b2;
sqrtABC=sqrt(A.^2+B.^2+C.^2);
dek(num3)=sqrtABC(num3)./lk; fm=max(Smin,(lk^3).*sqrtABC);
dekdx(num3,1)=((-b3(num3).*B(num3)+b2(num3).*C(num3))*lk^2+(sqrtABC(num3).^2).*(x2-x1))./fm(num3);
dekdx(num3,2)=(( b3(num3).*A(num3)-b1(num3).*C(num3))*lk^2+(sqrtABC(num3).^2).*(y2-y1))./fm(num3);
dekdx(num3,3)=((-b2(num3).*A(num3)+b1(num3).*B(num3))*lk^2+(sqrtABC(num3).^2).*(z2-z1))./fm(num3);
dekdx(num3,4)=(( a3(num3).*B(num3)-a2(num3).*C(num3))*lk^2-(sqrtABC(num3).^2).*(x2-x1))./fm(num3);
dekdx(num3,5)=((-a3(num3).*A(num3)+a1(num3).*C(num3))*lk^2-(sqrtABC(num3).^2).*(y2-y1))./fm(num3);
dekdx(num3,6)=(( a2(num3).*A(num3)-a1(num3).*B(num3))*lk^2-(sqrtABC(num3).^2).*(z2-z1))./fm(num3);
end

