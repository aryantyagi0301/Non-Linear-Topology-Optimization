% Geometric Nonlinear Topology Optimization base on SIMP
function SIMP_3D_GNTO(DL,DW,DH,nelx,nely,nelz,volfrac,rmin)
p=1;
dp=0.1;   % The increment of the penal coefficient
E=3e9;    % Young's modulus
nu=0.4; % Poisson ratio
xthreshold=0.01; % Low density threshold
force=200 ;  % Magnitude of external force
IterMax=16;  % Maximum number of Newton-Raphson iterations
pmin=1e-4;
parent_dir_name ='GNTO results';
if  exist(parent_dir_name,"dir") 
rmdir(parent_dir_name, 's') % Delete the existed file
end
mkdir(parent_dir_name);
% Calculate nodal coordinates of all elements
EL=DL/nelx; EW=DW/nely; EH=DH/nelz;
EA=EL*EW*EH; % Element volume
[ey1,ex1,ez1] = meshgrid(EW * (0 : nely),EL * (0: nelx),EH * (0: nelz));
Nodes=[ex1(:) ey1(:) ez1(:)];
% DOFs for the loading force for the slender beam
loadnid=(nelx+1)*(nely+1)*nelz+(nelx/2+1):(nelx+1):(nelx+1)*(nely+1)*(nelz+1);
loaddof = 3*loadnid;
% DOFs for the fixed ends
fixednidL =(nelx+1)*(nely+1)*nelz/2+1:(nelx+1):(nelx+1)*(nely+1)*(nelz/2+1);
fixednidR=fixednidL+nelx;
fixednid =[fixednidL fixednidR];
Fixeddofs = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2];
% Prepare for geometrically nonlinear FEA
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
Alldofs=1:ndof;
Freedofs=setdiff(Alldofs,Fixeddofs);
numSteps = 10; % Number of load steps (adjust based on nonlinearity)
ForceStep = sparse(loaddof,1,force/numSteps,ndof,1); % Incremental force 
% ForceVector = sparse(loaddof,1,force,ndof,1);   % External force vector
nodegrd = reshape(1:(nelx+1)*(nely+1),nelx+1,nely+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:);
edofMat=repmat(edofVec,1,24)+repmat([-2 -1 0 1 2 3 3*nelx+[4 5 6 1 2 3] ...
    3*(nely+1)*(nelx+1)+[-2 -1 0 1 2 3 3*nelx+[4 5 6 1 2 3]]] ,nele,1);
Eles= repmat(nodeids(:),1,8)+repmat([0 1 nelx+[2 1] ...
    (nely+1)*(nelx+1)+[0 1 nelx+[2 1]]],nele,1); % Nodal number corresponding to each element
% initial density values
x=volfrac*ones(nele,1);
[H,Hs] =prepareFilter(nele,nelx,nely,nelz,EL,rmin);% Prepare for the filter
[Nodes2Ele]=NodeSurroundingElement(Nodes,Eles); % Calculate the center and nodal coordinates of all elements
% Maximum and minimum values for the design variables
xmin    = 0.001*ones(nele,1);
xmax    = ones(nele,1);
% MMA parameters
m=1;   nn=nele;
low = xmin;
upp = xmax;
c=1000*ones(m,1);
a0=1;
[d,a]=deal(zeros(m,1));
xy00=x;
[xold1,xold2]=deal(xy00);
Loop=1; change=1; maxloop=201; Lamada=zeros(ndof,1);
while change>0.00001 && Loop<maxloop
    % The update for the penal coefficient
    if(Loop~=0)
        p=p+dp;
        fprintf(['\t\t *** Successful loop converge ***   p= ',num2str(p,3),' \n']);
        if p>2.998
            p=3;
        end
    end
    % Show the figure for the optimized result
    clf; showDensity(Nodes,Eles,xy00(:)); 
    FileName=[parent_dir_name,'\Fig_',int2str(Loop), '.png'];     saveas(gcf,FileName);
    % Find nodes surrounded by low-density elements
    [uselessNode]=NodeSurroundWithVoidElement(Nodes2Ele, xy00,xthreshold);
    % Geometrically nonlinear FEA + Residual force sensitivity with the element density
%     profile clear
%     profile on
    U = zeros(ndof,1); 

for step = 1:numSteps
    % incremental load
    CurrentForce = ForceStep * step;
    
    [U,GKF,dRdpe] = GNFEA(ndof,nele,Freedofs,IterMax,CurrentForce,Nodes,Eles,edofMat,E,nu,xy00,p,uselessNode,pmin,U);
    
    fprintf(' Load step %d/%d completed \n', step, numSteps);
end
% [U,GKF,dRdpe]=GNFEA(ndof,nele,Freedofs,IterMax,ForceVector,Nodes,Eles,edofMat,E,nu,xy00,p,uselessNode,pmin);
%     profile viewer
    % Lamada(Freedofs,:)=GKF(Freedofs,Freedofs)\ForceVector(Freedofs,:);%K*lamada=F
    Lamada(Freedofs,:)=GKF(Freedofs,Freedofs)\CurrentForce(Freedofs,:);
    dcdpe=full(-Lamada'*dRdpe)';
    % objective function, the sensitivity and the filter
    % Comp=ForceVector'*U; % objective
    Comp=CurrentForce'*U;
    df0dx= H*(xy00(:).*dcdpe(:))./Hs./max(1e-3,xy00(:));
    df0dx(:)=df0dx(:)./max(abs(df0dx(:)));
    fval=sum(xy00)-(nelx*nely*nelz*volfrac);
    dfdx=reshape(EA*ones(1,nelx*nely*nelz),1,nelx*nely*nelz);
    dfdx=H*(xy00(:).*dfdx(:))./Hs./xy00(:);
    %MMA
    xval=xy00;
%     xold=xy00;
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,nn,Loop,xval,xmin,xmax,xold1,xold2,...
        Comp,df0dx,fval,dfdx,low,upp,a0,a,c,d);
    xold2 = xold1;
    xold1 = xy00;
    xy00 = xmma;
    change=max(max(abs( xy00-xold1)));
    disp([' It.: ' sprintf('%4i\t',Loop) ' Obj.: ' sprintf('%6.3f\t',Comp) ' Vol.: ' ...
        sprintf('%6.4f\t',mean(xy00)) 'ch.:' sprintf('%6.4f\t',change)]); 
    Loop = Loop + 1; 
end
end
%=== Geometrically nonlinear FEA + Residual force sensitivity with the element density===
% function [U,GKF,dRdpe]=GNFEA(ndof,nele,Freedofs,IterMax,ForceVector,Nodes,Elements,edofMat,E,nu,x,p,uselessNode,pmin)
function [U,GKF,dRdpe] = GNFEA(ndof,nele,Freedofs,IterMax,ForceVector,Nodes,Elements,edofMat,E,nu,x,p,uselessNode,pmin,U)
% U=zeros(ndof,1);
loop = 0; % Iteration number
du=zeros(ndof,1);  % Displacement increment
ResidualForceMax=max(abs(ForceVector(Freedofs)));
ResidualForce = zeros(size(ForceVector));
D0=LinearElasticD(E,nu);  % Given constitutive tensor
TOL=1e-3;% Newtonian iteration tolerance
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
        ElementNodeCoordinate=Nodes(Elements(i,:),:);
        edof=edofMat(i,:);
        ElementDisplacement0=U(edof);
        ElementDisplacement=reshape(ElementDisplacement0,Dimension,8);
        % Calculate all Gaussian points in each element
        tmp=zeros(24);
        for LX=1:2, for LY=1:2, for LZ=1:2
                    E1=GaussCoordinate(LX); E2=GaussCoordinate(LY); E3=GaussCoordinate(LZ);
                    [dNdx, JacobiDET] = ShapeFunction([E1 E2 E3], ElementNodeCoordinate);
                    FAC=GaussWeight(LX)*GaussWeight(LY)*GaussWeight(LZ)*JacobiDET;
                    F=ElementDisplacement*dNdx' + eye(3);  %Calculate the deformation gradient
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
    % ResidualForce = ForceVector-Fint;
    ResidualForce = ResidualForce + ForceVector - Fint;
    du(Freedofs,:) = GKF(Freedofs,Freedofs)\ResidualForce(Freedofs,:);
    U = U + du;
    if  loop>0
        uselessDOF=[3*uselessNode-1;3*uselessNode];% Low density elements nodal DOFs
        Freedofs=setdiff(Freedofs,uselessDOF);  % Update freedoms DOFs for convergence
        ResidualForceMax=max(abs(ResidualForce(Freedofs)));
        % disp([' N-R process.' sprintf('%4i\t:',loop) 'Residual: ' sprintf('%6.3f\t',full(ResidualForceMax)) ]);
        % disp([' Load step: ' sprintf('%2i/%2i\t',step,numSteps) ' N-R Iteration: ' sprintf('%4i\t',loop) ' Residual: ' sprintf('%6.3f\t',full(ResidualForceMax)) ]);

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
%===Obtain elements and nodes during deformation large displacement===
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
[uselessNode]=find(max_a<xmin); % Find node numbers surrounded by low densities
end
% Calculate the derivation of shape function and Jacobi matrix
function [dNdx, JacobiDET] = ShapeFunction(GaussPoint, ElementNode)
ParentNodes=[-1  1  1 -1 -1  1  1 -1;
    -1 -1  1  1 -1 -1  1  1;
    -1 -1 -1 -1  1  1  1  1];
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
Jacobi = ParentNDerivative*ElementNode; % Calculate Jacobi matrix
JacobiDET = det(Jacobi);
% JacobiINV=inv(Jacobi);
% dNdx=JacobiINV*ParentNDerivative;
dNdx=Jacobi\ParentNDerivative;
end
