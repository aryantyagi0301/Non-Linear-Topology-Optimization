% Show the figure for the optimized result
function showDensity(Nodes,Elements,rho)
HighDensityElementCount=0;
for i=1:1:size(Elements,1)
    if (rho(i) > 0.1)
        HighDensityElementCount=HighDensityElementCount+1;
    end
end
AllFaces=zeros(HighDensityElementCount*6,4); 
FaceColor=zeros(HighDensityElementCount*6,3);
faceID=0;
for i=1:1:size(Elements,1) 
    if (rho(i) > 0.1)
        faces_matrix= [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
        nodeID=Elements(i,:); 
        for j=1:6
            faceID=faceID+1;
            AllFaces(faceID,:)=nodeID(faces_matrix(j,:));
        end
        FaceColor(faceID-5:faceID,:)=0.2+0.8*(1-rho(i));
    end
end
patch('vertices', Nodes,'faces',AllFaces,'FaceVertexCData',FaceColor,'FaceColor','flat', 'facealpha',1);
view(3);
axis equal
axis off
hold on
end