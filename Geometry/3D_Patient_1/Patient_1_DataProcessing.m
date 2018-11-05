%%%%%%%%%%%%%%%%%%%%%%%%%%
% MANAGE PATIENT GEOMETRY
% INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

% READ .CSV DATA FILE

filename = [pwd,'\morphology\centerlines.csv'];
DATA = csvread(filename);

minX = min(DATA(:,1));
minY = min(DATA(:,2));
minZ = min(DATA(:,3));
minR = min(DATA(:,4));

maxX = max(DATA(:,1));
maxY = max(DATA(:,2));
maxZ = max(DATA(:,3));
maxR = max(DATA(:,4));

% COMPUTE JUMPS IN CENTERLINE INFORMATION
% Compute the jumps in the dataset. If the jump is gretater than 0.1 mm,
% the precision of the centerline data, then we have another centerline.

len = size(DATA,1);

X     = DATA(1:end-1,1:3);
Xnext = DATA(2:end,1:3);
jump  = (Xnext - X).^2;
jump  = sum(jump,2);
jump  = jump.^(1/2);

I = find(jump > 0.5);

% EXTRACT CENTERLINES
% The dataset contains a list of points of the various centerlines in the
% arterial network. We first have to extract each individual segment, to
% separate the simulation domains.

lineStruct = cell(size(I,1)+1,1);

for ii = 1:size(I,1)+1
    if(ii == 1)
        lineStruct{1} = DATA(1:I(1),:);
    elseif(ii == size(I,1)+1)
        lineStruct{end} = DATA(I(end) + 2:end,:);
    else
        lineStruct{ii} = DATA(I(ii-1)+2:I(ii),:);
    end
end

% PLOT OBTAINED CENTERLINES

figure;
for ii = 1:size(lineStruct,1)
    plot3(lineStruct{ii}(:,1),lineStruct{ii}(:,2),lineStruct{ii}(:,3),'.');
    hold on
end
grid on

% IDENTIFY CROSS POINTS 
% We compute the intersection points and its plane properties in each
% bifurcation of the arterial network.

interpStruct = cell(size(lineStruct,1),4);

for ii = 1:size(lineStruct,1)
    
    param = linspace(0,1,length(lineStruct{ii}(:,1)));
    
    interpX = griddedInterpolant(param,lineStruct{ii}(:,1));
    interpY = griddedInterpolant(param,lineStruct{ii}(:,2));
    interpZ = griddedInterpolant(param,lineStruct{ii}(:,3));
    interpR = griddedInterpolant(param,lineStruct{ii}(:,4));
    
    interpStruct{ii,1} = interpX;
    interpStruct{ii,2} = interpY;
    interpStruct{ii,3} = interpZ;
    interpStruct{ii,4} = interpR;
    
end

% Biff  = cell(size(lineStruct,1));
Coord = cell(size(lineStruct,1),1);

for ii = 1:size(lineStruct,1)    
    N = 1e3;
    param = linspace(0,1,N);

    X1 = interpStruct{ii,1}(param);
    Y1 = interpStruct{ii,2}(param);
    Z1 = interpStruct{ii,3}(param);
    
    Coord{ii} = [X1',Y1',Z1'];
end

% EXTRACT THE ARTERY BRANCHES
% This step is necessary since all the branches have the same starting
% point and we want to separate each individual domain to simulate.

tic;

InterMapX  = cell(size(lineStruct,1));
InterMapY  = cell(size(lineStruct,1));

for ii = 1:size(lineStruct,1)
    for jj = 1:size(lineStruct,1)
        
        if (jj > ii)
            
            X1 = Coord{ii}(:,1);
            Y1 = Coord{ii}(:,2);
            Z1 = Coord{ii}(:,3);
            X2 = Coord{jj}(:,1);
            Y2 = Coord{jj}(:,2);
            Z2 = Coord{jj}(:,3);
            
            [X0,Z0] = intersections(X1,Z1,X2,Z2,0);
            InterMapX{ii,jj}(:,1) = X0;
            InterMapX{ii,jj}(:,2) = Z0;
            
            [Y0,Z0] = intersections(Y1,Z1,Y2,Z2,0);
            InterMapY{ii,jj}(:,1) = Y0;
            InterMapY{ii,jj}(:,2) = Z0;
            
            Xaux = InterMapX{ii,jj};
            Yaux = InterMapY{ii,jj};
            
            for kk = 1:10
                
                dif  = ([X0(end-kk+1),Z0(end-kk+1)] - [X0(end-kk),Z0(end-kk)]).^2;
                dif  = sum(dif,2);
                distance  = dif.^(1/2);
                
                if(distance < 0.05)
                    InterPtsX{ii,jj} = [X0(end-kk+1),Z0(end-kk+1)];
                    break
                end
                
                dif  = ([Y0(end-kk+1),Z0(end-kk+1)] - [Y0(end-kk),Z0(end-kk)]).^2;
                dif  = sum(dif,2);
                distance  = dif.^(1/2);
                
                if(distance < 0.05)
                    InterPtsY{ii,jj} = [Y0(end-kk+1),Z0(end-kk+1)];
                    break
                end
                
            end

            figure;
            plot(X1,Z1,X2,Z2); hold on
            plot(InterMapX{ii,jj}(end-10:end,1),InterMapX{ii,jj}(end-10:end,2),'k.')
            plot(InterPtsX{ii,jj}(1),InterPtsX{ii,jj}(2),'ro')
            title(['XZ ',num2str([ii jj])])
            figure;
            plot(Y1,Z1,Y2,Z2); hold on
            plot(InterMapY{ii,jj}(end-10:end,1),InterMapY{ii,jj}(end-10:end,2),'k.')
            title(['YZ ',num2str([ii jj])])
        end
       
    end
end

time = toc;
disp(['Time Find Intersection : ',num2str(time),' [sec]'])

% Lag = [1 8000 1 1 1 1];
% Range = [N-1, N-Lag(2), N-Lag(3), N-Lag(4), N-Lag(5), N-Lag(6)];
% 
% for ii = 1:size(lineStruct,1)
%     aux = Coord{ii};
%     Coord{ii} = aux(Lag(ii):Lag(ii) + Range(ii),:);   
% end

% for ii = size(lineStruct,1):2
%     aux = Coord{ii};
%     Coord{ii} = aux(Biff{ii-1,ii}(end):end,:);
% end

% for ii = 2:size(lineStruct,1)
%     aux = Coord{ii};
%     Coord{ii} = aux(Biff{1,ii}(end):end,:);
% end
% 
% for ii = 3:size(lineStruct,1)
%     aux = Coord{ii};
%     Coord{ii} = aux(Biff{2,ii}(end):end,:);
% end

% PLOT SELECTED BRANCHES
% Plot the individual domain separated from the original network.

m = 3;
n = 2;
figure;
ColorSet = varycolor(10);

for ii = 1:size(lineStruct,1)
    subplot(m,n,ii)
    plot3(Coord{ii}(:,1),Coord{ii}(:,2),Coord{ii}(:,3),'.', 'Color', ColorSet(ii,:));
    hold on
    scatter3(Coord{ii}(1,1),Coord{ii}(1,2),Coord{ii}(1,3),...
             'MarkerEdgeColor','k','MarkerFaceColor',ColorSet(ii,:));
    scatter3(Coord{ii}(end,1),Coord{ii}(end,2),Coord{ii}(end,3),...
             'MarkerEdgeColor','k','MarkerFaceColor',ColorSet(ii,:));
    xlim([minX maxX])
    ylim([minY maxY])
    zlim([minZ maxZ])
    grid on
end
