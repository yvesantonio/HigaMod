%% Flow 1
 
% %%%%%%%%%%%%%%%%%%%
% % FLOW PROPERTIES %
% %%%%%%%%%%%%%%%%%%%
% 
% R         = 1;
% numbNodes = 100;
% Y         = linspace(0,R,numbNodes);
% T         = 0.25;
% w         = 2 * pi / T;
% rho       = 1080; 
% mu        = 20;
% p0        = 0;
% pVect     = [0.78 1.32 -0.74];
% phiVect   = [-0.01 -1.45 -0.46];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % TIME SIMULATION PARAMETERS %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% numbIt = 100;
% dT     = 1/numbIt;
% time   = 0:dT:1;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SURFACE DISCRETIZATION PARAMETERS %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% numbTh = 100;
% radius = 2 * Y - R;
% theta  = linspace(0,pi,numbTh);
% [RR,TT] = meshgrid(radius,theta);
% XX = RR .* cos(TT);
% YY = RR .* sin(TT);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % COMPUTE FLOW FOR EACH TIME STEP %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% U  = zeros(numbNodes,numbIt);
% sU = zeros(numbNodes,numbTh,numbIt);
% for ii = 1:numbIt
%     U(:,ii) = 5e4 * WomersleyFlow(Y,time(ii),rho,mu,R,p0,pVect,w,phiVect);
%     for jj = 1:numbTh
%         sU(:,jj,ii) = U(:,ii);
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% % PLOT PROFILE IN TIME %
% %%%%%%%%%%%%%%%%%%%%%%%%
% 
% % for jj = 1:numbIt
% %     
% %     figure(1)
% %     plot(radius,U(:,jj),'-b');
% %     grid on;
% %     
% %     xMin = -1;
% %     xMax = 1;
% %     yMin = min(min(U));
% %     yMax = -min(min(U));
% %     
% %     axis([xMin xMax yMin yMax])
% %     
% %     F(jj) = getframe(gcf) ;
% %     drawnow
% %     
% % end
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % CREATE VIDEO WITH 1 FPS %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % writerObj = VideoWriter('myVideo.avi');
% % writerObj.FrameRate = 10;
% % 
% % 
% % % set the seconds per image
% % % open the video writer
% % 
% % open(writerObj);
% % 
% % % write the frames to the video
% % 
% % for i=1:length(F)
% %     
% %     % convert the image to a frame
% %     
% %     frame = F(i) ;    
% %     writeVideo(writerObj,frame);
% %     
% % end
% % 
% % % close the writer object
% % 
% % close(writerObj);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% % PLOT SURFACE IN TIME %
% %%%%%%%%%%%%%%%%%%%%%%%%
% 
% for jj = 1:numbIt
%     
%     figure(2)  
%     surf(XX',YY',sU(:,:,jj),sU(:,:,jj));
%     colorbar;
%     
%     xMin = -1;
%     xMax = 1;
%     yMin = -1;
%     yMax = 1;
%     zMin = min(min(U));
%     zMax = -min(min(U));
%     cMin = min(min(U));
%     cMax = -min(min(U));
%     
%     xlim([xMin xMax])
%     ylim([yMin yMax])
%     zlim([zMin zMax])
%     caxis([cMin cMax])
%     
%     axis off
%     grid off
%     
%     F(jj) = getframe(gcf) ;
%     drawnow
%     
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % CREATE VIDEO WITH 1 FPS %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% writerObj = VideoWriter('myVideo.avi');
% writerObj.FrameRate = 10;
% 
% 
% % set the seconds per image
% % open the video writer
% 
% open(writerObj);
% 
% % write the frames to the video
% 
% for i=1:length(F)
%     
%     % convert the image to a frame
%     
%     frame = F(i) ;    
%     writeVideo(writerObj,frame);
%     
% end
% 
% % close the writer object
% 
% close(writerObj);

%% Flow 2


Wn = table2array(pmod2);
FFTn = fft(Wn);

Idx = linspace(1,length(FFTn),length(FFTn))';
Re  = real(FFTn);
Im  = imag(FFTn);
A   = [Idx Re Im];

fileID = fopen('FFTBloodFlow.txt','w');
fprintf(fileID,'%d %e %e\n',A');
fclose(fileID);


MODn = abs(FFTn);
ANGn = angle(FFTn)/pi;

figure
plot (MODn,'-o');
grid on

CUTn = FFTn(1:50);
Rec  = ifft(FFTn);

figure
plot(Rec)

%%%%%%%%%%%%%%%%%%%
% FLOW PROPERTIES %
%%%%%%%%%%%%%%%%%%%

R         = 1;
numbNodes = 100;
Y         = linspace(0,R,numbNodes);
T         = 1;
w         = 2 * pi / T;
rho       = 1080; 
mu        = 3.5;
p0        = 0;
% phiVect   = ANGn(1:300);
% pVect     = FFTn(1:300); 

phiVect = [-4.027 -6.509 -1.913 -1.461 -0.074];
pVect   = [0.29244 -0.5908 0.2726 0.198 0.1124];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME SIMULATION PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numbIt = 100;
dT     = 1/numbIt;
time   = 0:dT:1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SURFACE DISCRETIZATION PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numbTh = 100;
radius = 2 * Y - R;
theta  = linspace(0,pi,numbTh);
[RR,TT] = meshgrid(radius,theta);
XX = RR .* cos(TT);
YY = RR .* sin(TT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE FLOW FOR EACH TIME STEP %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U  = zeros(numbNodes,numbIt);
sU = zeros(numbNodes,numbTh,numbIt);
for ii = 1:numbIt
    U(:,ii) = 5e4 * abs(WomersleyFlow(Y,time(ii),rho,mu,R,p0,pVect,w,phiVect));
    for jj = 1:numbTh
        sU(:,jj,ii) = U(:,ii);
    end
end

figure
plot(time(2:end),U(50,:))

%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT PROFILE IN TIME %
%%%%%%%%%%%%%%%%%%%%%%%%

% for jj = 1:numbIt
%     
%     figure(1)
%     plot(radius,U(:,jj),'-b');
%     grid on;
%     
%     xMin = -1;
%     xMax = 1;
%     yMin = min(min(U));
%     yMax = -min(min(U));
%     
%     axis([xMin xMax yMin yMax])
%     
%     F(jj) = getframe(gcf) ;
%     drawnow
%     
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % CREATE VIDEO WITH 1 FPS %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% writerObj = VideoWriter('myVideo.avi');
% writerObj.FrameRate = 10;
% 
% 
% % set the seconds per image
% % open the video writer
% 
% open(writerObj);
% 
% % write the frames to the video
% 
% for i=1:length(F)
%     
%     % convert the image to a frame
%     
%     frame = F(i) ;    
%     writeVideo(writerObj,frame);
%     
% end
% 
% % close the writer object
% 
% close(writerObj);

%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SURFACE IN TIME %
%%%%%%%%%%%%%%%%%%%%%%%%

for jj = 1:numbIt
    
    figure(2)  
    s = surf(XX',YY',sU(:,:,jj),sU(:,:,jj));
    colorbar;
    
    xMin = -1;
    xMax = 1;
    yMin = -1;
    yMax = 1;
    zMin = min(min(U));
    zMax = -min(min(U));
    cMin = min(min(U));
    cMax = -min(min(U));
    
    xlim([xMin xMax])
    ylim([yMin yMax])
    zlim([zMin zMax])
    caxis([cMin cMax])
    
    axis off
    grid off
    
    F(jj) = getframe(gcf) ;
    drawnow
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE VIDEO WITH 1 FPS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

writerObj = VideoWriter('myVideo.avi');
writerObj.FrameRate = 10;


% set the seconds per image
% open the video writer

open(writerObj);

% write the frames to the video

for i=1:length(F)
    
    % convert the image to a frame
    
    frame = F(i) ;    
    writeVideo(writerObj,frame);
    
end

% close the writer object

close(writerObj);