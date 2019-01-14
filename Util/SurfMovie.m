function SurfMovie(II)
num =  size(II,1);
grid = size(II,2);
for i=1:16
    J(:,:)=double(II(i,:,:));
    surf(J,'FaceColor','red','EdgeColor','none')
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    axis([0 grid 0 grid 0 max(max(max(II)))])
    axis off
    title(['Oscillation phase angle, \phi= ',num2str(i*360/num)])
    camlight left;
    lighting phong
    pause(0.5)
end
end