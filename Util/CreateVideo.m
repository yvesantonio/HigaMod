workingDir = '/Users/YvesBarbosa/Documents/PhD Courses/ANEDPII/LABORATORI/LAB6';

for ii = 1:119
    fileName = ['Sol_at_',num2str(ii),'.pdf'];
    imageNames{ii} = fileName;
    eps2xxx(fullfile(workingDir,imageNames{ii}),{'png'});
    fileName = ['Sol_at_',num2str(ii),'.png'];
    imageNames{ii} = fileName;
end

outputVideo = VideoWriter(fullfile(workingDir,'shuttle_out.avi'));
outputVideo.FrameRate = 5;
open(outputVideo)

for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir,imageNames{ii}));
   writeVideo(outputVideo,img)
end

close(outputVideo)