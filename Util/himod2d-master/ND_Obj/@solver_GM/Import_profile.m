function []=Import_profile(obj,num)
filetoread='profile_at_interface';
data = importdata(filetoread);
obj.profile{num}=data(:,2:3);
end
