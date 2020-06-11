% schistosome readme.m
% datastore @ octerian: /mnt/spaldingdata/nate/octerineDataStore/planarian/
%% date: April 7, 2020
% working through continuity equation
% sample the data over many frames 
% data - featurevector list in frame of average velocity
%% background update
% cat and sort the frames to better approximate the background
% the min value is 'very' zero.  The nth (n=20) is a small number and does
% not contain the noise of 'very' zero.
%%
project.save('schistosome');