%% Import Daten von ASCII
% to move over from using the mapping toolbox this needs (obviously) a
% major rewrite

filename = 'D:\02_Documents\04_Projekte\Couchtisch\Matlab\Daten\10m.asc';
delimiterIn = ' ';
headerlinesIn = 6;

data = importdata(filename, delimiterIn, headerlinesIn);