clear
close all

date = 20220617;
file = 924;

DataFolder = ['Z:\Data\WG3\' num2str(date)];

image = imread(fullfile(DataFolder, [num2str(file) '.tif']));
figure(1);
imagesc(image)
axis image
r_list = [195 100 30];
for i = 1:3
viscircles([250,250],r_list(i),'LineWidth',0.1)
end
line([0 500], [250 250], 'Color', 'k');
line([250 250], [0 500], 'Color', 'k');
[XYG]=ginput(1);
dX = (XYG(1)-250)*0.329;
dY = (250-XYG(2))*0.329;