%% Usage
% this code will produce images in the current folder,
% make sure to run this code in working folder (not the folder it
% lives in)

%% initialise mtex
folderpath = fileparts(mfilename('fullpath'));
if ~ismember(folderpath, strsplit(path, pathsep))
	startup
end
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outofPlane');
pfAnnotations = @(varargin) [];
setMTEXpref('pfAnnotations',pfAnnotations);

%% simulate bcc texture
CS = crystalSymmetry('cubic', [3.31 3.31 3.31],[90*degree 90*degree 90*degree]);
shape = crystalShape.cube(CS);

planes = [1 1 -2
	1 1 2];
directions = [1 1 1
	-1 -1 1];
ori = simluatePoleFigures(planes(1,:), directions(1,:), CS);
mtexColorMap WhiteJet
savePoleFigures(planes(1,:), directions(1,:));

ori = simluatePoleFigures(planes(2,:), directions(2,:), CS);
mtexColorMap WhiteJet
savePoleFigures(planes(2,:), directions(2,:));

ori = simluatePoleFigures(planes, directions, CS);
mtexColorMap WhiteJet
savePoleFigures(planes, directions);

%% simulate TMAZ bcc texture
CS = crystalSymmetry('cubic', [3.31 3.31 3.31],[90*degree 90*degree 90*degree]);
shape = crystalShape.cube(CS);
planes = [1 1 0];
directions = [0 0 1];
ori = simluatePoleFigures(planes, directions, CS);
mtexColorMap WhiteJet
savePoleFigures(planes, directions);

%% plot crystal shape
figure
plot(ori * shape)
view(0,90)
axis on
xlabel('x')
ylabel('y')
zlabel('z')
box on

%% simulate hcp texture
CS = crystalSymmetry('6/mmm',[2.95 2.95 4.68], [90 90 120]*degree);
shape = crystalShape.hex(CS);
planes = [1 0 -1 0];
directions = [1 1 -2 0];
ori = simluatePoleFigures(planes, directions, CS);
mtexColorMap WhiteJet
savePoleFigures(planes, directions);

%% CWZ bcc -> hcp
csParent = crystalSymmetry('cubic',[3.31 3.31 3.31],'mineral','Ti (beta)');
csChild = crystalSymmetry('6/mmm',[2.95 2.95 4.68], [90*degree 90*degree 120*degree],'mineral','Ti (alpha)');
beta2alpha = orientation.Burgers(csParent,csChild);

planes = [1 1 -2
	1 1 2];
directions = [1 1 1
	-1 -1 1];
oriParent = orientation.byMiller(planes, directions, csParent);
oriParents = oriParent.symmetrise;

oriChilds = oriParents * inv(beta2alpha);

plotPoleFigures(oriParents, csParent, 'OD', 'TD')
mtexColorMap WhiteJet
saveFigure(char("{11-2}[111] parent pole figures.png"));

plotPoleFigures(oriChilds, csChild, 'OD', 'TD')
mtexColorMap WhiteJet
saveFigure(char("{11-2}[111] child pole figures.png"));

%% TMAZ bcc -> hcp
csParent = crystalSymmetry('cubic',[3.31 3.31 3.31],'mineral','Ti (beta)');
csChild = crystalSymmetry('6/mmm',[2.95 2.95 4.68], [90*degree 90*degree 120*degree],'mineral','Ti (alpha)');
beta2alpha = orientation.Burgers(csParent,csChild);

planes = [1 1 0];
directions = [0 0 1];
oriParent = orientation.byMiller(planes, directions, csParent);
oriParents = oriParent.symmetrise;
oriChilds = oriParents * inv(beta2alpha);

plotPoleFigures(oriParents, csParent, 'OD', 'TD')
mtexColorMap WhiteJet
saveFigure(char("{110}[001] parent pole figures.png"));

plotPoleFigures(oriChilds, csChild, 'OD', 'TD')
mtexColorMap WhiteJet
saveFigure(char("{110}[001] child pole figures.png"));

%%
csParent = crystalSymmetry('cubic',[3.31 3.31 3.31],'mineral','Ti (beta)');
csChild = crystalSymmetry('6/mmm',[2.95 2.95 4.68], [90*degree 90*degree 120*degree],'mineral','Ti (alpha)');
beta2alpha = orientation.Burgers(csParent,csChild);

planes = [1 1 0];
directions = [0 0 1];
planes = [1 1 -2
	1 1 2];
directions = [1 1 1
	-1 -1 1];
oriParent = orientation.byMiller(planes, directions, csParent);
mParent = Miller({1 1 0}, {1 1 1}, csParent);
for i = 1:length(oriParent)
	plotPDF(oriParent(i), mParent, 'label', round2Miller(oriParent(i)))
	hold on
end

mChild = Miller({0 0 0 1},{1 1 -2 0}, csChild);
oriChildren = variants(beta2alpha, oriParent, mChild);
plotPDF(oriChildren, mChild)

%% add crystal shape
shape = crystalShape.hex(CS);
plotPDF(ori, Miller({1 1 0},CS), 'contour')
hold on
plot(ori*CS, ori*CS*shape*0.2,'add2all')

%% functions
function ori = getOrientations(planes, directions, CS)
for i = 1:size(planes, 1)
	plane = planes(i,:);
	direction = directions(i,:);
	ori(i) = orientation.byMiller(plane, direction, CS);
end
end

function ori = simluatePoleFigures(planes, directions, CS)
% ori = getOrientations(planes, directions, CS);
ori = orientation.byMiller(planes, directions, CS);
plotPoleFigures(ori, CS, 'OD', 'TD')
end

function savePoleFigures(planes, directions)
str = "";
for i = 1:size(planes, 1)
	plane = planes(i,:);
	direction = directions(i,:);
	str = str + "(" + getMillerString(plane) + ")["...
	+ getMillerString(direction) + "]";
	if i < size(planes, 1)
		str = str + "+";
	end
end
saveFigure(char(str + " pole figures.png"));
end

function str = getMillerString(hkl)
str = strrep(num2str(hkl), ' ','');
end



