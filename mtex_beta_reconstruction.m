%% Usage
% this code will produce images in the current folder,
% make sure to run this code in working folder (not the folder it
% lives in)

global filename
%% initialise mtex
folderpath = fileparts(mfilename('fullpath'));
if ~ismember(folderpath, strsplit(path, pathsep))
	startup
end
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outofPlane');
pfAnnotations = @(varargin) [];
setMTEXpref('pfAnnotations',pfAnnotations);

%% Import EBSD Data
file = dir('*.cpr');
assert(length(file) == 1);
filename = file.name;
filepath = file.folder;
disp("processing " + filename + "...");
if ~isa(filename, 'char')
    return
end
[~,filename,file_ext] = fileparts(filename);

CS = {...
    'notIndexed',...
    crystalSymmetry('m-3m', [3.2 3.2 3.2], 'mineral', 'bcc', 'color', 'r'),... [0.53 0.81 0.98]
    crystalSymmetry('6/mmm', [3 3 4.7], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'hcp', 'color', 'b'),... [0.56 0.74 0.56]
    'notIndexed',...
    % 	crystalSymmetry('6/mmm', [7.5 7.5 5.2], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Si3Ti5', 'color', 'g'),...
    };
ebsd = EBSD.load(fullfile(filepath, [filename,file_ext]),CS,'interface','crc',...
    'convertSpatial2EulerReferenceFrame');

%% Main
ebsd_cleaned = cleaning(ebsd);

rot = rotation.byAxisAngle(vector3d.X, 90*degree); % rotate first to align X//OD, Z//FD
ebsd_cleaned = rotate(ebsd_cleaned, rot, 'keepXY');

plotEBSD(ebsd_cleaned)
angular_tol = 1.5*degree;
reconstruction = betaReconstruction(ebsd_cleaned, angular_tol);
parentGrains = reconstruction.parentGrains;
parentEBSD = reconstruction.parentEBSD;
plotParent(parentGrains, parentEBSD);

rot = rotation.byAxisAngle(vector3d.Z, -19*degree);
ebsd_rotated = rotate(ebsd_cleaned, rot, 'keepXY');
plotEBSD_rotated(ebsd_rotated);

parentEBSD_rotated = rotate(parentEBSD, rot, 'keepXY');
parentGrains_rotated = rotate(parentGrains, rot, 'keepXY');
plotParent_rotated(parentGrains_rotated, parentEBSD_rotated);

close all


function ebsd_cleaned = cleaning(ebsd)
global filename
fullfilename = [filename, ' cleaned.mat'];
if exist(fullfilename, 'file')
    ebsd_cleaned = load(fullfilename).ebsd_cleaned;
    return
end
angular_tol = 10*degree;
ebsd_cleaned = ebsd('indexed'); % consider only indexed data
[grains, ebsd_cleaned.grainId, ebsd_cleaned.mis2mean] = calcGrains(ebsd_cleaned,'angle',angular_tol);
ebsd_cleaned(grains(grains.grainSize<5)) = [];
[grains, ebsd_cleaned.grainId] = calcGrains(ebsd_cleaned,'angle',angular_tol);
grains = smooth(grains,5);

F = halfQuadraticFilter;
F.alpha = 0.25;
ebsd_cleaned = smooth(ebsd_cleaned,F,'fill',grains);

save(fullfilename, 'ebsd_cleaned')
end

function plotEBSD(ebsd)
global filename
%% plot IPFs
fullfilename = [filename, ' IPFZ micronbar.png'];
if ~exist(fullfilename, 'file')
    figure
    plot(ebsd('bcc'),ebsd('bcc').orientations)
    hold on
    plot(ebsd('hcp'),ebsd('hcp').orientations)
    saveFigure(fullfilename);
end

fullfilename = [filename, ' IPFZ.png'];
if ~exist(fullfilename, 'file')
    figure
    plot(ebsd('bcc'),ebsd('bcc').orientations,'micronbar','off')
    hold on
    plot(ebsd('hcp'),ebsd('hcp').orientations,'micronbar','off')
    saveFigure(fullfilename);
end

fullfilename = [filename, ' IPFY.png'];
if ~exist(fullfilename, 'file')
    figure
    if ~isempty(ebsd('bcc'))
        ipfKey = ipfColorKey(ebsd('bcc'));
        ipfKey.inversePoleFigureDirection = vector3d.Y;
        colors = ipfKey.orientation2color(ebsd('bcc').orientations);
        plot(ebsd('bcc'),colors,'micronbar','off')
    end
    hold on
    ipfKey = ipfColorKey(ebsd('hcp'));
    ipfKey.inversePoleFigureDirection = vector3d.Y;
    colors = ipfKey.orientation2color(ebsd('hcp').orientations);
    plot(ebsd('hcp'),colors,'micronbar','off')
    saveFigure(fullfilename);
end

%% plot phase map
fullfilename = [filename, ' phase.png'];
if ~exist(fullfilename, 'file')
    figure
    plot(ebsd('bcc'),'micronbar','off')
    hold on
    plot(ebsd('hcp'),'micronbar','off')
    legend off
    saveFigure(fullfilename);
end

%% plot pole figures
fullfilename = [filename, ' alpha pole figures.png'];
if ~exist(fullfilename, 'file') && ~isempty(ebsd('hcp'))
    plotPoleFigures(ebsd('hcp').orientations, ebsd('hcp').CS, 'X', 'Y')
	saveFigure(fullfilename);
end
fullfilename = [filename, ' beta pole figures.png'];
if ~exist(fullfilename, 'file') && ~isempty(ebsd('bcc'))
	plotPoleFigures(ebsd('bcc').orientations, ebsd('bcc').CS, 'X', 'Y')
	saveFigure(fullfilename);
end

%% plot KAM
fullfilename = [filename, ' KAM.png'];
if ~exist(fullfilename, 'file')
    figure
    [~,ebsd.grainId] = calcGrains(ebsd('indexed'));
    kam = ebsd.KAM('threshold',2.5*degree) ./ degree;
    plot(ebsd,kam,'micronbar','off')
    mtexColorbar
    mtexColorMap BlueJet
    saveFigure(fullfilename);
end
end

function plotEBSD_rotated(ebsd)
global filename
%% plot IPFs
% if ~exist([filename, ' IPFZ rotated.png'], 'file')
%     figure
%     if ~isempty(ebsd('bcc'))
%         plot(ebsd('bcc'),ebsd('bcc').orientations,'micronbar','off')
%     end
%     hold on
%     plot(ebsd('hcp'),ebsd('hcp').orientations,'micronbar','off')
% 
%     saveFigure([filename, ' IPFZ rotated.png']);
% end

if ~exist([filename, ' IPFY rotated.png'], 'file')
    figure
    if ~isempty(ebsd('bcc'))
        ipfKey = ipfColorKey(ebsd('bcc'));
        ipfKey.inversePoleFigureDirection = vector3d.Y;
        colors = ipfKey.orientation2color(ebsd('bcc').orientations);
        plot(ebsd('bcc'),colors,'micronbar','off')
    end
    hold on
    ipfKey = ipfColorKey(ebsd('hcp'));
    ipfKey.inversePoleFigureDirection = vector3d.Y;
    colors = ipfKey.orientation2color(ebsd('hcp').orientations);
    plot(ebsd('hcp'),colors,'micronbar','off')

    saveFigure([filename, ' IPFY rotated.png']);
end

%% plot pole figures
fullfilename = [filename, ' alpha pole figures rotated.png'];
if ~exist(fullfilename, 'file') && ~isempty(ebsd('hcp'))
	plotPoleFigures(ebsd('hcp').orientations, ebsd('hcp').CS, 'OD', 'TD')
	saveFigure(fullfilename);
end
fullfilename = [filename, ' beta pole figures rotated.png'];
if ~exist(fullfilename, 'file') && ~isempty(ebsd('bcc'))
    plotPoleFigures(ebsd('bcc').orientations, ebsd('bcc').CS, 'OD', 'TD')
	saveFigure(fullfilename);
end
end

function reconstruction = betaReconstruction(ebsd, angular_tol)
global filename
if exist([filename, ' reconstruction.mat'], 'file')
    reconstruction = load([filename, ' reconstruction.mat']).reconstruction;
    return
end
% figure
beta2alpha = orientation.Burgers(ebsd('bcc').CS,ebsd('hcp').CS);
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'threshold',angular_tol);
reconstruction = parentGrainReconstructor(ebsd, grains);
reconstruction.p2c = beta2alpha;

reconstruction.calcGraph('threshold',angular_tol)
reconstruction.calcVariants('threshold',angular_tol)
reconstruction.clusterGraph('numIter',3)
% reconstruction.plotGraph
reconstruction.calcParentFromGraph

% ipfKey = ipfColorKey(ebsd('bcc'));
% ipfKey.inversePoleFigureDirection = vector3d.Z;
% color = ipfKey.orientation2color(reconstruction.parentGrains.meanOrientation);
% plot(reconstruction.parentGrains, color)

reconstruction.mergeSimilar('threshold',angular_tol)
% color = ipfKey.orientation2color(reconstruction.parentGrains.meanOrientation);
% plot(reconstruction.parentGrains,color)

reconstruction.mergeInclusions('maxSize',500)
% color = ipfKey.orientation2color(reconstruction.parentGrains.meanOrientation);
% plot(reconstruction.parentGrains,color)
% saveFigure(fullfile(output_path, [filename, ' prior beta IPFZ.png']));
save([filename, ' reconstruction.mat'], 'reconstruction')
end

function plotParent(grains, ebsd)
global filename
%% plot prior beta IPFZ
if ~exist([filename, ' prior beta IPFZ.png'], 'file')
    figure
    ipfKey = ipfColorKey(ebsd('bcc'));
    ipfKey.inversePoleFigureDirection = vector3d.Z;
    color = ipfKey.orientation2color(grains.meanOrientation);
    plot(grains,color,'micronbar','off')
    saveFigure([filename, ' prior beta IPFZ.png']);
end
%% plot prior beta IPFZ with crystal
if ~exist([filename, ' prior beta IPFZ with crystal.png'], 'file')
    figure
    plot(grains,color,'micronbar','off')
    hold on
    isBig = grains.grainSize > 100;
    shape = crystalShape.cube(ebsd('bcc').CS);
    plot(grains(isBig),0.5*shape, 'facecolor', [0.53 0.81 0.98])
    saveFigure([filename, ' prior beta IPFZ with crystal.png']);
end

%% plot prior beta IPFY
if ~exist([filename, ' prior beta IPFY.png'], 'file')
    figure
    ipfKey = ipfColorKey(ebsd('bcc'));
    ipfKey.inversePoleFigureDirection = vector3d.Y;
    color = ipfKey.orientation2color(grains.meanOrientation);
    plot(grains,color,'micronbar','off')
    saveFigure([filename, ' prior beta IPFY.png']);
end

%% plot prior beta IPFY with crystal
if ~exist([filename, ' prior beta IPFY with crystal.png'], 'file')
    figure
    plot(grains,color,'micronbar','off')
    hold on
    isBig = grains.grainSize > 100;
    shape = crystalShape.cube(ebsd('bcc').CS);
    rot = rotation.byAxisAngle(vector3d.X, 90*degree);
    plot(grains(isBig).centroid + rot * grains(isBig).meanOrientation * shape * 0.5 * sqrt(grains(isBig).area), 'facecolor', [0.53 0.81 0.98])
    saveFigure([filename, ' prior beta IPFY with crystal.png']);
end

%% plot prior beta pole figures
if ~exist([filename, ' prior beta pole figures.png'], 'file')
    plotPoleFigures(ebsd('bcc').orientations, ebsd('bcc').CS, 'X', 'Y')
    saveFigure([filename, ' prior beta pole figures.png']);
end
end

function plotParent_rotated(grains, ebsd)
global filename
%% plot prior beta IPFZ
% if ~exist([filename, ' prior beta IPFZ rotated.png'], 'file')
%     figure
%     ipfKey = ipfColorKey(ebsd('bcc'));
%     ipfKey.inversePoleFigureDirection = vector3d.Z;
%     color = ipfKey.orientation2color(grains.meanOrientation);
%     plot(grains,color,'micronbar','off')
%     saveFigure([filename, ' prior beta IPFZ rotated.png']);
% end

%% plot prior beta IPFZ with crystal
if ~exist([filename, ' prior beta IPFZ with crystal rotated.png'], 'file')
    figure
	ipfKey = ipfColorKey(ebsd('bcc'));
    ipfKey.inversePoleFigureDirection = vector3d.Z;
    color = ipfKey.orientation2color(grains.meanOrientation);
    plot(grains,color,'micronbar','off')
    hold on
    isBig = grains.grainSize > 100;
    shape = crystalShape.cube(ebsd('bcc').CS);
    plot(grains(isBig),0.5*shape, 'facecolor', [0.53 0.81 0.98])
    saveFigure([filename, ' prior beta IPFZ with crystal rotated.png']);
end

%% plot prior beta IPFY
if ~exist([filename, ' prior beta IPFY rotated.png'], 'file')
    figure
    ipfKey = ipfColorKey(ebsd('bcc'));
    ipfKey.inversePoleFigureDirection = vector3d.Y;
    color = ipfKey.orientation2color(grains.meanOrientation);
    plot(grains,color,'micronbar','off')
    saveFigure([filename, ' prior beta IPFY rotated.png']);
end

%% plot prior beta IPFY with crystal
if ~exist([filename, ' prior beta IPFY with crystal rotated.png'], 'file')
    figure
    ipfKey = ipfColorKey(ebsd('bcc'));
    ipfKey.inversePoleFigureDirection = vector3d.Y;
    color = ipfKey.orientation2color(grains.meanOrientation);
    plot(grains,color,'micronbar','off')
    hold on
    isBig = grains.grainSize > 100;
    shape = crystalShape.cube(ebsd('bcc').CS);
    rot = rotation.byAxisAngle(vector3d.X, 90*degree);
    plot(grains(isBig).centroid + rot * grains(isBig).meanOrientation * shape * 0.5 * sqrt(grains(isBig).area), 'facecolor', [0.53 0.81 0.98])
    saveFigure([filename, ' prior beta IPFY with crystal rotated.png']);
end

%% plot prior beta pole figures
if ~exist([filename, ' prior beta pole figures rotated.png'], 'file')
    plotPoleFigures(ebsd('bcc').orientations, ebsd('bcc').CS, 'OD', 'TD')
    saveFigure([filename, ' prior beta pole figures rotated.png']);
end
end
