%% Line Processing Algorithm
% EL image goes through a series of image processing steps to detect the
% orientation and length of microcracks using the Hough Transform. The
% orientation and length is saved in a structure array. The percent of
% black pixels for each image is saved in an array.

tic
import mlreportgen.dom.*

% number of solar cells
N = 10; 
% cell to save images
d = cell(N,1);
% array to save percent of black pixels in an image
percent_black = zeros(N,1);
% structure to save length (mm) and orientation (deg) of crack 
line_data = struct();

% iterate through each image
for i=1:N

filename = strcat('rig_v',num2str(i),'_1s.tiff');
A = imread(filename);

% crop image
A = imcrop(A, [752.5 472.5 2416 2416]);

% convert to grayscale
gs = rgb2gray(A);

% apply median filtering to remove noise
filteredA = medfilt2(gs, [5 5]);

% normalize the image to account for varying lighting conditions
filteredA = imadjust(filteredA,stretchlim(filteredA));

% perform adaptive binarization
mask = imbinarize(filteredA, 'adaptive', 'ForegroundPolarity', 'dark', 'Sensitivity', 0.5);

% uncomment to detect vertical crack in validation set
mask_closed = imclose(mask,strel('line',20,90));

% count the percent of black pixels using a small area of interest
targetSize = [2217 2217];
r = centerCropWindow2d(size(mask),targetSize);
J = imcrop(mask,r); 
% calculates percent of black pixels
percent_black(i)=nnz(~J)/nnz(J)*100;

% perform canny edge detection
% outline = edge(mask,'Canny',.05);

% uncomment to detect vertical crack
outline = edge(mask_closed,'Canny',.05);

% apply hough transform
[H,T,R] = hough(outline);

% Identify Peaks in Hough Transform 
% neighborhood suppresses lines that appear on top of one another
hPeaks =  houghpeaks(H,200,'NHoodSize',[351 11]); 

% extract lines from hough transform and peaks 
hLines = houghlines(outline,T,R,hPeaks, 'FillGap', 700);

% overlay lines
[linePos,markerPos] = getVizPosArray(hLines);

% create cell to save lengths
lengths = cell(size(linePos,1), 1);

% create double array to store angles
angles = zeros(size(linePos,1), 1);

% iterate through all the lines detected and save length and orientation
for k = 1:length(hLines)
    check = cat(2,hLines(k).point1,hLines(k).point2);
    for m = 1:size(linePos,1)
        if ismember(linePos(m, 1:4), check) == [1 1 1 1]
            % If want to find angle that goes with endpoints, they have the
            % same index in their respective array.
            lengths{m, 1} = norm(hLines(k).point1 - hLines(k).point2) * 0.05; % 0.05 converts from px to mm
            angles(m) = hLines(k).theta;
        end
    end
end

% save length and angles into structure array
line_data(i).lengths = lengths;
line_data(i).angles = angles;
    
% overlay lines detected on original cropped image
if size(markerPos,1) ~= 0 
    lineFrame = insertShape(A,'Line',linePos,'Color','red','LineWidth',15);
    outFrame = insertObjectAnnotation(lineFrame,...
            'circle',markerPos,'','Color','yellow','LineWidth',5);
else
    outFrame = A;
end
  
% display montage of 
d{i} = outFrame; 
end 
montage(d,"size",[2 5]);
toc


function [linePos,markerPos] = getVizPosArray(hLines)

% Copyright 2015-2016 The MathWorks, Inc.
% Get line and marker positions (output) of houghlines (input) formatted to
% use directly with insertShape and insertObjectAnnotation functions.

% Modified to ensure that edges of solar cell are not detected
linePos = [];
markerPos = [];
markerSize = 10;

for lidx = 1:length(hLines)
    % only looks at lines that are greater than 270 px in length
    if norm(hLines(lidx).point1 - hLines(lidx).point2) > 270 
        % vertical cracks detected 
        if (hLines(lidx).theta > -5 && hLines(lidx).theta < 5) 
            % x coordinates of region of interest (roi)
            xv = [710; 1709; 1709; 710; 710;];
            % y coordinates of roi
            yv = [1; 1; 2417; 2417; 1;];
            % x and y coordinates of enpoints detected
            xq = [hLines(lidx).point1(1); hLines(lidx).point2(1)];
            yq = [hLines(lidx).point1(2); hLines(lidx).point2(2)];
            % finds whether endpoints detected lie within or on the roi
            [in,on] = inpolygon(xq,yq,xv,yv);
            % only save the enpoints of the vertical lines that are within
            % the roi (ignores vertical edges of solar cells)
            if numel(xq(in&~on)) ~= 0 && numel(yq(in&~on)) ~= 0                
                linePos = [linePos; hLines(lidx).point1 hLines(lidx).point2];  % Get linePos in [x1 y1;x2 y2...] format.   
                markerPos = [markerPos;...
                                [hLines(lidx).point1 markerSize];...
                                [hLines(lidx).point2 markerSize]]; % Get markerPos in [x1 markerSize;y1 markerSize...] format.  
            end
        % horizontal lines detected 
        elseif hLines(lidx).theta >= -90 && hLines(lidx).theta <= -85 || hLines(lidx).theta >= 85 && hLines(lidx).theta <= 90 
            % only save the horizontal lines within the roi
            yv = [310; 2109; 2109; 310; 310;];
            xv = [1; 1; 2417; 2417; 1;];
            xq = [hLines(lidx).point1(1); hLines(lidx).point2(1)];
            yq = [hLines(lidx).point1(2); hLines(lidx).point2(2)];
            [in,on] = inpolygon(xq,yq,xv,yv);
            if numel(xq(in&~on)) ~= 0 && numel(yq(in&~on)) ~= 0
                linePos = [linePos; hLines(lidx).point1 hLines(lidx).point2];  % Get linePos in [x1 y1;x2 y2...] format.   
                markerPos = [markerPos;...
                            [hLines(lidx).point1 markerSize];...
                            [hLines(lidx).point2 markerSize]]; % Get markerPos in [x1 markerSize;y1 markerSize...] format. 
            end
        else
            linePos = [linePos; hLines(lidx).point1 hLines(lidx).point2];  % Get linePos in [x1 y1;x2 y2...] format.   
            markerPos = [markerPos;...
                            [hLines(lidx).point1 markerSize];...
                            [hLines(lidx).point2 markerSize]]; % Get markerPos in [x1 markerSize;y1 markerSize...] format.  
        end
    end
end
end