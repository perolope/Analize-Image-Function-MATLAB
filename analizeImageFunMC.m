function [results1,results2,results3] = analizeImageFunMC(ImageFile,XC1,YC1,TC1,XC2,YC2,TC2,varargin)
%   Welcome to the function analizeImage, that's what it does:
%   
%   If you have a third Channel you should add it manually writting:
%   [results1,results2,results3] = analizeImageFunMC(ImageFile,XC1,YC1,TC1,XC2,YC2,TC2,'XC3',XC3,'YC3',YC3,'TC3',TC3);
%   
%   INPUTS:
%   XC1: X coordinate of the microscopy points (in nm) CHANNEL 1
%   YC1: Y coordinate of the microscopy points (in nm) CHANNEL 1
%   TC1: T time variable of each localitzation         CHANNEL 1
%   XC2: X coordinate of the microscopy points (in nm) CHANNEL 2
%   YC2: Y coordinate of the microscopy points (in nm) CHANNEL 2
%   TC2: T time variable of each localitzation         CHANNEL 2
%   XC3: X coordinate of the microscopy points (in nm) CHANNEL 3
%   YC3: Y coordinate of the microscopy points (in nm) CHANNEL 3
%   TC3: T time variable of each localitzation         CHANNEL 3
%   ImageFile: Name of the image file
%
%   OPTIONAL INPUTS:
%   PixSize: pixel size of camera (in nm), default value 160nm
%
%   OUTPUTS:
%   Results1: array with number of the area, localitzations, area and
%   density CHANNEL 1
%   Results2: array with number of the area, localitzations, area and
%   density CHANNEL 2
%   Results3: array with number of the area, localitzations, area and
%   density CHANNEL 3
%   MaskPoints: logical matrix showing 1 for points in selected areas (raws) 
%   results.txt: file with the same as the array Results

%************************************************************************
% Input check:

p = inputParser; %init parser object
validChar = @(x) ischar(x);
validNum = @(x) isreal(x); 

%Default Values
defaultPixSize=160.4; %PixSize
defaultXC3=0;
defaultYC3=0;
defaultTC3=0;

%Required and optional values

addRequired(p,'ImageFile',validChar);
addRequired(p,'XC1',validNum);
addRequired(p,'YC1',validNum);
addRequired(p,'TC1',validNum);
addRequired(p,'XC2',validNum);
addRequired(p,'YC2',validNum);
addRequired(p,'TC2',validNum);

% Channel 3 optional

addParameter(p,'XC3', defaultXC3, validNum);
addParameter(p,'YC3', defaultYC3, validNum);
addParameter(p,'TC3', defaultTC3, validNum);

addParameter(p,'PixSize', defaultPixSize, validNum);

%read input values:
parse(p, ImageFile, XC1, YC1, TC1, XC2, YC2, TC2, varargin{:});

%assign parsed values

ImageFile = p.Results.ImageFile;
XC1 = p.Results.XC1;
YC1 = p.Results.YC1;
TC1 = p.Results.TC1;
XC2 = p.Results.XC2;
YC2 = p.Results.YC2;
TC2 = p.Results.TC2;
XC3 = p.Results.XC3;
YC3 = p.Results.YC3;
TC3 = p.Results.TC3;
PixSize = p.Results.PixSize;

%PIETRO NOTE: add a check to make sure coordinates are positive values, it's
%important since we will initialize stuff to -1

%**********************************************************************

%   Plot the image of the cell and the scrollpane menu for zoom:

fm = imread(ImageFile); %read the image with full 16-bit depth
if size(fm,3)==3 % check if we have a grey image or RGB, if it's RGB we change the channel
    fmg = rgb2gray(fm);
else
    fmg = fm;
end
fm2 = imadjust(fmg); 
%%%%    1. Create a scroll panel for interactive image navigation
hFig = figure('Toolbar','none','Menubar','none');  % figure window withot default toolbar and menubar
hIm = imshow(fm2); %show image in figure window
hSP = imscrollpanel(hFig,hIm); %generate scrollpanel bars in figure
set(hSP,'Units','normalized','Position',[0 .1 1 .9]) %options of scroll panel

%%% 2. Add a Magnification Box and an Overview tool.
hMagBox = immagbox(hFig,hIm);
pos = get(hMagBox,'Position'); % box at bottom left of main figure with magnif
set(hMagBox,'Position',[0 0 pos(3) pos(4)]) %settings of magnif box
imoverview(hIm) %create another window with Overwiew of image to select zoom area

%%% Following instructions are for move real time magnification
%%% square:
%%% 3. Get the scroll panel API to programmatically control the view.
api = iptgetapi(hSP);
% 5. View the top left corner of the image.
api.setVisibleLocation(0.5,0.5)
% 6. Change the magnification to the value that just fits.
api.setMagnification(api.findFitMag())

imcontrast(hIm); %control the contrast

% here function stops to allow user modulate contrast:

msgbox('use square in Overview window to zoom-in, press any kew to continue'); 
pause; 


% Coordinates are converted to pixel units for the overlap:

% Center the image after the alingment of the Nikon Software
convx=1.65*PixSize;
convy=1.45*PixSize;
XC1=XC1+convx;
YC1=YC1+convy;
XC2=XC2+convx;
YC2=YC2+convy;
XC3=XC3+convx;
YC3=YC3+convy;


XM1=XC1/PixSize;
YM1=YC1/PixSize;
XM2=XC2/PixSize;
YM2=YC2/PixSize;
XM3=XC3/PixSize;
YM3=YC3/PixSize;

%************************************************************************

%   Ask number of areas
while true
    prompt = {'Enter number of regions to select: (max=20)'}; %text of the input
    dlgtitle = 'Analize Image';  %title of the panel
    dims = [1 55];  %dimensions of the text field
    definput = {'1'}; %default value to 1 region
    maxpoly = 20;
    answer = inputdlg(prompt,dlgtitle,dims,definput);  %collect the user-introduced answer
    nArea = str2double(answer{1}); %convert the answer from string to number
    if rem(nArea,1)==0 && nArea<maxpoly %if integer value and lower than max
        break;
    else
        msgbox('You should enter an integer positive less than the maximum, press anykey to continue');
        pause;
    end
end

%************************************************************************

%   Area selection:

hold on; %plot the areas at the same as the cell

% Initializations.

results1=ones(nArea,4); %array results 1
results2=ones(nArea,4); %array results 2
results3=ones(nArea,4); %array results 3
totarea=0; %initialize the totarea
labp=ones(nArea,2); %place of label number for each area
MaxPolygonVertices = 200; %intialize here max ginput points
MaskPoints1=false(nArea,length(XM1)); %initialize mask matrix: each raw is logical mask for points in corresponding area
MaskPoints2=false(nArea,length(XM2));
MaskPoints3=false(nArea,length(XM3));

%p2=plot(XM,YM,'.g');set(p2,'markersize',1); % plot points outside
hold on;

for i=1:nArea 
    cont=1; %counter polygon vertices
    Vx=double.empty(0,MaxPolygonVertices);  %vector Vx empty, maximum of ginput points 200
    Vy=double.empty(0,MaxPolygonVertices);  %vector Vy empty, maximum of ginput points 200
    while true
        [pxp,pyp,ExitButton] = ginput(1); %select one point
        Vx(cont)=pxp; %add the points to the polygon vector 
        Vy(cont)=pyp;
        plot(Vx(1:cont),Vy(1:cont),'.-r','HandleVisibility','off'); %plot the points at the moment 
        if (cont>2 && ExitButton==3) || (cont == MaxPolygonVertices)   %condition to break while-loop: right click on mouse or max number of vertices
              
              plot(Vx,Vy,'-y','HandleVisibility','off');plot([Vx(end),Vx(1)],[Vy(end),Vy(1)],'-y','HandleVisibility','off'); %plot polygon area, also with polyshape command
              area=polyarea(Vx,Vy)*(PixSize^2)*(10^(-6)); %calculates the area (nm^2) convert from Pix^2 to um^2
              
              %Channel 1
              
              [in] = inpolygon(XM1,YM1,Vx,Vy); %count the points inside and on the border of the polygon
              points1=numel(XM1(in)); %count the number of points inside the cell       
              density1=points1/(area); %density (points/um^2) 
              in=in';
              MaskPoints1(i,:)=in; %fill MaskPoint raw corresponding to i-th area
              
              %Channel 2
              
              [in] = inpolygon(XM2,YM2,Vx,Vy); %count the points inside and on the border of the polygon
              points2=numel(XM2(in)); %count the number of points inside the cell       
              density2=points2/(area); %density (points/um^2) 
              in=in';
              MaskPoints2(i,:)=in; %fill MaskPoint raw corresponding to i-th area
              
              %Channel 3
              
              [in] = inpolygon(XM3,YM3,Vx,Vy); %count the points inside and on the border of the polygon
              points3=numel(XM3(in)); %count the number of points inside the cell       
              density3=points3/(area); %density (points/um^2) 
              in=in';
              MaskPoints3(i,:)=in; %fill MaskPoint raw corresponding to i-th area
              
              totarea=totarea+area; %sum of all the areas for calculate the density outside              
              

              %Save the results
              results1(i,1)=i; 
              results1(i,2)=points1;
              results1(i,3)=area;
              results1(i,4)=density1;
              
              results2(i,1)=i; 
              results2(i,2)=points2;
              results2(i,3)=area;
              results2(i,4)=density2;
              
              results3(i,1)=i; 
              results3(i,2)=points3;
              results3(i,3)=area;
              results3(i,4)=density3;
              
              break;
        end
        
        % add the label points on the label points matrix
        labp(i,1)=pxp+3;
        labp(i,2)=pyp+3;
        
        cont=cont+1; %allows the condition on the first iteration 
    end
end

% Count total number of points within areas:

% Channel 1

TotInPoints1 = nnz(MaskPoints1); %total points within selected areas as non-zero elements of mask, NB: points are overcounted if overlap between aeas 
AreasPoints1 = sum(MaskPoints1,1); %sum over columns of mask: for each point, in how many areas is counted
TotInPointsCorr1 = nnz(AreasPoints1); %total points within selected areas corrected for possible overlap between areas

% Channel 2

TotInPoints2 = nnz(MaskPoints2);
AreasPoints2 = sum(MaskPoints2,1);
TotInPointsCorr2 = nnz(AreasPoints2);

% Channel 3

TotInPoints3 = nnz(MaskPoints3);
AreasPoints3 = sum(MaskPoints3,1);
TotInPointsCorr3 = nnz(AreasPoints3);

% Number and density of points outside areas:

areaout=max(XM1)*max(YM1)*PixSize^2*10^(-6)-totarea;

% For each channel

pointsout1=length(XM1)-TotInPointsCorr1;
DensityOut1=pointsout1/(areaout);

pointsout2=length(XM2)-TotInPointsCorr2;
DensityOut2=pointsout2/(areaout);

pointsout3=length(XM3)-TotInPointsCorr3;
DensityOut3=pointsout3/(areaout);

% display some results:

PercPointsIn1 = 100*(TotInPointsCorr1/length(XM1));
PercPointsOut1 = 100*(pointsout1/length(XM1));
PercOverlapPoints1 = 100*((TotInPoints1-TotInPointsCorr1)/length(XM1));

PercPointsIn2 = 100*(TotInPointsCorr2/length(XM2));
PercPointsOut2 = 100*(pointsout2/length(XM2));
PercOverlapPoints2 = 100*((TotInPoints2-TotInPointsCorr2)/length(XM2));

PercPointsIn3 = 100*(TotInPointsCorr3/length(XM3));
PercPointsOut3 = 100*(pointsout3/length(XM3));
PercOverlapPoints3 = 100*((TotInPoints3-TotInPointsCorr3)/length(XM3));

disp(strcat( ['Percentage of localizations inside selected areas of Channel 1: ',num2str(PercPointsIn1)]));
disp(strcat( ['Percentage of localizations outside selected areas of Channel 1: ',num2str(PercPointsOut1)]));
disp(strcat( ['Percentage of localizations in multiple areas of Channel 1: ',num2str(PercOverlapPoints1)]));

disp(strcat( ['Percentage of localizations inside selected areas of Channel 2: ',num2str(PercPointsIn2)]));
disp(strcat( ['Percentage of localizations outside selected areas of Channel 2: ',num2str(PercPointsOut2)]));
disp(strcat( ['Percentage of localizations in multiple areas of Channel 2: ',num2str(PercOverlapPoints2)]));

disp(strcat( ['Percentage of localizations inside selected areas of Channel 3: ',num2str(PercPointsIn3)]));
disp(strcat( ['Percentage of localizations outside selected areas of Channel 3: ',num2str(PercPointsOut3)]));
disp(strcat( ['Percentage of localizations in multiple areas of Channel 3: ',num2str(PercOverlapPoints3)]));


% Plot of results
p1=plot(XM1,YM1,'.g','DisplayName','Points Outside');set(p1,'markersize',1); % plot points outside
p2=plot(XM2,YM2,'.g','HandleVisibility','off');set(p2,'markersize',1);
p3=plot(XM3,YM3,'.g','HandleVisibility','off');set(p3,'markersize',1);
for i=1:nArea
    text(labp(i,1),labp(i,2),num2str(i),'Color','Black','FontSize',7); %add label
    if i==1
        plotmask1=MaskPoints1(i,:);
        p1=plot(XM1(plotmask1),YM1(plotmask1),'.r','DisplayName','Channel 1');set(p1,'markersize',1); % plot points inside with legend
        plotmask2=MaskPoints2(i,:);
        p2=plot(XM2(plotmask2),YM2(plotmask2),'.b','DisplayName','Channel 2');set(p2,'markersize',1);
        plotmask3=MaskPoints3(i,:);
        p3=plot(XM3(plotmask3),YM3(plotmask3),'.y','DisplayName','Channel 3');set(p3,'markersize',1);
    end
    plotmask1=MaskPoints1(i,:);
    p1=plot(XM1(plotmask1),YM1(plotmask1),'.r','HandleVisibility','off');set(p1,'markersize',1); % plot points inside without legend
    plotmask2=MaskPoints2(i,:);
    p2=plot(XM2(plotmask2),YM2(plotmask2),'.b','HandleVisibility','off');set(p2,'markersize',1);
    plotmask3=MaskPoints3(i,:);
    p3=plot(XM3(plotmask3),YM3(plotmask3),'.y','HandleVisibility','off');set(p3,'markersize',1);
end
legend('show','Location','northeastoutside')

% Retrieve selected coordinates in original scale (nm):

% Channel 1

SelMask1 = (AreasPoints1 ~= 0); %create a mask: 1 if point appear in one of the selected areas, 0 otherwise
SelX1 = XC1(SelMask1); %selected original coordinates in nm
SelY1 = YC1(SelMask1);
SelT1 = TC1(SelMask1);
SelCoords1 = [SelX1 SelY1];

% Channel 2
SelMask2 = (AreasPoints2~= 0);
SelX2 = XC2(SelMask2);
SelY2 = YC2(SelMask2);
SelT2 = TC2(SelMask2);
SelCoords2 = [SelX2 SelY2];

% Channel 3
SelMask3 = (AreasPoints3~= 0);
SelX3 = XC3(SelMask3);
SelY3 = YC3(SelMask3);
SelT3 = TC3(SelMask3);
SelCoords3 = [SelX3 SelY3];


%************************************************************************


%************************************************************************

% % Save the results at the file
% 
file = fopen('results.txt','w');
fprintf(file,'%12s %12s %12s %12s\n','Polygon num','Points','Area (nm^2)','Density (loc/um^2)');
fprintf(file,'%12s\n', 'Results Channel 1');
fprintf(file,'%12d %12d %12.2f %12.4e\n',(results1)');
fprintf(file,'%12s\n', 'Results Channel 2');
fprintf(file,'%12d %12d %12.2f %12.4e\n',(results2)');
fprintf(file,'%12s\n', 'Results Channel 3');
fprintf(file,'%12d %12d %12.2f %12.4e\n',(results3)');
fprintf(file,'%12s %12s %12s\n','Density out CH1:','Density out CH2:','Density out CH3:');
fprintf(file,'%16.4e %16.4e %16.4e\n',[DensityOut1,DensityOut2,DensityOut3]');
fclose(file);
hold off;

%************************************************************************

%   Density map of selected localizations

%plot selected localizations in original scale
figure(4); hold on;
p1=plot(SelX1,SelY1,'.b');set(p1,'markersize',1,'DisplayName','Channel 1'); % plot points outside
p2=plot(SelX2,SelY2,'.r');set(p2,'markersize',1,'DisplayName','Channel 2');
p3=plot(SelX3,SelY3,'.y');set(p3,'markersize',1,'DisplayName','Channel 3');
set(gca,'ydir','reverse'); %this will flip plot to make it look like previous figure
legend('show','Location','northeastoutside'); axis equal;
hold off;

%define edges for binning
BinSize=1000; % bin-size in nm
minX=min(XC1); minY=min(YC1);
maxX=max(XC1); maxY=max(YC1);
smallestCoord = min([minX,minY]);
biggestCoord = min([maxX, maxY]);
roundsmallestCoord=(floor(smallestCoord/BinSize))*BinSize; %smallest XY-value rounded by defect for BinSize
roundbiggestCoord=(ceil(biggestCoord/BinSize))*BinSize; %biggest XY-value rounded by excess for BinSize
edg=(roundsmallestCoord:BinSize:roundbiggestCoord); %edges values for binning
edges=cell(1,2); %transfer same edges in cell array, for X and Y, required by hist3
edges{1}=edg;
edges{2}=edg;

figure(5);
% binning and plot the density map:
[Occurr, BinCent]=hist3(SelCoords1,'Edges',edges,'CdataMode','auto'); %returns counts for each binning square, and coords of bin centers
subplot(1,3,1);
hist3(SelCoords1,'Edges',edges,'CdataMode','auto'); %plot the binned data
xlabel('X')
ylabel('Y')
title('Channel 1')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
set(gca,'ydir','reverse')
axis tight;
pbaspect([1 1 1]) %to have square graph
colorbar
view(2)

%Channel 2

minX=min(XC2); minY=min(YC2);
maxX=max(XC2); maxY=max(YC2);
smallestCoord = min([minX,minY]);
biggestCoord = min([maxX, maxY]);
roundsmallestCoord=(floor(smallestCoord/BinSize))*BinSize; %smallest XY-value rounded by defect for BinSize
roundbiggestCoord=(ceil(biggestCoord/BinSize))*BinSize; %biggest XY-value rounded by excess for BinSize
edg=(roundsmallestCoord:BinSize:roundbiggestCoord); %edges values for binning
edges=cell(1,2); %transfer same edges in cell array, for X and Y, required by hist3
edges{1}=edg;
edges{2}=edg;

% binning and plot the density map:
[Occurr, BinCent]=hist3(SelCoords2,'Edges',edges,'CdataMode','auto'); %returns counts for each binning square, and coords of bin centers
subplot(1,3,2);
hist3(SelCoords2,'Edges',edges,'CdataMode','auto'); %plot the binned data
xlabel('X')
ylabel('Y')
title('Channel 2')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
set(gca,'ydir','reverse')
axis tight;
pbaspect([1 1 1]) %to have square graph
colorbar
view(2)

%Channel 3

minX=min(XC3); minY=min(YC3);
maxX=max(XC3); maxY=max(YC3);
smallestCoord = min([minX,minY]);
biggestCoord = min([maxX, maxY]);
roundsmallestCoord=(floor(smallestCoord/BinSize))*BinSize; %smallest XY-value rounded by defect for BinSize
roundbiggestCoord=(ceil(biggestCoord/BinSize))*BinSize; %biggest XY-value rounded by excess for BinSize
edg=(roundsmallestCoord:BinSize:roundbiggestCoord); %edges values for binning
edges=cell(1,2); %transfer same edges in cell array, for X and Y, required by hist3
edges{1}=edg;
edges{2}=edg;

% binning and plot the density map:
[Occurr, BinCent]=hist3(SelCoords3,'Edges',edges,'CdataMode','auto'); %returns counts for each binning square, and coords of bin centers
subplot(1,3,3);
hist3(SelCoords3,'Edges',edges,'CdataMode','auto'); %plot the binned data
xlabel('X')
ylabel('Y')
title('Channel 3')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
set(gca,'ydir','reverse')
axis tight;
pbaspect([1 1 1]) %to have square graph
colorbar
view(2)
%save density map matrix as TIFF image: 
% Occ=uint16(Occurr); %convert the bin matrix Occur into uint16 format
% imwrite(Occ,'DensMap.tif','TIFF');  %create corresponding 16-bit TIFF


%************************************************************************
end

