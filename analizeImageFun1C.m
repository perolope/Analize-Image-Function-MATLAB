function [results, MaskPoints, cellarray] = analizeImageFun1C(X,Y,T,ImageFile,varargin)
%   Welcome to the function analizeImage, that's what it does:
%   
%   
%   INPUTS:
%   X: X coordinate of the microscopy points (in nm)
%   Y: Y coordinate of the microscopy points (in nm)
%   T: T time variable of each localitzation
%   ImageFile: Name of the image file
%
%   OPTIONAL INPUTS:
%   PixSize: pixel size of camera (in nm), default value 160nm
%
%   OUTPUTS:
%   Results: array with number of the area, localitzations, area and
%   density
%   MaskPoints: logical matrix showing 1 for points in selected areas (raws) 
%   results.txt: file with the same as the array Results
%   cellarray: cellarray with X,Y,T of each ROI for then do time trace

%************************************************************************
% Input check:

p = inputParser; %init parser object
validChar = @(x) ischar(x);
validNum = @(x) isreal(x); 

%Default Values
defaultPixSize=160.4; %PixSize

%Required and optional values
addRequired(p,'X',validNum);
addRequired(p,'Y',validNum);
addRequired(p,'T',validNum);
addRequired(p,'ImageFile',validChar);

addParameter(p,'PixSize', defaultPixSize, validNum);

%read input values:
parse(p, X, Y, T, ImageFile, varargin{:});

%assign parsed values
X = p.Results.X;
Y = p.Results.Y;
T = p.Results.T;
ImageFile = p.Results.ImageFile;
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

% test of correct the drift of the image
convx=1.65*PixSize;
convy=1.45*PixSize;
X=X+convx;
Y=Y+convy;

XM=X/PixSize;
YM=Y/PixSize;

%cellarray={[XM,YM,T]}; %cell array needed for time traces
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

results=ones(nArea,4); %array results
totarea=0; %initialize the totarea
labp=ones(nArea,2); %place of label number for each area
MaxPolygonVertices = 200; %intialize here max ginput points
MaskPoints=false(nArea,length(XM)); %initialize mask matrix: each raw is logical mask for points in corresponding area
cellarray=cell(1,nArea); %create the cell array to fill it later

hold on;

for i=1:nArea 
    cont=1; %counter polygon vertices
    Vx=double.empty(0,MaxPolygonVertices);  %vector px filled with -1, maximum of ginput points 200
    Vy=double.empty(0,MaxPolygonVertices);  %vector py filled with -1, maximum of ginput points 200
    while true
        [pxp,pyp,ExitButton] = ginput(1); %select one point
        Vx(cont)=pxp; %add the points to the polygon vector 
        Vy(cont)=pyp;
        plot(Vx(1:cont),Vy(1:cont),'.-r','HandleVisibility','off'); %plot the points at the moment 
        if (cont>2 && ExitButton==3) || (cont == MaxPolygonVertices)   %condition to break while-loop: right click on mouse or max number of vertices
              [in] = inpolygon(XM,YM,Vx,Vy); %count the points inside and on the border of the polygon
              in=in';

              plot(Vx,Vy,'-y','HandleVisibility','off');plot([Vx(end),Vx(1)],[Vy(end),Vy(1)],'-y','HandleVisibility','off'); %plot polygon area, also with polyshape command
              
              points=numel(XM(in)); %count the number of points inside the cell
              area=polyarea(Vx,Vy)*(PixSize^2)*(10^(-6)); %calculates the area (nm^2) convert from Pix^2 to um^2
              density=points/(area); %density (points/um^2) 
              totarea=totarea+area; %sum of all the areas for calculate the density outside
              MaskPoints(i,:)=in; %fill MaskPoint raw corresponding to i-th area
                
              cellarray{i}=[XM(in),YM(in),T(in)]; %fill the cell array with X,Y,T
              
              %Save the results
              results(i,1)=i; 
              results(i,2)=points;
              results(i,3)=area;
              results(i,4)=density;
              
              break;
        end
        
        % add the label points on the label points matrix
        labp(i,1)=pxp+3;
        labp(i,2)=pyp+3;
        
        cont=cont+1; %allows the condition on the first iteration 
    end
end

% Count total number of points within areas:

TotInPoints = nnz(MaskPoints); %total points within selected areas as non-zero elements of mask, NB: points are overcounted if overlap between aeas 
AreasPoints = sum(MaskPoints,1); %sum over columns of mask: for each point, in how many areas is counted
TotInPointsCorr = nnz(AreasPoints); %total points within selected areas corrected for possible overlap between areas

% Number and density of points outside areas:

area2=max(XM)*max(YM)*PixSize^2*10^(-6)-totarea;
pointsout=length(XM)-TotInPointsCorr;
DensityOut=pointsout/(area2);

% display some results:

PercPointsIn = 100*(TotInPointsCorr/length(XM));
PercPointsOut = 100*(pointsout/length(XM));
PercOverlapPoints = 100*((TotInPoints-TotInPointsCorr)/length(XM));
disp(strcat( ['Percentage of localizations inside selected areas: ',num2str(PercPointsIn)]));
disp(strcat( ['Percentage of localizations outside selected areas: ',num2str(PercPointsOut)]));
disp(strcat( ['Percentage of localizations in multiple areas: ',num2str(PercOverlapPoints)]));


% Plot of results
p2=plot(XM,YM,'.g','DisplayName','Points Outside');set(p2,'markersize',1); % plot points outside
for i=1:nArea
    text(labp(i,1),labp(i,2),num2str(i),'Color','Black','FontSize',7); %add label
    plotmask=MaskPoints(i,:);
    if i==1
        p1=plot(XM(plotmask),YM(plotmask),'.y','DisplayName','Points Inside');set(p1,'markersize',1); % plot points inside
    end
    p1=plot(XM(plotmask),YM(plotmask),'.y','HandleVisibility','off');set(p1,'markersize',1); % plot points inside
end
%legend('show','Location','northeastoutside');
% Retrieve selected coordinates in original scale (nm):

SelMask = (AreasPoints ~= 0); %create a mask: 1 if point appear in one of the selected areas, 0 otherwise
SelX = X(SelMask); %selected original coordinates in nm
SelY = Y(SelMask);
SelT = T(SelMask);
SelCoords = [SelX SelY];
%************************************************************************


%************************************************************************

% % Save the results at the file
% 
ImageName=erase(ImageFile,'.tif'); %erase the part .tif
Name=strcat(ImageName,'_RESULTS','.txt'); %create the name of the results file
file = fopen(Name,'w');
fprintf(file,'%12s %12s %12s %12s\n','Polygon num','Points','Area (nm^2)','Density (loc/um^2)');
fprintf(file,'%12d %12d %12.2f %12.4e\n',(results)');
fprintf(file,'%12s\n','Density outside:');
fprintf(file,'%12.4e\n',DensityOut);
fclose(file);
hold off;

%************************************************************************

%   Density map of selected localizations

%plot selected localizations in original scale
figure(4);
p2=plot(SelX,SelY,'.b');set(p2,'markersize',1); axis equal; % plot points outside
set(gca,'ydir','reverse'); %this will flip plot to make it look like previous figure

%define edges for binning
BinSize=1000; % bin-size in nm
minX=min(X); minY=min(Y);
maxX=max(X); maxY=max(Y);
smallestCoord = min([minX,minY]);
biggestCoord = min([maxX, maxY]);
roundsmallestCoord=(floor(smallestCoord/BinSize))*BinSize; %smallest XY-value rounded by defect for BinSize
roundbiggestCoord=(ceil(biggestCoord/BinSize))*BinSize; %biggest XY-value rounded by excess for BinSize
edg=(roundsmallestCoord:BinSize:roundbiggestCoord); %edges values for binning
edges=cell(1,2); %transfer same edges in cell array, for X and Y, required by hist3
edges{1}=edg;
edges{2}=edg;

% binning and plot the density map:
[Occurr, BinCent]=hist3(SelCoords,'Edges',edges,'CdataMode','auto'); %returns counts for each binning square, and coords of bin centers
figure(5);
hist3(SelCoords,'Edges',edges,'CdataMode','auto'); %plot the binned data
xlabel('X')
ylabel('Y')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
set(gca,'ydir','reverse')
axis tight;
pbaspect([1 1 1]) %to have square graph
colorbar
view(2)

% %save density map matrix as TIFF image: 
% Occ=uint16(Occurr); %convert the bin matrix Occur into uint16 format
% imwrite(Occ,'DensMap.tif','TIFF');  %create corresponding 16-bit TIFF


% %************************************************************************
end

