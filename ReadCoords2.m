function [XCoords647,YCoords647,TCoords647,XCoordsFid,YCoordsFid,TCoordsFid]=ReadCoords(FileName, InputType, varargin)

        %%%%%   READCOORDS: reads raw data from SMLM
        %%%%% read raw data (TXT or CSV) obtained from SMLM analysis in
        %%%%% NIS-elements (Nikon N-STORM) or ONI software, and extract
        %%%%% coordinates of interest for further processing. It generates
        %%%%% txt file of 3 columns containing XYT coordinates to use for
        %%%%% further processing.
        %%%%%   
        
        %------------------------------------------------------------------
        % INPUTS:
        % FileName: name of the file with extension, e.g. 'MyFile.txt'
        % InputType: denotes the type of file, type 'N-STORM' for Nikon
        % software TXT files, or 'ONI' for CSV files from ONI.
        % 
        % N.B.
        % only two-channels can be read: main-channel (named 647) and
        % fiducial markers channels (named Fid). File format should be
        % checked.
        %
        % N-STORM files (TXT) are supposed to be 26-column, columns 4-5-13 are
        % read as X-Y-T.        %
        % ONI files (CSV) are supposed to be 11-columns, columns 3-4-2 are
        % read as X-Y-T. 
        %
        %------------------------------------------------------------------
        % OUTPUTS:
        % X-Y-TCoords647: X, Y, T(frames) coordinates of localization in
        % the main channel, named 647, expressed in nanometers
        %
        % X-Y-TCoordsFid: X, Y, T(frames) coordinates of localization in
        % the second channel (typically fiducial markers), expressed in
        % nanometers
        % 
        % 
        %------------------------------------------------------------------
        
p = inputParser; %init parser object
validChar = @(x) ischar(x);
validNum = @(x) isreal(x);
% here add something to read the coordinate file


%define defaults values for optional param:
defaultSTORMname = '405/647'; % main channel name for N-STORM data
defaultONIname = 1 ; % main channel name for ONI
defaultSTORMref = 'Bead Drift Correction'; % ref channel name for N-STORM data
defaultONIref = 0 ; % ref channel name for ONI

%define required and optional input parameters:
addRequired(p,'FileName',validChar);
addRequired(p,'InputType',validChar); 

addParameter(p,'STORMname', defaultSTORMname, validChar);
addParameter(p,'STORMref', defaultSTORMref, validChar);
addParameter(p,'ONIname', defaultONIname, validNum);
addParameter(p,'ONIref', defaultONIref, validNum);

%read input values:
parse(p, FileName, InputType, varargin{:});

%assign the parsed values:
FileName = p.Results.FileName; % file name
InputType = p.Results.InputType;
STORMname = p.Results.STORMname;
STORMref = p.Results.STORMref;
ONIname = p.Results.ONIname;
ONIref = p.Results.ONIref;

switch InputType
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%------ TXT file from N-STORM --------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 'N-STORM'

        %%%----- Read text file (check format)---
        disp('Importing N-STORM data...');
        fileID = fopen(FileName,'r');
        DataIn = textscan(fileID,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','Delimiter','\t','HeaderLines',1); % Saving data as a cell, check the numb of columns 
        fclose(fileID);

        %%%--- Retrieve useful information (channel, corrected XY location and frame number or the localization)
        XCoords = cell2mat(DataIn(4)); %% colum 2 X non drift corrected, 4 X drift corrected
        YCoords = cell2mat(DataIn(5)); %% colum 3 Y non drift corrected, 5 Y drift corrected
        TCoords = cell2mat(DataIn(13)); %% column 13 Frame
        Channel = DataIn(1);
        clear DataIn;

        %%%--- Split channel based on first column and generate output
        %%%--- main channel is 647 (first column of TXT file)
        Channel647=strcmp(Channel{1},STORMname);%Check in the .txt if 647 or 405/647
        XCoords647=XCoords(Channel647);
        YCoords647=YCoords(Channel647);
        TCoords647=TCoords(Channel647);
        ChannelFid=strcmp(Channel{1},STORMref);%Check in the .txt if 647 or other name
        XCoordsFid=XCoords(ChannelFid);
        YCoordsFid=YCoords(ChannelFid);
        TCoordsFid=TCoords(ChannelFid);
        
        %%%--- Export a txt file with XYT coords
        n=length(XCoords647);
        nref=length(XCoordsFid);
        A=zeros(n,3);
        Ref=zeros(nref,3);
        A(:,1)=XCoords647;
        A(:,2)=YCoords647;
        A(:,3)=TCoords647;
        Ref(:,1)=XCoordsFid;
        Ref(:,2)=YCoordsFid;
        Ref(:,3)=TCoordsFid;
        save XYTcoordinates.txt A -ascii
        save XYTref.txt Ref -ascii
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%------ CSV file from ONI ------------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 'ONI'
        
        %%%----- Read text file (check format)---
        disp('Importing ONI data...');
        fileID = fopen(FileName,'r');
        DataIn = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f %f','Delimiter',',','HeaderLines',1); % Saving data as a cell, check the numb of columns 
        fclose(fileID);
        
        %%%--- Retrieve useful information (channel, corrected XY location and frame number or the localization)
        XCoords = cell2mat(DataIn(3)); %% colum 3 X coords
        YCoords = cell2mat(DataIn(4)); %% colum 4 Y coords
        TCoords = cell2mat(DataIn(2)); %% column 2 Frame coords
        Channel = DataIn(1);
        clear DataIn;
        
        
        %%%--- Split channel based on first column and generate output
        Channel647= (Channel{1}==ONIname); % 647-channel  with 1
        XCoords647=XCoords(Channel647);
        YCoords647=YCoords(Channel647);
        TCoords647=TCoords(Channel647);
        ChannelFid= (Channel{1}==ONIref);%C Fiducial-channels denoted with 0
        XCoordsFid=XCoords(ChannelFid);
        YCoordsFid=YCoords(ChannelFid);
        TCoordsFid=TCoords(ChannelFid);
        
        %%%--- Export a txt file with XYT coords
        n=length(XCoords647);
        nref=length(XCoordsFid);
        A=zeros(n,3);
        Ref=zeros(nref,3);
        A(:,1)=XCoords647;
        A(:,2)=YCoords647;
        A(:,3)=TCoords647;
        Ref(:,1)=XCoordsFid;
        Ref(:,2)=YCoordsFid;
        Ref(:,3)=TCoordsFid;
        save XYTcoordinates.txt A -ascii
        save XYTref.txt Ref -ascii
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%------ CSV file from Thunderstorm ---------------------------------
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     
    case 'THUNDER'
    
        %%%----- Read text file (check format)---
        disp('Importing Thunder data...');
        fileID = fopen(FileName,'r');
        DataIn = textscan(fileID,'%f %f %f %f %f %f %f','Delimiter',',','HeaderLines',1); % Saving data as a cell, check the numb of columns 
        fclose(fileID);
        
        %%%--- Retrieve useful information (channel, corrected XY location and frame number or the localization)
        XCoords = cell2mat(DataIn(3));
        YCoords = cell2mat(DataIn(4));
        TCoords = cell2mat(DataIn(2));
        XCoordsFid = 0; 
        YCoordsFid = 0; 
        TCoordsFid = 0;
%         Channel = DataIn(1);
        %clear DataIn;
        
        
        %%%--- Split channel based on first column and generate output
        %Channel647= (Channel{1}==1); % 647-channel denoted with 1
        XCoords647=XCoords;
        YCoords647=YCoords;
        TCoords647=TCoords;
%         ChannelFid= (Channel{1}==0);%C Fiducial-channels denoted with 0
%         XCoordsFid=XCoords(ChannelFid);
%         YCoordsFid=YCoords(ChannelFid);
%         TCoordsFid=TCoords(ChannelFid);
%         
    otherwise
        disp('invalid InputType!');
end
end