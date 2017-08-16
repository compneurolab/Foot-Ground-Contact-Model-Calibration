function CreateFootSpringMesh_v3
% This file will determine foot spring positions relative to the heel and toe
% markers, which can be used to create a foot ground contact model.
% Must specify
%   Wand trace TRC file
%   Wand marker names
%   Distance from wand tip to each marker on the wand (optional)
%   Foot marker names (Heel, Toe, Medial toe, Lateral toe)
%   Name of the subject specific model file used in later simulations.
%   Output file names
% Note: all gaps in trc file data must be filled in with 0s for this
% program to function correctly

import org.opensim.modeling.*
import java.lang.*

close all

%% Settings

%load model to determine model spring locations
%osimModel = 'Patient5_optModel4_GPOPS.osim'
osimModel = 'Patient4_optModel.osim'

model=Model(osimModel);

FootSide = 'L';

TraceTRC = ['Patient4_Wand' FootSide '_mokka.trc'];

Range = 500:1000; % Time range to determine foot experimental marker locations

% Give transform to convert coordinates to SIMM coordinates
% Enter [1 2 3] to keep same order
AxisTransform = [2 3 1];

WandNames = {'Wand.End'
    'Wand.Tip'
    'Wand.Mid1'
    'Wand.Mid2'};

% Specify index of marker in names listed above
WandEnd = 1;
WandTip = 2; %optional if  2 midpoint markers are specified
WandMid1 = 0; %optional, enter 0 to ignore
WandMid2 = 0; %optional, enter 0 to ignore

ToeMarker = [FootSide '.Toe'];
HeelMarker = [FootSide '.Heel'];
ToeLateral = [FootSide '.Toe.Lateral'];
ToeMedial = [FootSide '.Toe.Medial'];
MidfootSup = [FootSide '.Midfoot.Superior'];
MidfootLat = [FootSide '.Midfoot.Lateral'];

% Specify numer of rows of springs and columns of springs along the length
% and width of the foot
numLength = 12;
numWidth = 6;

%% Load data
% Load TRC file with marker data
[MarkerData, MarkerNames, SampleRate, Time] = LoadTRCFile(TraceTRC);

nframes = numel(Time);
MarkerData = MarkerData/1000;

MarkerDataOrig = MarkerData;

%% Transform Data
%Transform marker data
for i = 1:numel(MarkerNames)
    MarkerData(:,(i-1)*3+(1:3)) = MarkerData(:,(i-1)*3+abs(AxisTransform)).*(ones(nframes,1)*sign(AxisTransform));
end

%% Marker Indexing
%Determine index of markers within TRC file
HeelMarkerIndex = find(strcmp(MarkerNames,HeelMarker),1);
ToeMarkerIndex = find(strcmp(MarkerNames,ToeMarker),1);
ToeLateralIndex = find(strcmp(MarkerNames,ToeLateral),1);
ToeMedialIndex = find(strcmp(MarkerNames,ToeMedial),1);
MidfootSupIndex = find(strcmp(MarkerNames,MidfootSup),1);
MidfootLatIndex = find(strcmp(MarkerNames,MidfootLat),1);

WandIndex = cell(length(WandNames),1);
for i = 1:length(WandNames)
    WandIndex{i} = find(strcmp(MarkerNames,WandNames{i}),1);
end

%% Assign marker locations
%Determine marker locations
HeelMarkerLoc = mean(MarkerData(Range,(HeelMarkerIndex-1)*3+1:(HeelMarkerIndex-1)*3+3));
ToeMarkerLoc = mean(MarkerData(Range,(ToeMarkerIndex-1)*3+1:(ToeMarkerIndex-1)*3+3));
ToeLateralLoc = mean(MarkerData(Range,(ToeLateralIndex-1)*3+1:(ToeLateralIndex-1)*3+3));
ToeMedialLoc = mean(MarkerData(Range,(ToeMedialIndex-1)*3+1:(ToeMedialIndex-1)*3+3));
MidfootSupLoc = mean(MarkerData(Range,(MidfootSupIndex-1)*3+1:(MidfootSupIndex-1)*3+3));
MidfootLatLoc = mean(MarkerData(Range,(MidfootLatIndex-1)*3+1:(MidfootLatIndex-1)*3+3));

WandMarkerLocs = cell(length(WandNames),1);
for i = 1:length(WandNames)
    WandMarkerLocs{i} = MarkerData(:,(WandIndex{i}-1)*3+1:(WandIndex{i}-1)*3+3);
end

%Specify specific position of critical wand markers for calculating wand
%tip location
if WandTip > 0
    WandEndLoc = WandMarkerLocs{WandEnd};
end

WandTipLoc = WandMarkerLocs{WandTip};

if WandMid1 > 0 && WandMid2>0 && WandEnd == 0
    WandEndLoc = mean(cat(3,WandMarkerLocs{WandMid1},WandMarkerLocs{WandMid2}),3);
    WandMidLoc = 0;
    NoMid = 1;
elseif WandMid1 > 0 && WandMid2>0
    WandMidLoc = mean(cat(3,WandMarkerLocs{WandMid1},WandMarkerLocs{WandMid2}),3);
    NoMid = 0;
elseif WandMid1 > 0
    WandMidLoc = WandMarkerLocs{WandMid1};
    NoMid = 0;
elseif WandMid2 > 0
    WandMidLoc = WandMarkerLocs{WandMid2};
    NoMid = 0;
else
    WandMidLoc = 0;
    NoMid = 1;
end

%% Find foot trace

if ~NoMid
    
    WandLine(:,:,1) = (WandTipLoc-WandEndLoc)./mean(sqrt(sum((WandTipLoc-WandEndLoc).^2,2)));
    WandLine(:,:,2) = (WandTipLoc-WandMidLoc)./mean(sqrt(sum((WandTipLoc-WandMidLoc).^2,2)));
    WandLine(:,:,3) = (WandMidLoc-WandEndLoc)./mean(sqrt(sum((WandMidLoc-WandEndLoc).^2,2)));
    
    WandLine = mean(WandLine,3);
else
    WandLine(:,:,1) = (WandTipLoc-WandEndLoc)./mean(sqrt(sum((WandTipLoc-WandEndLoc).^2,2)));
end

optParams.WandTipLoc = WandTipLoc;
optParams.WandMidLoc = WandMidLoc;
optParams.WandEndLoc = WandEndLoc;
optParams.NoMid = NoMid;
optParams.nframes = nframes;
optParams.WandLine = WandLine;

options = optimset('display','iter');
TipDist = lsqnonlin(@costFuncOutline, [.085 .2 .25],[],[],options,optParams)

% TipDist = [.085 .2 .25];

if ~NoMid
    
    TipLoc(:,:,1) = TipDist(1)*WandLine+WandTipLoc;
    TipLoc(:,:,2) = TipDist(2)*WandLine+WandMidLoc;
    TipLoc(:,:,3) = TipDist(3)*WandLine+WandEndLoc;
    
else
    
    TipLoc(:,:,1) = TipDist(1)*WandLine+WandTipLoc;
    TipLoc(:,:,2) = TipDist(2)*WandLine+WandEndLoc;
    
end

TipLoc = TipLoc(:,:,1);%mean(TipLoc,3);

% Remove time frames were wand is not in contact with ground
if strcmp('r',lower(FootSide))
    TipLoc(TipLoc(:,2)>0.005,:) = [];
elseif strcmp('l',lower(FootSide))
    TipLoc(TipLoc(:,2)>0.005,:) = [];
end

%% Sort wand tip locations based on angle away from MidFoot axis make consistent polygon

MidFootLoc = (HeelMarkerLoc+ToeMarkerLoc+ToeMedialLoc+ToeLateralLoc)/4;

MidFootLoc(2) = 0;
% ToeMarkerLoc(2) = 0;
% HeelMarkerLoc(2) = 0;
% ToeMedialLoc(2) = 0;
% ToeLateralLoc(2) = 0;

MidAxis = ToeMarkerLoc-MidFootLoc;
MidAxis = MidAxis/norm(MidAxis);
nframes = length(TipLoc);
MidVecs = (TipLoc-ones(nframes,1)*MidFootLoc);
MidVecs = MidVecs./((sum(MidVecs.^2,2).^0.5)*ones(1,3));

crossProd = sum(cross(MidVecs,ones(nframes,1)*MidAxis,2),2);
dotProd = dot(MidVecs,ones(nframes,1)*MidAxis,2);

Angles = atan2(crossProd,dotProd);

% figure
% plot(Angles)

[AngleSorted,order] = sort(Angles);
TipLoc = TipLoc(order,:);
% TipLoc(:,1) = csaps(AngleSorted,TipLoc(:,1),.99,AngleSorted);
% TipLoc(:,3) = csaps(AngleSorted,TipLoc(:,3),.99,AngleSorted);

%% Filter Tip Locations to get smooth outline curve

TraceIndex = linspace(0,1,length(TipLoc))';

[~,TipLoc] = spaps(TraceIndex,TipLoc',0.000007,1);
TipLoc = TipLoc';

%% Transform data to use foot trace principle axes as x and y coordinates
TipLoc = TipLoc(:,1:3);

covariance = cov(TipLoc(:,[1 3]));
[coeffs,D] = eig(covariance);
coeffs = -coeffs(:,[2 1]);

% coeffs = pca(TipLoc(:,[1 3]));
TipLoc(:,[1 3]) = TipLoc(:,[1 3])*coeffs;

MinTipLoc = min(TipLoc);
MaxTipLoc = max(TipLoc);

HeelMarkerLoc([1 3]) = HeelMarkerLoc([1 3])*coeffs;
ToeMarkerLoc([1 3]) = ToeMarkerLoc([1 3])*coeffs;
ToeLateralLoc([1 3]) = ToeLateralLoc([1 3])*coeffs;
ToeMedialLoc([1 3]) = ToeMedialLoc([1 3])*coeffs;
MidfootSupLoc([1 3]) = MidfootSupLoc([1 3])*coeffs;
MidfootLatLoc([1 3]) = MidfootLatLoc([1 3])*coeffs;

TipLoc(:,[1 3]) = TipLoc(:,[1 3])-ones(length(TipLoc),1)*MinTipLoc(:,[1 3]);
HeelMarkerLoc(:,[1 3]) = HeelMarkerLoc(:,[1 3])-MinTipLoc(:,[1 3]);
ToeMarkerLoc(:,[1 3]) = ToeMarkerLoc(:,[1 3])-MinTipLoc(:,[1 3]);
ToeLateralLoc(:,[1 3]) = ToeLateralLoc(:,[1 3])-MinTipLoc(:,[1 3]);
ToeMedialLoc(:,[1 3]) = ToeMedialLoc(:,[1 3])-MinTipLoc(:,[1 3]);
MidfootSupLoc(:,[1 3]) = MidfootSupLoc(:,[1 3])-MinTipLoc(:,[1 3]);
MidfootLatLoc(:,[1 3]) = MidfootLatLoc(:,[1 3])-MinTipLoc(:,[1 3]);

%% Define foot axes

MidFootLoc = (HeelMarkerLoc+ToeMarkerLoc+ToeMedialLoc+ToeLateralLoc)/4;

MidFootLoc(2) = 0;
% ToeMarkerLoc(2) = 0;
% HeelMarkerLoc(2) = 0;
% ToeMedialLoc(2) = 0;
% ToeLateralLoc(2) = 0;

% Determine normalized axes
LengthAxis = [1 0 0];
WidthAxis = [0 0 1];
MidAxis = ToeMarkerLoc-MidFootLoc;
MidAxis = MidAxis/norm(MidAxis);
ToesAxisApprox = ToeLateralLoc-ToeMedialLoc;

%% Plot foot outline

figure
plot(TipLoc(:,1),TipLoc(:,3),'kx-')
hold on
plot([HeelMarkerLoc(1) ToeMarkerLoc(1)],[HeelMarkerLoc(3) ToeMarkerLoc(3)],'mo-')
plot([ToeLateralLoc(1) ToeMedialLoc(1)],[ToeLateralLoc(3) ToeMedialLoc(3)],'mo-')
plot(MidFootLoc(1), MidFootLoc(3),'mo')
plot(ToeLateralLoc(1),ToeLateralLoc(3),'rx','LineWidth',3)
plot(ToeMedialLoc(1),ToeMedialLoc(3),'bx','LineWidth',3)
plot(HeelMarkerLoc(1),HeelMarkerLoc(3),'rx','LineWidth',3)
plot(ToeMarkerLoc(1),ToeMarkerLoc(3),'bx','LineWidth',3)

%% Create spring locations
% creates a rectangular grid on the foot using axes drawn from the heel to
% the toe and from the medial toe to lateral toe

offset = -0.005;
% Define grid spread along foot
if strcmp('r',lower(FootSide))
    SpringLocsML = linspace(-offset,MaxTipLoc(3)-MinTipLoc(3)+offset,numWidth)';
    SpringLocsAP = linspace(-offset,MaxTipLoc(1)-MinTipLoc(1)+offset,numLength)';
elseif strcmp('l',lower(FootSide))
    SpringLocsML = linspace(-offset,MaxTipLoc(3)-MinTipLoc(3)+offset,numWidth)';
    SpringLocsAP = linspace(-offset,MaxTipLoc(1)-MinTipLoc(1)+offset,numLength)';
else
    error('\nSpecify foot side as right or left.\n')
end

SpringLocs = zeros(numLength*numWidth,3);
% Create spring locations by row
for i = 1:numLength
    SpringLocs((i-1)*numWidth+1:i*numWidth,:) = SpringLocsML*WidthAxis...
        +ones(numWidth,1)*(SpringLocsAP(i,:)*LengthAxis);
end

plot(SpringLocs(:,1), SpringLocs(:,3),'dg')

%% Determine which springs are within foot trace region

% Find which springs are within polygon from trace
inSprings = inpolygon(SpringLocs(:,1),SpringLocs(:,3),TipLoc(:,1),TipLoc(:,3));

SpringLocs = SpringLocs(inSprings,:);

plot(SpringLocs(:,1), SpringLocs(:,3),'dr')

numSprings = size(SpringLocs,1);

%% Define orientation of shoe trace in opensim bodies
% define positions of markers in body frames from osim file for use in
% calculating rotation matrices
if strcmp('r',lower(FootSide))
    ToeJoint = model.getBodySet.get(String('toes_r')).getJoint;
    ToeAxisLoc = char(ToeJoint.get_location_in_parent);
    ToeAxisLoc = str2num(ToeAxisLoc(2:end));
    ToeAxisOrientation = char(ToeJoint.get_orientation_in_parent); % Angle of rotation for joint in parent body in opensim
    ToeAxisOrientation = str2num(ToeAxisOrientation(2:end));
    ToeRotationAxis = [0 0 1]; % rotates about z axis
    MarkerPlateLoc = char(model.getBodySet.get(String('calcn_r_markerPlate')).getJoint.get_location_in_parent);
    MarkerPlateLoc = str2num(MarkerPlateLoc(2:end));
    HeelMarkerB = char(model.getMarkerSet.get(String(HeelMarker)).getOffset);
    HeelMarkerB = str2num(HeelMarkerB(2:end));
    MidfootLatB = char(model.getMarkerSet.get(String(MidfootLat)).getOffset);
    MidfootLatB = str2num(MidfootLatB(2:end));
    MidfootSupB = char(model.getMarkerSet.get(String(MidfootSup)).getOffset);
    MidfootSupB = str2num(MidfootSupB(2:end));
    ToeMarkerB = char(model.getMarkerSet.get(String(ToeMarker)).getOffset);
    ToeMarkerB = str2num(ToeMarkerB(2:end));
    ToeLateralB = char(model.getMarkerSet.get(String(ToeLateral)).getOffset);
    ToeLateralB = str2num(ToeLateralB(2:end));
    ToeMedialB = char(model.getMarkerSet.get(String(ToeMedial)).getOffset);
    ToeMedialB = str2num(ToeMedialB(2:end));
elseif strcmp('l',lower(FootSide))
    ToeJoint = model.getBodySet.get(String('toes_l')).getJoint;
    ToeAxisLoc = char(ToeJoint.get_location_in_parent);
    ToeAxisLoc = str2num(ToeAxisLoc(2:end));
    ToeAxisOrientation = char(ToeJoint.get_orientation_in_parent); % Angle of rotation for joint in parent body in opensim
    ToeAxisOrientation = str2num(ToeAxisOrientation(2:end));
    ToeRotationAxis = [0 0 -1]; % rotates about z axis
    MarkerPlateLoc = char(model.getBodySet.get(String('calcn_l_markerPlate')).getJoint.get_location_in_parent);
    MarkerPlateLoc = str2num(MarkerPlateLoc(2:end));
    HeelMarkerB = char(model.getMarkerSet.get(String(HeelMarker)).getOffset);
    HeelMarkerB = str2num(HeelMarkerB(2:end));
    MidfootLatB = char(model.getMarkerSet.get(String(MidfootLat)).getOffset);
    MidfootLatB = str2num(MidfootLatB(2:end));
    MidfootSupB = char(model.getMarkerSet.get(String(MidfootSup)).getOffset);
    MidfootSupB = str2num(MidfootSupB(2:end));
    ToeMarkerB = char(model.getMarkerSet.get(String(ToeMarker)).getOffset);
    ToeMarkerB = str2num(ToeMarkerB(2:end));
    ToeLateralB = char(model.getMarkerSet.get(String(ToeLateral)).getOffset);
    ToeLateralB = str2num(ToeLateralB(2:end));
    ToeMedialB = char(model.getMarkerSet.get(String(ToeMedial)).getOffset);
    ToeMedialB = str2num(ToeMedialB(2:end));
else
    error('\nSpecify foot side as right or left.\n')
end

% Use triad method to determine rotation matrix from global frame to body frame
r1Heel = MidfootLatLoc-HeelMarkerLoc;
r2Heel = MidfootSupLoc-HeelMarkerLoc;
R1Heel = MidfootLatB-HeelMarkerB;
R2Heel = MidfootSupB-HeelMarkerB;

mHeel = cross(r1Heel,r2Heel);
MHeel = cross(R1Heel,R2Heel);

r2Heel = cross(mHeel,r1Heel);
R2Heel = cross(MHeel,R1Heel);

r1Heel = (r1Heel/norm(r1Heel))';
r2Heel = (r2Heel/norm(r2Heel))';
R1Heel = (R1Heel/norm(R1Heel))';
R2Heel = (R2Heel/norm(R2Heel))';
mHeel = (mHeel/norm(mHeel))';
MHeel = (MHeel/norm(MHeel))';

Rgroundtoheel = [R1Heel R2Heel MHeel]*[r1Heel r2Heel mHeel]';
% Rgroundtotoe = Rgroundtoheel;%*[1 0 0; 0 1 0; ToeRotationAxis];

r1Toe = ToeLateralLoc-ToeMarkerLoc;
r2Toe = ToeMedialLoc-ToeMarkerLoc;
R1Toe = ToeLateralB-ToeMarkerB;
R2Toe = ToeMedialB-ToeMarkerB;

mToe = cross(r1Toe,r2Toe);
MToe = cross(R1Toe,R2Toe);

r2Toe = cross(mToe,r1Toe);
R2Toe = cross(MToe,R1Toe);

r1Toe = (r1Toe/norm(r1Toe))';
r2Toe = (r2Toe/norm(r2Toe))';
R1Toe = (R1Toe/norm(R1Toe))';
R2Toe = (R2Toe/norm(R2Toe))';
mToe = (mToe/norm(mToe))';
MToe = (MToe/norm(MToe))';

Rgroundtotoe = [R1Toe R2Toe MToe]*[r1Toe r2Toe mToe]';

% Determine position of toe axis origin in foot trace frame
ToeAxisLocTrace = Rgroundtoheel'*(ToeAxisLoc-MarkerPlateLoc-HeelMarkerB)'+HeelMarkerLoc';

% Determine orientation of toe axis in foot trace frame
AxisRotInOpenSimX = [1 0 0;...
    0 cos(ToeAxisOrientation(1)) -sin(ToeAxisOrientation(1));...
    0 sin(ToeAxisOrientation(1)) cos(ToeAxisOrientation(1))];
AxisRotInOpenSimY = [cos(ToeAxisOrientation(2)) 0 sin(ToeAxisOrientation(2));...
    0 1 0;...
    -sin(ToeAxisOrientation(2)) 0 cos(ToeAxisOrientation(2))];
AxisRotInOpenSimZ = [cos(ToeAxisOrientation(3)) -sin(ToeAxisOrientation(3)) 0;...
    sin(ToeAxisOrientation(3)) cos(ToeAxisOrientation(3)) 0;
    0 0 1];
ToeAxisRotMatOsim = AxisRotInOpenSimX*AxisRotInOpenSimY*AxisRotInOpenSimZ;
ToeRotationAxisTrace = Rgroundtoheel'*ToeAxisRotMatOsim*ToeRotationAxis'

%plot toe rotational axis
plot([ToeAxisLocTrace(1)+ToeRotationAxisTrace(1)/10; ToeAxisLocTrace(1); ToeAxisLocTrace(1)-ToeRotationAxisTrace(1)/10],...
    [ToeAxisLocTrace(3)+ToeRotationAxisTrace(3)/10; ToeAxisLocTrace(3); ToeAxisLocTrace(3)-ToeRotationAxisTrace(3)/10],'cx-')

%% Determine which segments contain springs using dot product to compute angles

% Define point on toe axis for determining spring segments
ToeAxisPoint = (ToeAxisLocTrace+ToeRotationAxisTrace/10)';

SpringVec = (SpringLocs-ones(numSprings,1)*ToeAxisPoint)./(sum((SpringLocs-ones(numSprings,1)*ToeAxisPoint).^2,2).^0.5*ones(1,3));
Angle = asin(cross(SpringVec,ones(numSprings,1)*ToeRotationAxisTrace',2));

Angle = Angle(:,2);

if strcmp('r',lower(FootSide))
    SpringsH = SpringLocs(Angle>0,:);
    SpringsT = SpringLocs(Angle<0,:);
elseif strcmp('l',lower(FootSide))
    SpringsH = SpringLocs(Angle<0,:);
    SpringsT = SpringLocs(Angle>0,:);
else
    error('\nSpecify foot side as right or left.\n')
end

plot(SpringsT(:,1),SpringsT(:,3),'db')

if strcmpi('r',FootSide)
	save spring_pos_for_polynomial_r.mat SpringsH SpringsT
elseif strcmpi('l',FootSide)
    save spring_pos_for_polynomial_l.mat SpringsH SpringsT
end

%% Transform spring locations into foot and toes frame

% transform spring positions into body frame
numHeel = length(SpringsH)
for i = 1:numHeel
    Position = SpringsH(i,:)-HeelMarkerLoc;
    Position = Rgroundtoheel*Position'+HeelMarkerB';
    SpringsH(i,:) = Position';
end

numToe = length(SpringsT)
for i = 1:numToe
    Position = SpringsT(i,:)-ToeMarkerLoc;
    Position = Rgroundtotoe*Position'+ToeMarkerB';
    SpringsT(i,:) = Position';
end

% figure
% % plot(SpringsH(:,1), SpringsH(:,3),'dg')
% plot(SpringsT(:,1), SpringsT(:,3),'dg')
% hold on
% plot([ToeMarkerB(1)],[ToeMarkerB(3)],'mo-')
% plot([ToeLateralB(1) ToeMedialB(1)],[ToeLateralB(3) ToeMedialB(3)],'mo-')
% 
% figure
% plot(SpringsH(:,1), SpringsH(:,3),'dg')
% hold on
% plot([HeelMarkerB(1)],[HeelMarkerB(3)],'mo-')
% plot([MidfootLatB(1) MidfootSupB(1)],[MidfootLatB(3) MidfootSupB(3)],'mo-')

if strcmp('r',lower(FootSide))
    save SpringLocsR.mat SpringsH SpringsT
elseif strcmp('l',lower(FootSide))
    save SpringLocsL.mat SpringsH SpringsT
else
    error('\nSpecify foot side as right or left.\n')
end

%% Now do the other foot

if strcmp('r',lower(FootSide))
    FootSide = 'L';
elseif strcmp('l',lower(FootSide))
    FootSide = 'R';
end 

TraceTRC = ['Patient4_Wand' FootSide '_mokka.trc'];

Range = 500:1000; % Time range to determine foot experimental marker locations

WandNames = {'Wand.End'
    'Wand.Tip'
    'Wand.Mid1'
    'Wand.Mid2'};

% Specify index of marker in names listed above
%NOTE If midpoint markers are not on the axis of the wand, then two must be
%specified
WandEnd = 1;
WandTip = 2; %optional if  2 midpoint markers are specified
WandMid1 = 0; %optional, enter 0 to ignore
WandMid2 = 0; %optional, enter 0 to ignore

ToeMarker = [FootSide '.Toe'];
HeelMarker = [FootSide '.Heel'];
ToeLateral = [FootSide '.Toe.Lateral'];
ToeMedial = [FootSide '.Toe.Medial'];
MidfootSup = [FootSide '.Midfoot.Superior'];
MidfootLat = [FootSide '.Midfoot.Lateral'];

%% Load data
% Load TRC file with marker data
[MarkerData, MarkerNames, SampleRate, Time] = LoadTRCFile(TraceTRC);

nframes = numel(Time);
MarkerData = MarkerData/1000;

MarkerDataOrig = MarkerData;

%% Transform Data
%Transform marker data
for i = 1:numel(MarkerNames)
    MarkerData(:,(i-1)*3+(1:3)) = MarkerData(:,(i-1)*3+abs(AxisTransform)).*(ones(nframes,1)*sign(AxisTransform));
end

%% Marker Indexing
%Determine index of markers within TRC file
HeelMarkerIndex = find(strcmp(MarkerNames,HeelMarker),1);
ToeMarkerIndex = find(strcmp(MarkerNames,ToeMarker),1);
ToeLateralIndex = find(strcmp(MarkerNames,ToeLateral),1);
ToeMedialIndex = find(strcmp(MarkerNames,ToeMedial),1);
MidfootSupIndex = find(strcmp(MarkerNames,MidfootSup),1);
MidfootLatIndex = find(strcmp(MarkerNames,MidfootLat),1);

WandIndex = cell(length(WandNames),1);
for i = 1:length(WandNames)
    WandIndex{i} = find(strcmp(MarkerNames,WandNames{i}),1);
end

%% Assign marker locations
%Determine marker locations
HeelMarkerLoc = mean(MarkerData(Range,(HeelMarkerIndex-1)*3+1:(HeelMarkerIndex-1)*3+3));
ToeMarkerLoc = mean(MarkerData(Range,(ToeMarkerIndex-1)*3+1:(ToeMarkerIndex-1)*3+3));
ToeLateralLoc = mean(MarkerData(Range,(ToeLateralIndex-1)*3+1:(ToeLateralIndex-1)*3+3));
ToeMedialLoc = mean(MarkerData(Range,(ToeMedialIndex-1)*3+1:(ToeMedialIndex-1)*3+3));
MidfootSupLoc = mean(MarkerData(Range,(MidfootSupIndex-1)*3+1:(MidfootSupIndex-1)*3+3));
MidfootLatLoc = mean(MarkerData(Range,(MidfootLatIndex-1)*3+1:(MidfootLatIndex-1)*3+3));

WandMarkerLocs = cell(length(WandNames),1);
for i = 1:length(WandNames)
    WandMarkerLocs{i} = MarkerData(:,(WandIndex{i}-1)*3+1:(WandIndex{i}-1)*3+3);
end

%Specify specific position of critical wand markers for calculating wand
%tip location
if WandTip > 0
    WandEndLoc = WandMarkerLocs{WandEnd};
end

WandTipLoc = WandMarkerLocs{WandTip};

if WandMid1 > 0 && WandMid2>0 && WandEnd == 0
    WandEndLoc = mean(cat(3,WandMarkerLocs{WandMid1},WandMarkerLocs{WandMid2}),3);
    WandMidLoc = 0;
    NoMid = 1;
elseif WandMid1 > 0 && WandMid2>0
    WandMidLoc = mean(cat(3,WandMarkerLocs{WandMid1},WandMarkerLocs{WandMid2}),3);
    NoMid = 0;
elseif WandMid1 > 0
    WandMidLoc = WandMarkerLocs{WandMid1};
    NoMid = 0;
elseif WandMid2 > 0
    WandMidLoc = WandMarkerLocs{WandMid2};
    NoMid = 0;
else
    WandMidLoc = 0;
    NoMid = 1;
end

%% Find foot trace

clear WandLine

if ~NoMid
    
    WandLine(:,:,1) = (WandTipLoc-WandEndLoc)./mean(sqrt(sum((WandTipLoc-WandEndLoc).^2,2)));
    WandLine(:,:,2) = (WandTipLoc-WandMidLoc)./mean(sqrt(sum((WandTipLoc-WandMidLoc).^2,2)));
    WandLine(:,:,3) = (WandMidLoc-WandEndLoc)./mean(sqrt(sum((WandMidLoc-WandEndLoc).^2,2)));
    
    WandLine = mean(WandLine,3);
else
    WandLine(:,:,1) = (WandTipLoc-WandEndLoc)./mean(sqrt(sum((WandTipLoc-WandEndLoc).^2,2)));
end

optParams.WandTipLoc = WandTipLoc;
optParams.WandMidLoc = WandMidLoc;
optParams.WandEndLoc = WandEndLoc;
optParams.NoMid = NoMid;
optParams.nframes = nframes;
optParams.WandLine = WandLine;

options = optimset('display','iter');
TipDist = lsqnonlin(@costFuncOutline, [.085 .2 .25],[],[],options,optParams)

% TipDist = [.085 .2 .25];

TipLocOtherFoot = TipLoc;
clear TipLoc

if ~NoMid
    
    TipLoc(:,:,1) = TipDist(1)*WandLine+WandTipLoc;
    TipLoc(:,:,2) = TipDist(2)*WandLine+WandMidLoc;
    TipLoc(:,:,3) = TipDist(3)*WandLine+WandEndLoc;
    
else
    
    TipLoc(:,:,1) = TipDist(1)*WandLine+WandTipLoc;
    TipLoc(:,:,2) = TipDist(2)*WandLine+WandEndLoc;
    
end

TipLoc = TipLoc(:,:,1);%mean(TipLoc,3);

% Remove time frames were wand is not in contact with ground
if strcmp('r',lower(FootSide))
    TipLoc(TipLoc(:,2)>0.005,:) = [];
elseif strcmp('l',lower(FootSide))
    TipLoc(TipLoc(:,2)>0.005,:) = [];
end

%% Sort wand tip locations based on angle away from MidFoot axis make consistent polygon

MidFootLoc = (HeelMarkerLoc+ToeMarkerLoc+ToeMedialLoc+ToeLateralLoc)/4;

MidFootLoc(2) = 0;
% ToeMarkerLoc(2) = 0;
% HeelMarkerLoc(2) = 0;
% ToeMedialLoc(2) = 0;
% ToeLateralLoc(2) = 0;

MidAxis = ToeMarkerLoc-MidFootLoc;
MidAxis = MidAxis/norm(MidAxis);
nframes = length(TipLoc);
MidVecs = (TipLoc-ones(nframes,1)*MidFootLoc);
MidVecs = MidVecs./((sum(MidVecs.^2,2).^0.5)*ones(1,3));

crossProd = sum(cross(MidVecs,ones(nframes,1)*MidAxis,2),2);
dotProd = dot(MidVecs,ones(nframes,1)*MidAxis,2);

Angles = atan2(crossProd,dotProd);

% figure
% plot(Angles)

[AngleSorted,order] = sort(Angles);
TipLoc = TipLoc(order,:);
% TipLoc(:,1) = csaps(AngleSorted,TipLoc(:,1),.99,AngleSorted);
% TipLoc(:,3) = csaps(AngleSorted,TipLoc(:,3),.99,AngleSorted);

%% Filter Tip Locations to get smooth outline curve

TraceIndex = linspace(0,1,length(TipLoc))';

[~,TipLoc] = spaps(TraceIndex,TipLoc',0.00005,1);
TipLoc = TipLoc';

%% Transform data to use foot trace principle axes as x and y coordinates
TipLoc = TipLoc(:,1:3);

covariance = cov(TipLoc(:,[1 3]));
[coeffs,D] = eig(covariance);
coeffs = -coeffs(:,[2 1]);

% coeffs = pca(TipLoc(:,[1 3]));
TipLoc(:,[1 3]) = TipLoc(:,[1 3])*coeffs;
HeelMarkerLoc([1 3]) = HeelMarkerLoc([1 3])*coeffs;
ToeMarkerLoc([1 3]) = ToeMarkerLoc([1 3])*coeffs;
ToeLateralLoc([1 3]) = ToeLateralLoc([1 3])*coeffs;
ToeMedialLoc([1 3]) = ToeMedialLoc([1 3])*coeffs;
MidfootSupLoc([1 3]) = MidfootSupLoc([1 3])*coeffs;
MidfootLatLoc([1 3]) = MidfootLatLoc([1 3])*coeffs;

MinTipLoc = min(TipLoc);
MaxTipLoc = max(TipLoc);

TipLoc(:,[1 3]) = TipLoc(:,[1 3])-ones(length(TipLoc),1)*MinTipLoc(:,[1 3]);
HeelMarkerLoc(:,[1 3]) = HeelMarkerLoc(:,[1 3])-MinTipLoc(:,[1 3]);
ToeMarkerLoc(:,[1 3]) = ToeMarkerLoc(:,[1 3])-MinTipLoc(:,[1 3]);
ToeLateralLoc(:,[1 3]) = ToeLateralLoc(:,[1 3])-MinTipLoc(:,[1 3]);
ToeMedialLoc(:,[1 3]) = ToeMedialLoc(:,[1 3])-MinTipLoc(:,[1 3]);
MidfootSupLoc(:,[1 3]) = MidfootSupLoc(:,[1 3])-MinTipLoc(:,[1 3]);
MidfootLatLoc(:,[1 3]) = MidfootLatLoc(:,[1 3])-MinTipLoc(:,[1 3]);

%% Define foot axes

MidFootLoc = (HeelMarkerLoc+ToeMarkerLoc+ToeMedialLoc+ToeLateralLoc)/4;

MidFootLoc(2) = 0;
% ToeMarkerLoc(2) = 0;
% HeelMarkerLoc(2) = 0;
% ToeMedialLoc(2) = 0;
% ToeLateralLoc(2) = 0;

% Determine normalized axes
LengthAxis = [1 0 0];
WidthAxis = [0 0 1];
MidAxis = ToeMarkerLoc-MidFootLoc;
MidAxis = MidAxis/norm(MidAxis);

ToesAxisApprox = ToeLateralLoc-ToeMedialLoc;

%% Plot foot outline

TipLoc = TipLocOtherFoot;
TipLoc(:,3) = -TipLoc(:,3)+max(TipLoc(:,3));

figure
plot(TipLoc(:,1),TipLoc(:,3),'kx-')
hold on
plot([HeelMarkerLoc(1) ToeMarkerLoc(1)],[HeelMarkerLoc(3) ToeMarkerLoc(3)],'mo-')
plot([ToeLateralLoc(1) ToeMedialLoc(1)],[ToeLateralLoc(3) ToeMedialLoc(3)],'mo-')
plot(MidFootLoc(1), MidFootLoc(3),'mo')
plot(ToeLateralLoc(1),ToeLateralLoc(3),'rx','LineWidth',3)
plot(ToeMedialLoc(1),ToeMedialLoc(3),'bx','LineWidth',3)
plot(HeelMarkerLoc(1),HeelMarkerLoc(3),'rx','LineWidth',3)
plot(ToeMarkerLoc(1),ToeMarkerLoc(3),'bx','LineWidth',3)

%% Create spring locations
% creates a rectangular grid on the foot using axes drawn from the heel to
% the toe and from the medial toe to lateral toe

SpringLocs = zeros(numLength*numWidth,3);
% Create spring locations by row
for i = 1:numLength
    SpringLocs((i-1)*numWidth+1:i*numWidth,:) = SpringLocsML*WidthAxis...
        +ones(numWidth,1)*(SpringLocsAP(i,:)*LengthAxis);
end

plot(SpringLocs(:,1), SpringLocs(:,3),'dg')

%% Determine which springs are within foot trace region

% Find which springs are within polygon from trace
inSprings = inpolygon(SpringLocs(:,1),SpringLocs(:,3),TipLoc(:,1),TipLoc(:,3));

SpringLocs = SpringLocs(inSprings,:);

plot(SpringLocs(:,1), SpringLocs(:,3),'dr')

numSprings = size(SpringLocs,1);

%% Define orientation of shoe trace in opensim bodies
% define positions of markers in body frames from osim file for use in
% calculating rotation matrices
if strcmp('r',lower(FootSide))
    ToeJoint = model.getBodySet.get(String('toes_r')).getJoint;
    ToeAxisLoc = char(ToeJoint.get_location_in_parent);
    ToeAxisLoc = str2num(ToeAxisLoc(2:end));
    ToeAxisOrientation = char(ToeJoint.get_orientation_in_parent); % Angle of rotation for joint in parent body in opensim
    ToeAxisOrientation = str2num(ToeAxisOrientation(2:end));
    ToeRotationAxis = [0 0 1]; % rotates about z axis
    MarkerPlateLoc = char(model.getBodySet.get(String('calcn_r_markerPlate')).getJoint.get_location_in_parent);
    MarkerPlateLoc = str2num(MarkerPlateLoc(2:end));
    HeelMarkerB = char(model.getMarkerSet.get(String(HeelMarker)).getOffset);
    HeelMarkerB = str2num(HeelMarkerB(2:end));
    MidfootLatB = char(model.getMarkerSet.get(String(MidfootLat)).getOffset);
    MidfootLatB = str2num(MidfootLatB(2:end));
    MidfootSupB = char(model.getMarkerSet.get(String(MidfootSup)).getOffset);
    MidfootSupB = str2num(MidfootSupB(2:end));
    ToeMarkerB = char(model.getMarkerSet.get(String(ToeMarker)).getOffset);
    ToeMarkerB = str2num(ToeMarkerB(2:end));
    ToeLateralB = char(model.getMarkerSet.get(String(ToeLateral)).getOffset);
    ToeLateralB = str2num(ToeLateralB(2:end));
    ToeMedialB = char(model.getMarkerSet.get(String(ToeMedial)).getOffset);
    ToeMedialB = str2num(ToeMedialB(2:end));
elseif strcmp('l',lower(FootSide))
    ToeJoint = model.getBodySet.get(String('toes_l')).getJoint;
    ToeAxisLoc = char(ToeJoint.get_location_in_parent);
    ToeAxisLoc = str2num(ToeAxisLoc(2:end));
    ToeAxisOrientation = char(ToeJoint.get_orientation_in_parent); % Angle of rotation for joint in parent body in opensim
    ToeAxisOrientation = str2num(ToeAxisOrientation(2:end));
    ToeRotationAxis = [0 0 -1]; % rotates about z axis
    MarkerPlateLoc = char(model.getBodySet.get(String('calcn_l_markerPlate')).getJoint.get_location_in_parent);
    MarkerPlateLoc = str2num(MarkerPlateLoc(2:end));
    HeelMarkerB = char(model.getMarkerSet.get(String(HeelMarker)).getOffset);
    HeelMarkerB = str2num(HeelMarkerB(2:end));
    MidfootLatB = char(model.getMarkerSet.get(String(MidfootLat)).getOffset);
    MidfootLatB = str2num(MidfootLatB(2:end));
    MidfootSupB = char(model.getMarkerSet.get(String(MidfootSup)).getOffset);
    MidfootSupB = str2num(MidfootSupB(2:end));
    ToeMarkerB = char(model.getMarkerSet.get(String(ToeMarker)).getOffset);
    ToeMarkerB = str2num(ToeMarkerB(2:end));
    ToeLateralB = char(model.getMarkerSet.get(String(ToeLateral)).getOffset);
    ToeLateralB = str2num(ToeLateralB(2:end));
    ToeMedialB = char(model.getMarkerSet.get(String(ToeMedial)).getOffset);
    ToeMedialB = str2num(ToeMedialB(2:end));
else
    error('\nSpecify foot side as right or left.\n')
end

% Use triad method to determine rotation matrix from global frame to body frame
r1Heel = MidfootLatLoc-HeelMarkerLoc;
r2Heel = MidfootSupLoc-HeelMarkerLoc;
R1Heel = MidfootLatB-HeelMarkerB;
R2Heel = MidfootSupB-HeelMarkerB;

mHeel = cross(r1Heel,r2Heel);
MHeel = cross(R1Heel,R2Heel);

r2Heel = cross(mHeel,r1Heel);
R2Heel = cross(MHeel,R1Heel);

r1Heel = (r1Heel/norm(r1Heel))';
r2Heel = (r2Heel/norm(r2Heel))';
R1Heel = (R1Heel/norm(R1Heel))';
R2Heel = (R2Heel/norm(R2Heel))';
mHeel = (mHeel/norm(mHeel))';
MHeel = (MHeel/norm(MHeel))';

Rgroundtoheel = [R1Heel R2Heel MHeel]*[r1Heel r2Heel mHeel]';
% Rgroundtotoe = Rgroundtoheel;%*[1 0 0; 0 1 0; ToeRotationAxis];

r1Toe = ToeLateralLoc-ToeMarkerLoc;
r2Toe = ToeMedialLoc-ToeMarkerLoc;
R1Toe = ToeLateralB-ToeMarkerB;
R2Toe = ToeMedialB-ToeMarkerB;

mToe = cross(r1Toe,r2Toe);
MToe = cross(R1Toe,R2Toe);

r2Toe = cross(mToe,r1Toe);
R2Toe = cross(MToe,R1Toe);

r1Toe = (r1Toe/norm(r1Toe))';
r2Toe = (r2Toe/norm(r2Toe))';
R1Toe = (R1Toe/norm(R1Toe))';
R2Toe = (R2Toe/norm(R2Toe))';
mToe = (mToe/norm(mToe))';
MToe = (MToe/norm(MToe))';

Rgroundtotoe = [R1Toe R2Toe MToe]*[r1Toe r2Toe mToe]';

% Determine position of toe axis origin in foot trace frame
ToeAxisLocTrace = Rgroundtoheel'*(ToeAxisLoc-MarkerPlateLoc-HeelMarkerB)'+HeelMarkerLoc';

% Determine orientation of toe axis in foot trace frame
AxisRotInOpenSimX = [1 0 0;...
    0 cos(ToeAxisOrientation(1)) -sin(ToeAxisOrientation(1));...
    0 sin(ToeAxisOrientation(1)) cos(ToeAxisOrientation(1))];
AxisRotInOpenSimY = [cos(ToeAxisOrientation(2)) 0 sin(ToeAxisOrientation(2));...
    0 1 0;...
    -sin(ToeAxisOrientation(2)) 0 cos(ToeAxisOrientation(2))];
AxisRotInOpenSimZ = [cos(ToeAxisOrientation(1)) -sin(ToeAxisOrientation(1)) 0;...
    sin(ToeAxisOrientation(1)) cos(ToeAxisOrientation(1)) 0;
    0 0 1];
ToeAxisRotMatOsim = AxisRotInOpenSimX*AxisRotInOpenSimY*AxisRotInOpenSimZ;
ToeRotationAxisTrace = Rgroundtoheel'*ToeAxisRotMatOsim*ToeRotationAxis'

%plot toe rotational axis
plot([ToeAxisLocTrace(1)+ToeRotationAxisTrace(1)/10; ToeAxisLocTrace(1); ToeAxisLocTrace(1)-ToeRotationAxisTrace(1)/10],...
    [ToeAxisLocTrace(3)+ToeRotationAxisTrace(3)/10; ToeAxisLocTrace(3); ToeAxisLocTrace(3)-ToeRotationAxisTrace(3)/10],'cx-')

%% Determine which segments contain springs using dot product to compute angles

% Define point on toe axis for determining spring segments
ToeAxisPoint = (ToeAxisLocTrace+ToeRotationAxisTrace/10)';

SpringVec = (SpringLocs-ones(numSprings,1)*ToeAxisPoint)./(sum((SpringLocs-ones(numSprings,1)*ToeAxisPoint).^2,2).^0.5*ones(1,3));
Angle = asin(cross(SpringVec,ones(numSprings,1)*ToeRotationAxisTrace',2));

Angle = Angle(:,2);

if strcmp('r',lower(FootSide))
    SpringsH = SpringLocs(Angle>0,:);
    SpringsT = SpringLocs(Angle<0,:);
elseif strcmp('l',lower(FootSide))
    SpringsH = SpringLocs(Angle<0,:);
    SpringsT = SpringLocs(Angle>0,:);
else
    error('\nSpecify foot side as right or left.\n')
end

plot(SpringsT(:,1),SpringsT(:,3),'db')

if strcmpi('r',FootSide)
	save spring_pos_for_polynomial_r.mat SpringsH SpringsT
elseif strcmpi('l',FootSide)
    save spring_pos_for_polynomial_l.mat SpringsH SpringsT
end

%% Transform spring locations into foot and toes frame

% transform spring positions into body frame
numHeel = length(SpringsH)
for i = 1:numHeel
    Position = SpringsH(i,:)-HeelMarkerLoc;
    Position = Rgroundtoheel*Position'+HeelMarkerB';
    SpringsH(i,:) = Position';
end

numToe = length(SpringsT)
for i = 1:numToe
    Position = SpringsT(i,:)-ToeMarkerLoc;
    Position = Rgroundtotoe*Position'+ToeMarkerB';
    SpringsT(i,:) = Position';
end

% figure
% % plot(SpringsH(:,1), SpringsH(:,3),'dg')
% plot(SpringsT(:,1), SpringsT(:,3),'dg')
% hold on
% plot([ToeMarkerB(1)],[ToeMarkerB(3)],'mo-')
% plot([ToeLateralB(1) ToeMedialB(1)],[ToeLateralB(3) ToeMedialB(3)],'mo-')
% 
% figure
% plot(SpringsH(:,1), SpringsH(:,3),'dg')
% hold on
% plot([HeelMarkerB(1)],[HeelMarkerB(3)],'mo-')
% plot([MidfootLatB(1) MidfootSupB(1)],[MidfootLatB(3) MidfootSupB(3)],'mo-')

if strcmp('r',lower(FootSide))
    save SpringLocsR.mat SpringsH SpringsT
elseif strcmp('l',lower(FootSide))
    save SpringLocsL.mat SpringsH SpringsT
else
    error('\nSpecify foot side as right or left.\n')
end

function cost = costFuncOutline(guess,optParams)
% Calculates the point distance from each marker to the wand tip by
% minimizing the height of the wand tip over all time frames.

WandTipLoc = optParams.WandTipLoc;
WandMidLoc = optParams.WandMidLoc;
WandEndLoc = optParams.WandEndLoc;
NoMid = optParams.NoMid;
nframes = optParams.nframes;
WandLine = optParams.WandLine;

if ~NoMid
    
    TipLoc(:,:,1) = guess(1)*WandLine+WandTipLoc;
    TipLoc(:,:,2) = guess(2)*WandLine+WandMidLoc;
    TipLoc(:,:,3) = guess(3)*WandLine+WandEndLoc;
    
else
    
    TipLoc(:,:,1) = guess(1)*WandLine+WandTipLoc;
    TipLoc(:,:,2) = guess(2)*WandLine+WandEndLoc;
    
end

cost = TipLoc(:,2,:);

function [MarkerData, MarkerNames, SampleRate, Time] = LoadTRCFile(inFile)

% Read 5 header lines from input .trc file
fid = fopen(inFile);
header1 = fgetl(fid);
header2 = fgetl(fid);
header3 = fgetl(fid);
header4 = fgetl(fid);
header5 = fgetl(fid);

% Read marker data as one long column
data = fscanf(fid, '%f');
fclose(fid);

% Reshape marker data based on number of rows and
% columns of marker data specified in the input data
% file, with two extra columns for frame and time
info = sscanf(header3,'%f %f %d %d');
nrows = info(3,1);
ncols = info(4,1)*3;
data = reshape(data, ncols+2, nrows)';

SampleRate = info(1,1);

MarkerNames = strrep(header4,'Patient 4:','');
MarkerNames = textscan(MarkerNames,'%s');
MarkerNames = MarkerNames{1};
MarkerNames(1:2) = [];

Time = data(:,2);
MarkerData = data(:,3:end);
