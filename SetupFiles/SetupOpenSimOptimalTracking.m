function auxdata = SetupOpenSimOptimalTracking(GPOPSResult)
% If GPOPSResult = 1, it will load the result of a GPOPS simulation instead
% of IK and ID results. This simulation result comes from folder 3_1,
% dynamically consistant IK-ID.

Initial_Guess = 1;

SpringConstants=1;

close all

%% Load experimental data and osim analysis data

% Specify names of files containing necessary data
MkrDataFile = 'Patient4_TMGait_0pt5_01_left_02.trc'; % Must be TRC file
GRFDataFile = 'Patient4_TMGait_0pt5_01_left_02_GRF.mot'; % Must be mot or sto file
IKDataFile = 'Patient4_TMGait_0pt5_01_left_02_ik.mot'; % Must be mot or sto file
IDDataFile = 'Patient4_TMGait_0pt5_01_left_02_id.sto'; % Must be mot or sto file

% Load data from files
[MkrData, MkrNames, MkrRate, MkrTime] = readTRCFile(MkrDataFile);
[GRFData, GRFLabels, GRFTime] = readMotFile(GRFDataFile);
[IKData, IKLabels, IKTime] = readMotFile(IKDataFile);
[IDData, ~, IDTime] = readMotFile(IDDataFile);

% Specify files containing faster speed data
MkrDataFile = 'Patient4_TMGait_0pt8_01_left_18.trc';
GRFDataFile = 'Patient4_TMGait_0pt8_01_left_18_GRF.mot'; % Must be mot or sto file
IKDataFile = 'Patient4_TMGait_0pt8_01_left_18_ik.mot'; % Must be mot or sto file
IDDataFile = 'Patient4_TMGait_0pt8_01_left_18_id.sto'; % Must be mot or sto file

% Load data from files
[MkrDataFast, MkrNamesFasat, MkrRateFast, MkrTimeFast] = readTRCFile(MkrDataFile);
[GRFDataFast, ~, GRFTimeFast] = readMotFile(GRFDataFile);
[IKDataFast, ~, IKTimeFast] = readMotFile(IKDataFile);
[IDDataFast, ~, IDTimeFast] = readMotFile(IDDataFile);

%% Format data

% Remove extra time frames from the data
MkrData = MkrData(21:121,:)/1000;
MkrTime = MkrTime(21:121,:);
GRFData = GRFData(21:121,:);
GRFTime = GRFTime(21:121,:);
IKData = IKData(21:121,:);
IKTime = IKTime(21:121,:);
IDData = IDData(21:121,7:end);
IDTime = IDTime(21:121,:);

MkrDataFast = MkrDataFast(21:121,:)/1000;
MkrTimeFast = MkrTimeFast(21:121,:);
GRFDataFast = GRFDataFast(21:121,:);
GRFTimeFast = GRFTimeFast(21:121,:);
IKDataFast = IKDataFast(21:121,:);
IKTimeFast = IKTimeFast(21:121,:);
IDDataFast = IDDataFast(21:121,7:end);
IDTimeFast = IDTimeFast(21:121,:);

% Convert IKData to radians
IKData(:,[1:3 7:end]) = IKData(:,[1:3 7:end])*pi/180;
IKDataFast(:,[1:3 7:end]) = IKDataFast(:,[1:3 7:end])*pi/180;

Torques = IDData;
q = IKData;

if GPOPSResult
    load GPOPSResults_0pt5_left_01_02.mat
end

MkrNamesOSIM = {'R.Shoulder' 'L.Shoulder' 'Neck' 'Chest' 'R.Elbow' 'R.Wrist' 'L.Elbow' 'L.Wrist' 'R.ASIS' 'L.ASIS'...
    'Sacral' 'R.Thigh.Superior' 'R.Thigh.Inferior' 'R.Thigh.Lateral' 'R.Shank.Superior' 'R.Shank.Inferior' 'R.Shank.Lateral'...
    'R.Heel' 'R.Midfoot.Superior' 'R.Midfoot.Lateral' 'R.Toe' 'L.Thigh.Superior' 'L.Thigh.Inferior' 'L.Thigh.Lateral'...
    'L.Shank.Superior' 'L.Shank.Lateral' 'L.Shank.Inferior' 'L.Heel' 'L.Midfoot.Superior' 'L.Midfoot.Lateral' 'L.Toe'};

% Determine MarkerOrder
MkrDataOrder = zeros(1,numel(MkrNamesOSIM));
for i = 1:numel(MkrNamesOSIM)
    MkrDataOrder(i) = find(ismember(MkrNames, MkrNamesOSIM{i}));
end
% Create new reordered marker matrix
MkrDataR = zeros(101,length(MkrDataOrder)*3);
for i = 1:length(MkrDataOrder)
    j = MkrDataOrder(i);
    MkrDataR(:,(i-1)*3+1:3*i) = MkrData(:,(j-1)*3+1:3*j);
end

MkrNames = MkrNames(MkrDataOrder);
MkrData = MkrDataR;

% Create left and right GRF matrices
GRFDataL = GRFData(:,1:6);
ECL = GRFData(1,13:15);
GRFDataR = GRFData(:,7:12);
ECR = GRFData(1,16:18);

GRFDataL(abs(GRFDataL(:,2))<10,:)=0;
GRFDataR(abs(GRFDataR(:,2))<10,:)=0;

time = MkrTime-MkrTime(1);

for i = 1:size(q,2)
   [~, qp(:,i)] = SplineSmooth(time,q(:,i),0.99999); 
   [~, qpp(:,i)] = SplineSmooth(time,qp(:,i),0.99999);
   [~, qppp(:,i)] = SplineSmooth(time,qpp(:,i),0.99999);
end

for i = 1:size(MkrData,2)
    [~, MkrDatap(:,i)] = SplineSmooth(time,MkrData(:,i),1);
end

%% Define moment calculation point (midfoot marker projected onto the ground)

% Determine heel and toe marker coordinates
FootMkrNames = {'R.Midfoot.Superior' 'L.Midfoot.Superior'};

FootMkrs = zeros(1,numel(FootMkrNames));
for i = 1:numel(FootMkrNames)
    FootMkrs(i) = find(ismember(MkrNames, FootMkrNames{i}));
end

% Define points to calculate joint moments about
MCR = MkrDataR(:,3*(FootMkrs(1)-1)+1:3*FootMkrs(1));%(MkrDataR(:,3*(FootMkrs(1)-1)+1:3*FootMkrs(1))+MkrDataR(:,3*(FootMkrs(2)-1)+1:3*FootMkrs(2)))/2;
MCL = MkrDataR(:,3*(FootMkrs(2)-1)+1:3*FootMkrs(2));%(MkrDataR(:,3*(FootMkrs(3)-1)+1:3*FootMkrs(3))+MkrDataR(:,3*(FootMkrs(4)-1)+1:3*FootMkrs(4)))/2;

MCR(:,2) = 0;
MCL(:,2) = 0;

% Transfer GRFs to moment centers
GRFDataR(:,4:6) = cross(ones(101,1)*ECR-MCR,GRFDataR(:,1:3),2)+GRFDataR(:,4:6);
GRFDataL(:,4:6) = cross(ones(101,1)*ECL-MCL,GRFDataL(:,1:3),2)+GRFDataL(:,4:6);

%% Load spring locations

load spring_pos_for_polynomial_r.mat

SpringsPoly_r = [SpringsH; SpringsT];

SpringsPoly_r(:,1) = SpringsPoly_r(:,1)-min(SpringsPoly_r(:,1))-(max(SpringsPoly_r(:,1))-min(SpringsPoly_r(:,1)))/2;
SpringsPoly_r(:,3) = SpringsPoly_r(:,3)-min(SpringsPoly_r(:,3))-(max(SpringsPoly_r(:,3))-min(SpringsPoly_r(:,3)))/2;

load spring_pos_for_polynomial_l.mat

SpringsPoly_l = [SpringsH; SpringsT];

SpringsPoly_l(:,1) = SpringsPoly_l(:,1)-min(SpringsPoly_l(:,1))-(max(SpringsPoly_l(:,1))-min(SpringsPoly_l(:,1)))/2;
SpringsPoly_l(:,3) = SpringsPoly_l(:,3)-min(SpringsPoly_l(:,3))-(max(SpringsPoly_l(:,3))-min(SpringsPoly_l(:,3)))/2;

% SpringsPoly_l(:,1) = SpringsPoly_r(:,1);
% SpringsPoly_l(:,3) = -SpringsPoly_r(:,3);

load SpringLocsR.mat

SpringsHR = SpringsH;
SpringsTR = SpringsT;

nSprings = length(SpringsHR)+length(SpringsTR);

% Define body index where springs are located
SpringBodyHR = 28;
SpringBodyTR = 39;

load SpringLocsL.mat

SpringsHL = SpringsH;
SpringsTL = SpringsT;

numSpringsBody = [length(SpringsHR) length(SpringsTR)];

nSprings = length(SpringsHR)+length(SpringsTR)+nSprings;

% Define body index where springs are located
SpringBodyHL = 27;
SpringBodyTL = 38;

numSpringsBody = [numSpringsBody length(SpringsHL) length(SpringsTL)];
numSprings = sum(numSpringsBody);

auxdata.latchvel = 5e-2;

%% Load spring constants

auxdata.SpringMat = [zeros(4,3); 0.0865527 0.03508 -0.0160007;...
    0.0805573 0.0391741 0.00583226; SpringsHR; SpringsTR; SpringsHL; SpringsTL];
auxdata.SpringBodyMat = [8;9;16;17;28;27;...
    SpringBodyHR*ones(length(SpringsHR),1);...
    SpringBodyTR*ones(length(SpringsTR),1);...
    SpringBodyHL*ones(length(SpringsHL),1);...
    SpringBodyTL*ones(length(SpringsTL),1)];
auxdata.Cval = 1e-2*ones(1,numSprings);
auxdata.mu = 0.5;

if SpringConstants == 1
    load SpringConstants.mat
    auxdata.Kval = Kval; % spring constants
    auxdata.Yval = Yval; % Spring height off of foot
end


%% Load EtM parameters that define activation to moment characteristics of muscles

load EtM_params_Patient4_v22f_optModel5_0pt5_l_01_02.mat
% load AlteredExcitations.mat
whichLoadsEMG = [1 2 4 5 6 8 9 11 12 13];
whichLoadsActuators = 1:25;
whichLoadsActuators(:,whichLoadsEMG) = [];
% IDData(:,whichLoadsEMG) = [];
nMusc = size(EMG,2);

% switch order of EMG and activation data to be consistent with OpenSim
% model
a = [a(:,36:70) a(:,1:35)];
EMG = [EMG(:,36:70) EMG(:,1:35)];

% Switch order of EMGs since left leg is listed first in file and right leg
% is calculated first in optimization.
Excitations = [Excitations(1:101,36:70) Excitations(1:101,1:35)];
EMGScaleOpt = [EMGScaleOpt(:,36:70) EMGScaleOpt(:,1:35)];
tdelay = [tdelay(2) tdelay(1)];

% Reference indicating which muscles cross which joints
muscRef = zeros(1,nMusc/2);
muscRef(1:17) = 1;
muscRef([18:20 22]) = 2;
muscRef([21 23:25]) = 3;
muscRef([26:27]) = 4;
muscRef([28:35]) = 5;
muscRef = [muscRef muscRef+5];

%% Scale variables using canonical scaling

tscale = 1;
tscale2 = tscale^2;
tscale3 = tscale^3;
Lscale = 1; %Note, cannot scale length becuase it requires changing osim file
Mscale = 1; %Note, cannot scale mass becuase it requires changing osim file

% Scale Time
time = time*tscale;

% Scale velocities
qp = qp/tscale;
% Scale accelerations
qpp = qpp/tscale2;
% Scale forces and torques
GRFDataR = GRFDataR/tscale2;
GRFDataL = GRFDataL/tscale2;
IDData = IDData/tscale2;
% Scale spring constants
% Kvals = Kvals/tscale2;
auxdata.latchvel = auxdata.latchvel/tscale;

%% Create curvature minimizing data splines for use in optimal control problem

[spGRFR, GRFRNewData] = spaps(time,GRFDataR',1/tscale2);
[spGRFL, GRFLNewData] = spaps(time,GRFDataL',1/tscale2);
[spMkrData, MkrNewData] = spaps(time,MkrData',.00001);
[spIK, IKNewData] = spaps(time,IKData',.01*pi/180);
[spID, IDNewData] = spaps(time,IDData',0.1/tscale2);
[spTorques, TorquesNew] = spaps(time,Torques',0.00001);
[spq, qNewData] = spaps(time,q',0.0000001);
[spqp, qpNewData] = spaps(time,qp',1/tscale);
[spEMG, EMGNewData] = spaps(time,Excitations',.001);
[spActs, ActsNewData] = spaps(time,a',.001);
[spMCR,MCRnew] = spaps(time,MCR',0.0000001);
[spMCL,MCLnew] = spaps(time,MCL',0.0000001);

%% Set output auxiliary data matrix

auxdata.Time = time; %time vector
auxdata.spGRFR = spGRFR; % right ground reactions spline
auxdata.spGRFL = spGRFL; % left ground reactions spline
auxdata.spMkrData = spMkrData; % marker data spline
auxdata.spIK = spIK; % IK data spline
auxdata.spID = spID; % ID data spline
auxdata.spTorques = spTorques;
auxdata.spq = spq; % joint angles spline
auxdata.spqp = spqp; % joint velocities spline
auxdata.SpringsPoly_r = SpringsPoly_r;
auxdata.SpringsPoly_l = SpringsPoly_l;
auxdata.SpringsHR = SpringsHR; % spring locations right heel
auxdata.SpringsTR = SpringsTR; % spring locations right toes
auxdata.SpringsHL = SpringsHL; % spring locations left heel
auxdata.SpringsTL = SpringsTL; % spring locations left toes
auxdata.SpringBodyHR = SpringBodyHR; % index number of right heel body in opensim
auxdata.SpringBodyTR = SpringBodyTR; % index number of right toes body in opensim
auxdata.SpringBodyHL = SpringBodyHL; % index number of left heel body in opensim
auxdata.SpringBodyTL = SpringBodyTL; % index number of left toes body in opensim
auxdata.GRFDataR = GRFDataR; % GRF data right
auxdata.GRFDataL = GRFDataL; % GRF data left
auxdata.MkrData = MkrData; % marker data
auxdata.IKData = IKData; % IK data
auxdata.IDData = IDData; % ID data
auxdata.Torques = Torques;
auxdata.q = q; % joint angles
auxdata.qp = qp; % joint velocities
auxdata.qpp = qpp; % joint accelerations
auxdata.qppp = qppp; % joint jerk
auxdata.MkrNames = MkrNames; % marker names
auxdata.numSprings = nSprings; % number of springs
auxdata.numSpringsBody = numSpringsBody; % number of springs on each foot body
auxdata.SpringBodyHR = SpringBodyHR; % index number of right heel body in opensim
auxdata.SpringBodyTR = SpringBodyTR; % index number of right toes body in opensim
auxdata.SpringBodyHL = SpringBodyHL; % index number of left heel body in opensim
auxdata.SpringBodyTL = SpringBodyTL; % index number of left toes body in opensim
auxdata.ECR = ECR; % right force plate electrical center
auxdata.ECL = ECL; % left force plate electrical center
auxdata.MCR = MCR; % right moment calculation point
auxdata.MCL = MCL; % left moment calculation point
auxdata.spMCR = spMCR; % right moment calculation point spline
auxdata.spMCL = spMCL; % left moment calculation point spline
auxdata.tscale = tscale; % time scale factor
auxdata.lmo = lmo; % muscle optimal fiber length
auxdata.lts = lts; % muscle tendon slack length
auxdata.Fmax = Fmax; % muscle maximum force
auxdata.pennAngle = pennAngle; % muscle pennation angle
auxdata.EtMcoefs = coefs; % coefficients defining muscle lengths and moment arms as functions of kinematics
auxdata.tact = tact; % muscle activation time constants
auxdata.tdelay = tdelay; % muscle time delay
auxdata.Anonlin = Anonlin; % muscle nonlinearity parameter
auxdata.EMG = EMG; % processed and normalized EMG data
auxdata.Acts = a;
auxdata.Excitations = Excitations;
auxdata.spEMG = spEMG; % EMG data spline
auxdata.spActs = spActs;
auxdata.whichLoadsEMG = whichLoadsEMG; % define which joints are actuated by EMG
auxdata.whichLoadsActuators = whichLoadsActuators; % which joints are actuated by torques
auxdata.muscRef = muscRef; % Defines functions of muscles (which joints each muscle crosses)
auxdata.EMGScaleOpt = EMGScaleOpt;

auxdata.qFast = IKDataFast;
auxdata.TorquesFast = IDDataFast;

function [Data, ColumnLabels, Time] = readMotFile(inFile)

%open file
fid = fopen(inFile);

%skip initial lines
line = 'asdf';
while ~strcmpi(line,'ENDHEADER')
    line = fgetl(fid);
    if length(line)>8
        line = line(1:9);
    end
end

% Find the labels on the columns
ColumnLabels = fgetl(fid);
ColumnLabels = textscan(ColumnLabels, '%t');
ColumnLabels = ColumnLabels{1}';

% Read numeric data
Data = fscanf(fid, '%f');

% Find length of data
nColumns = length(ColumnLabels);
nRows = length(Data)/nColumns;

% Reshape data into matrix
Data = reshape(Data, nColumns, nRows)';

% Extract time vector
Time = Data(:,1);
Data(:,1) = [];

% Remove time label from column labels
ColumnLabels(1) = [];

fclose(fid);

function [MarkerData, MarkerNames, SampleRate, Time] = readTRCFile(inFile)

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