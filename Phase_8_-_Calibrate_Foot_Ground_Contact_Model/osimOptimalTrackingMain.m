%---------------------------------------------------%
% Foot Motion Problem:                              %
%---------------------------------------------------%

clear all; close all;

Initial_Guess = 1;

% addpath('C:\Files\GPOPS-II New\nlp\snopt')

% Load model parameters determined by past calibration runs
cdorig = cd;
cd('../SetupFiles')
auxdata = SetupOpenSimOptimalTracking(0);
cd(cdorig)

% Load initial and final model pose and time
t0 = auxdata.Time(1);
tf = auxdata.Time(end);
Time = auxdata.Time;
q0 = auxdata.q(1,:);
qp0 = auxdata.qp(1,:);
qf = auxdata.q(end,:);
qpf = auxdata.qp(end,:);

% Set initial guesses for Loads (control) and joint angle trajectories
% (states)
LoadsGuess = auxdata.IDData;
qGuess = auxdata.q;
qpGuess = auxdata.qp;
qppGuess = auxdata.qpp;
qpppGuess = auxdata.qppp;

% Determine the range of the states
qrange = max(qGuess)-min(qGuess);
qprange = max(qpGuess)-min(qpGuess);
qpprange = max(qppGuess)-min(qppGuess);
qppprange = max(qpppGuess)-min(qpppGuess);

auxdata.activeDOFs = logical(1-(qrange == 0));
auxdata.activeDOFs(1:6) = [];

% Give locked DOFs a range of motion
qppprange(qrange == 0) = 40000*pi/180;
qpprange(qrange == 0) = 4000*pi/180;
qprange(qrange == 0) = 400*pi/180;
qrange(qrange == 0) = 60*pi/180;

% Set bounds on the states based on the previously calcualted range
qmin = min(qGuess)-2*qrange;  qmax  = max(qGuess)+2*qrange;
qpmin = min(qpGuess)-1.5; qpmax  = max(qpGuess)+1.5*qprange;
qppmin = min(qppGuess)-1*qpprange; qppmax  = max(qppGuess)+1*qpprange;
qpppmin = min(qpppGuess)-1*qppprange; qpppmax  = max(qpppGuess)+1*qppprange;

% Determing the range of the controls
LoadsRange = max(auxdata.IDData)-min(auxdata.IDData);

% Give bigger range to DOFs with small loads (for instance the toes)
LoadsRange(LoadsRange<1/(auxdata.tscale^2)) = 25/(auxdata.tscale^2);

% Set bounds on the controls
umin = min(auxdata.IDData)-2*LoadsRange; umax = max(auxdata.IDData)+2*LoadsRange;

% Set weights for marker errors (the model tracks markers from experimental
% data
auxdata.ErrorWeights = ones(1,numel(auxdata.MkrNames));
auxdata.ErrorWeights([18 21 28 31]) = 10;
auxdata.ErrorWeights([3:8]) = .1;

auxdata.ErrorWeights = repmat(auxdata.ErrorWeights,3,1);
auxdata.ErrorWeights = auxdata.ErrorWeights(:)';

% Initilize osim models in mex files
opensimPointKin_v3('Patient4_optModel5_GPOPS.osim');
% Control_OpenSimMPIv3('Patient4_optModel2_GPOPS.osim');

MkrErrorsAllowed = 2e-2*ones(1,numel(auxdata.MkrNames));
MkrErrorsAllowed([18:21 28:31]) = 1.5e-2;
MkrErrorsAllowed([3:8]) = 3e-2;
MkrErrorsAllowed = repmat(MkrErrorsAllowed, 3,1);% Duplicate for each dimension
MkrErrorsAllowed = MkrErrorsAllowed(:)';

auxdata.numSprings = sum(auxdata.numSpringsBody);

if Initial_Guess==1
% load previous solution for use as an initial guess
load solution_osimOptTracking_calibrateSprings29.mat
% load output_osimOptTracking_jerk.mat
end

auxdata.BeltSpeed = 0.5;

% solution = output.result.interpsolution;

% Setup scale factors
auxdata.tmin = t0;
auxdata.tmax = tf;
auxdata.statemin = [qmin,qpmin,qppmin];
auxdata.statemax = [qmax,qpmax,qppmax];
auxdata.controlmin = [qpppmin];
auxdata.controlmax = [qpppmax];
auxdata.integralmin = 0*(ones(1,93+6+6+23+31));
auxdata.integralmax = auxdata.tmax*([MkrErrorsAllowed ...
    1e0*[1e1 2e1 1e1]...
    1e0*[1e1 2e1 1e1]... 
    1e1*ones(1,3)...
    1e1*ones(1,3)...
    1*LoadsRange([1:6,8:13,15:end])...
    1*qppprange]).^2;
onescol = ones(length(auxdata.Time),1);
auxdata.pathmin = 1e0*[-1e-1*ones(1,3) -1e0*ones(1,3)];
auxdata.pathmax = 1e0*[1e-1*ones(1,3) 1e0*ones(1,3)];
%Specify bounds of spring constants
auxdata.parammin = [0 -2e3*ones(1,5) -0.01*ones(1,2)];
auxdata.parammax = [1e4 2e3*ones(1,5) 0.01*ones(1,2)];

if 1 %Run Opt    
%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%

nstates = length(auxdata.statemin);
ncontrols = length(auxdata.controlmin);
% nparams = length(auxdata.paramsmin);

iphase = 1;
bounds.phase(iphase).initialtime.lower = -.5;
bounds.phase(iphase).initialtime.upper = -.5;
bounds.phase(iphase).finaltime.lower = .5;
bounds.phase(iphase).finaltime.upper = .5;
% states are the joint angles, velocities, accelerations and the joint
% torques
bounds.phase(iphase).initialstate.lower = -.5*ones(1,nstates);%,vmin];
bounds.phase(iphase).initialstate.upper = .5*ones(1,nstates);%,vmax];
bounds.phase(iphase).finalstate.lower = -.5*ones(1,nstates);%,vmin];%[qf,qpf];
bounds.phase(iphase).finalstate.upper = .5*ones(1,nstates);%,vmax];%[qf,qpf];
bounds.phase(iphase).state.lower = -.5*ones(1,nstates);%,vmin];
bounds.phase(iphase).state.upper = .5*ones(1,nstates);%,vmax];
% Dynamics equations are forced to be satisfied by path constraints
bounds.phase(iphase).path.lower = -0.5*ones(1,6); % Bounds on pelvis residuals
bounds.phase(iphase).path.upper = 0.5*ones(1,6);
% Control with joint angle jerk and the time derivative of the joint loads
bounds.phase(iphase).control.lower = -.5*ones(1,ncontrols);
bounds.phase(iphase).control.upper = .5*ones(1,ncontrols);
bounds.phase(iphase).integral.lower = 0*(ones(size(auxdata.integralmin)));
bounds.phase(iphase).integral.upper = 1*(ones(size(auxdata.integralmax)));
bounds.parameter.lower = -2*ones(1,8);
bounds.parameter.upper = 2*ones(1,8);

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
% solution.phase.control(:,30:99) = solution.phase.control(:,30:99)/10;
% Set initial guess either from previous optimal control solution or from
% the auxilary data which was calculated from a different model
iphase = 1;
if Initial_Guess ==1;
    guess = solution;
else
    guess.phase(iphase).time    = ([auxdata.Time]-(auxdata.tmax+auxdata.tmin)/2)./(auxdata.tmax-auxdata.tmin);
    guess.phase(iphase).state   = ([qGuess qpGuess qppGuess]-onescol*(auxdata.statemax+auxdata.statemin)/2)./(onescol*(auxdata.statemax-auxdata.statemin));% auxdata.Loads];% auxdata.Loadsp]
    guess.phase(iphase).control = ([qpppGuess]-onescol*(auxdata.controlmax+auxdata.controlmin)/2)./(onescol*(auxdata.controlmax-auxdata.controlmin)); %zeros(size(auxdata.LoadsExp(:,:)));
%     guess.parameter = [0 0]; % Must be active for model calibration
    guess.phase(iphase).integral = (1e1-(auxdata.integralmax+auxdata.integralmin)/2)./(auxdata.integralmax-auxdata.integralmin);
    guess.parameter = ([3e3 zeros(1,5) 0 0]-(auxdata.parammax+auxdata.parammin)/2)./(auxdata.parammax-auxdata.parammin);
    %     guess.parameter = (parameterGuess-(auxdata.paramsmax+auxdata.paramsmin)/2)./(auxdata.paramsmax-auxdata.paramsmin);
end

%-------------------------------------------------------------------------%
%------------- Set up Event Constraints That Link Phases -----------------%
%-------------------------------------------------------------------------%
% bounds.eventgroup(1).lower = ;
% bounds.eventgroup(1).upper = ;

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%
%-------------------------------------------------------------------------%
setup.name = 'openSim Full Body Motion';
setup.functions.continuous = @osimOptimalTrackingContinuous;
auxdata.ContinuousFunc = setup.functions.continuous;
setup.functions.endpoint = @osimOptimalTrackingEndpoint;
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.nlp.solver ='snopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.tolerance = 1e-5;
setup.nlp.snoptoptions.tolerance = 1e-5;
setup.nlp.ipoptoptions.maxiterations = 100;
setup.derivatives.stepsize1 = 1e-9;
setup.derivatives.supplier = 'sparseFD';
setup.derivatives.derivativelevel = 'first';
setup.derivatives.dependencies = 'sparse';
% setup.scales.method='automatic-bounds';
mesh.method = 'hp-PattersonRao';
mesh.tolerance       = 1e-4;
mesh.maxiterations    = 0;
mesh.colpointsmin    = 2;
mesh.colpointsmax    = 10;
N                    = 5;
mesh.phase.colpoints = 10*ones(1,N);
mesh.phase.fraction  = ones(1,N)/N;
setup.method = 'RPM-integration';
setup.mesh = mesh;

%-------------------------------------------------------------------------%
%------------------------- Solve Problem Using GPOP2 ---------------------%
%-------------------------------------------------------------------------%
tic
output = gpops2(setup);
toc
solution = output.result.solution;
save solution_osimOptTracking_calibrateSprings29.mat solution
save output_osimOptTracking_calibrateSprings29.mat output

end

%--------------------------------------------------------------------------%
%------------------------------- Plot Solution ----------------------------%
%--------------------------------------------------------------------------%

%load output_osimOptTracking_calibrateSprings4.mat

tscale = 1;%auxdata.tscale;

input = solution;
input.phase.parameter = solution.parameter;
input.auxdata = auxdata;

% outputContinuous = FullBodyOptimalControlContinuous_experimentalGRFs_Markers(input);
outputContinuous = osimOptimalTrackingContinuous(input);

onescol = ones(size(solution.phase.state(:,1)));

solution.phase.time = (solution.phase.time)*(auxdata.tmax-auxdata.tmin)+(auxdata.tmax+auxdata.tmin)/2;
solution.phase.state = (solution.phase.state).*(onescol*(auxdata.statemax-auxdata.statemin))+onescol*(auxdata.statemax+auxdata.statemin)/2;
solution.phase.control = (solution.phase.control).*(onescol*(auxdata.controlmax-auxdata.controlmin))+onescol*(auxdata.controlmax+auxdata.controlmin)/2;
solution.parameter = (solution.parameter(1,:)).*((input.auxdata.parammax-input.auxdata.parammin))+(input.auxdata.parammax+input.auxdata.parammin)/2;

mkrDataSplined = interp1(auxdata.Time,auxdata.MkrData,solution.phase.time,'pchip');

mkrDataSplined = fnval(auxdata.spMkrData,solution.phase.time)';

error = mkrDataSplined-outputContinuous.markers;

% Calculate RMS error for each marker at each time frame
for i = 1:size(mkrDataSplined,2)/3
    RMSError(:,i) = sum(error(:,(i-1)*3+1:3*i).^2,2).^(1/2);
end

GRFsSplit = outputContinuous.GRFsSplit;
figure
plot(solution.phase.time,GRFsSplit(:,1:3)+GRFsSplit(:,4:6),'b')
hold on
plot(auxdata.Time,auxdata.GRFDataR(:,1:3),'r')

figure
plot(solution.phase.time,GRFsSplit(:,13:15)+GRFsSplit(:,16:18),'b')
hold on
plot(auxdata.Time,auxdata.GRFDataL(:,1:3),'r')

figure
plot(solution.phase.time,GRFsSplit(:,7:9)+GRFsSplit(:,10:12),'b')
hold on
plot(auxdata.Time,auxdata.GRFDataR(:,4:6),'r')

figure
plot(solution.phase.time,GRFsSplit(:,19:21)+GRFsSplit(:,22:24),'b')
hold on
plot(auxdata.Time,auxdata.GRFDataL(:,4:6),'r')

figure
plot(solution.phase.time, solution.phase.control(:,1:23))
hold on
plot(auxdata.Time, auxdata.qpp)

figure
for i = 1:29
    
    subplot(6,6,i)
    plot(solution.phase.time/solution.phase.time(end), solution.phase.state(:,i),'r')
    hold on
    plot(auxdata.Time/auxdata.Time(end), auxdata.q(:,i),'b')
    
end

figure
for i = 1:29
    
    subplot(6,6,i)
    plot(solution.phase.time/solution.phase.time(end), solution.phase.control(:,i),'r')
    hold on
    plot(auxdata.Time/auxdata.Time(end), auxdata.qpp(:,i),'b')
    
end

figure
for i = 1:23
    subplot(5,5,i)
    plot(solution.phase.time/solution.phase.time(end), outputContinuous.Torques(:,i),'r')
    hold on
    plot(auxdata.Time/auxdata.Time(end), auxdata.IDData(:,i),'b')
end

figure
plot3(auxdata.SpringsPoly_r(:,1),auxdata.SpringsPoly_r(:,3),outputContinuous.Kval(1:47),'x')
figure
plot3(auxdata.SpringsPoly_l(:,1),auxdata.SpringsPoly_l(:,3),outputContinuous.Kval(48:end),'rx')

ToeTorques = outputContinuous.Torques(:,[7 14]);

ToeKinematics = solution.phase.state(:,[13 20 42 49]);

Mat1 = [ones(51,1) ToeKinematics(:,[1 3])];
Mat2 = [ones(51,1) ToeKinematics(:,[2 4])];

ApproxTorques(:,1) = Mat1*(Mat1\ToeTorques(:,1));
ApproxTorques(:,2) = Mat2*(Mat2\ToeTorques(:,2));

[spIDnew] = spaps(solution.phase.time,outputContinuous.Torques',.5);

control = fnval(spIDnew,auxdata.Time)';

NewID = control(:,[8 9 11 12 13 1 2 4 5 6]);

% Get spring constants from continuous function
Kval = outputContinuous.Kval;
Yval = outputContinuous.Yval;

% Output spring constants
save SpringConstants.mat Kval Yval

save control.mat control

save NewIDLoads.mat NewID

outputIKMot('IK.mot',solution.phase.state(:,1:29),solution.phase.time)