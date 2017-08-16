function phaseout = osimOptimalTrackingContinuous(input)

% Define persistent variables so that splines are not conducted on each
% continuous function call
persistent MkrData GRFR GRFL ActsExp IDData MCR MCL q

% rescale everything back to original
nframes = length(input.phase.time);
onescol = ones(nframes,1);

% inputorig = input;

input.phase.time = (input.phase.time)*(input.auxdata.tmax-input.auxdata.tmin)+(input.auxdata.tmax+input.auxdata.tmin)/2;
input.phase.state = (input.phase.state).*(onescol*(input.auxdata.statemax-input.auxdata.statemin))+onescol*(input.auxdata.statemax+input.auxdata.statemin)/2;
input.phase.control = (input.phase.control).*(onescol*(input.auxdata.controlmax-input.auxdata.controlmin))+onescol*(input.auxdata.controlmax+input.auxdata.controlmin)/2;
input.phase.parameter = (input.phase.parameter(1,:)).*((input.auxdata.parammax-input.auxdata.parammin))+(input.auxdata.parammax+input.auxdata.parammin)/2;

%% Load parameters from the auxilary data
tscale = input.auxdata.tscale;

mu_s = 0; % Coefficient of static friction for foot-ground contact
mu_d = 1;%input.phase.parameter(1,9);%input.auxdata.mu;%parameters(1,1); % Coefficient of dynamic friction
mu_v = 0; % Coefficient of viscous friction
KvalFits = input.phase.parameter(1,1:6);%input.auxdata.Kval; % Spring stiffnesses for foot ground contact
Cval = input.auxdata.Cval; % Spring damping
Yval = input.phase.parameter(1,7:8);
SpringMat = input.auxdata.SpringMat;
SpringBodyMat = input.auxdata.SpringBodyMat;
ECR = input.auxdata.ECR;
ECL = input.auxdata.ECL;
numSpringsBody = input.auxdata.numSpringsBody;

% SpringMat(7:end,2) = SpringMat(7:end,2)+Yval(1);
SpringMat(7:6+input.auxdata.numSprings/2,2) = SpringMat(7:6+input.auxdata.numSprings/2,2)+Yval(1);
SpringMat(7+input.auxdata.numSprings/2:end,2) = SpringMat(7+input.auxdata.numSprings/2:end,2)+Yval(2);

% Transition velocity between static and dynamic friction
latchvel = input.auxdata.latchvel;

muscParams.lmo = input.auxdata.lmo;
muscParams.lts = input.auxdata.lts;
muscParams.Fmax = input.auxdata.Fmax;
muscParams.pennAngle = input.auxdata.pennAngle;
muscParams.EtMcoefs = input.auxdata.EtMcoefs;
muscParams.nframes = nframes;
muscParams.muscRef = input.auxdata.muscRef;

%% Calculate Spring Constant Values
Kval = zeros(1,input.auxdata.numSprings);

SpringsPoly_l = input.auxdata.SpringsPoly_l;
SpringsPoly_r = input.auxdata.SpringsPoly_r;

SpringFitsMatrix_r = [ones(length(SpringsPoly_r),1) SpringsPoly_r(:,1) ...
    SpringsPoly_r(:,3) SpringsPoly_r(:,1).*SpringsPoly_r(:,3) ...
    SpringsPoly_r(:,1).^2 SpringsPoly_r(:,3).^2];

SpringFitsMatrix_l = [ones(length(SpringsPoly_l),1) SpringsPoly_l(:,1) ...
    -SpringsPoly_l(:,3) -SpringsPoly_l(:,1).*SpringsPoly_l(:,3) ...
    SpringsPoly_l(:,1).^2 SpringsPoly_l(:,3).^2];

Kval(1:numSpringsBody(1)+numSpringsBody(2)) = SpringFitsMatrix_r*KvalFits(1:6)';
Kval(numSpringsBody(1)+numSpringsBody(2)+1:end) = SpringFitsMatrix_l*KvalFits(1:6)';

%% Initialize Variables

t = input.phase(1).time;

if ~exist('MkrData','var')
    q = fnval(input.auxdata.spq,t)';
    MkrData = fnval(input.auxdata.spMkrData,t)';
    GRFR = fnval(input.auxdata.spGRFR,t)';
    GRFL = fnval(input.auxdata.spGRFL,t)';
    IDData = fnval(input.auxdata.spID,t)';
    MCR = fnval(input.auxdata.spMCR,t)';
    MCL = fnval(input.auxdata.spMCL,t)';
elseif size(MkrData,1) ~= nframes
    q = fnval(input.auxdata.spq,t)';
    MkrData = fnval(input.auxdata.spMkrData,t)';
    GRFR = fnval(input.auxdata.spGRFR,t)';
    GRFL = fnval(input.auxdata.spGRFL,t)';
    IDData = fnval(input.auxdata.spID,t)';
    MCR = fnval(input.auxdata.spMCR,t)';
    MCL = fnval(input.auxdata.spMCL,t)';
end

%% Kinematics

% Set states
x = input.phase(1).state(:,1:31);
xp = input.phase(1).state(:,32:62);

% Calculate spring positions and velocities
[SpringPos, SpringVels] = opensimPointKin_v3(t,[x xp],SpringMat',SpringBodyMat);

% Store location of heel and toes body for tranferring GRFs
BodyPositions = SpringPos(:,1:12);
% MCR = SpringPos(:,13:15);
% MCL = SpringPos(:,16:18);
SpringPos(:,1:18) = [];
SpringVels(:,1:18) = [];

MCR(:,2) = 0;
MCL(:,2) = 0;

%% Calculate Ground Reaction Forces

% Calculate spring forces and moments about electrical centers
[FootGRFsRMC, FootGRFsLMC] = calcGroundReactions_v3(SpringPos, SpringVels, ...
    numSpringsBody, mu_s, mu_d, mu_v, Kval, Cval, MCR, MCL, latchvel, tscale,input.auxdata.BeltSpeed);

% Transfer moments to heel and toe origin for input into full body model
% This is necessary because the contact model ground reactions are
% calculated about some point in the room (force plate electrical centers),  
% and the model applies the loads at the origin of the foot and toes body
FootGRFsR = FootGRFsRMC;
FootGRFsL = FootGRFsLMC;

% Transfer to heel
FootGRFsR(:,7:9) = cross(MCR-BodyPositions(:,1:3),FootGRFsRMC(:,1:3),2)+FootGRFsRMC(:,7:9);
FootGRFsL(:,7:9) = cross(MCL-BodyPositions(:,7:9),FootGRFsLMC(:,1:3),2)+FootGRFsLMC(:,7:9);

% Transfer to toe
FootGRFsR(:,10:12) = cross(MCR-BodyPositions(:,4:6),FootGRFsRMC(:,4:6),2)+FootGRFsRMC(:,10:12);
FootGRFsL(:,10:12) = cross(MCL-BodyPositions(:,10:12),FootGRFsLMC(:,4:6),2)+FootGRFsLMC(:,10:12);

%% Dynamics

% Extract for input into full body dynamic model
Torques = zeros(nframes,25);%Set to zero because they will be determined by inverse dynamics

% Formulate inputs into opensim (states and controls w/ GRFs)
States = [x xp];
AppliedLoads = [Torques FootGRFsR FootGRFsL];
% xpp = input.phase(1).control(:,1:29);

States(isnan(States)) = 0;
States(isinf(States)) = 0;

% Set accelerations
xpp = input.phase.state(:,63:end);
% Set jerk
xppp = input.phase.control;

% Calculate state derivatives and marker locations
[ResidualLoads, mkrPos] = Control_OpenSimMPI_muscle_ID([x xp], xpp, AppliedLoads,t,[1 0 0]);

% Set rescaled derivative of state vector
phaseout(1).dynamics = (input.auxdata.tmax-input.auxdata.tmin)*([xp xpp xppp])./(onescol*(input.auxdata.statemax-input.auxdata.statemin));

phaseout(1).path = ([ResidualLoads(:,1:6)]-...
    ones(nframes,1)*(input.auxdata.pathmax+input.auxdata.pathmin)/2)...
    ./(ones(nframes,1)*(input.auxdata.pathmax-input.auxdata.pathmin));

Torques = ResidualLoads(:,7:end);

%% Cost Function

% Output GRFs split between segments used for creating plots after
% simulation finishes
phaseout.GRFsSplit = [FootGRFsRMC FootGRFsLMC];

% Calculate net GRFs at the electrical center for use in objective function
FootGroundReactionsRMC = FootGRFsRMC(:,[1:3 7:9])+FootGRFsRMC(:,[4:6 10:12]);
FootGroundReactionsLMC = FootGRFsLMC(:,[1:3 7:9])+FootGRFsLMC(:,[4:6 10:12]);

mkrError = (mkrPos-MkrData);
GRFError = [GRFR-FootGroundReactionsRMC GRFL-FootGroundReactionsLMC];
% activationError = activations-ActsExp;
TorqueErrors = Torques-IDData;
qerror = x-q;

% Set integrand, we penalize differences between model and experimental
% markers and differences between model and experimental ground reactions
% calculated at the electrical centers of the force plates, additionally we
% minimize the rate of change of the applied joint loads
phaseout(1).integrand = [mkrError...
    GRFError(:,[1:3 7:9]) ...
    GRFError(:,[4:6 10:12]) ...
    TorqueErrors(:,[1:6,8:13,15:end])...
    xppp].^2;

% Scale integrand
phaseout.integrand = (phaseout.integrand)./(onescol*(input.auxdata.integralmax-input.auxdata.integralmin));

% phaseout.EtMmoments = EtMmoments;
% Set more outputs for post-simulation analyses
phaseout(1).markers = mkrPos;
phaseout.Torques = Torques;
phaseout.Kval = Kval;
phaseout.Yval = Yval;
