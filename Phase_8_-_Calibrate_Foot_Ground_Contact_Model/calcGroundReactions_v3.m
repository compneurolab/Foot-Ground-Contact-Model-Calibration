function [FspringsR, FspringsL] = calcGroundReactions(SpringPos, SpringVels,...
    numSpringsBody, mu_s, mu_d, mu_v, Kval, Cval, ECR, ECL, latchvel,tscale,BeltSpeed)
% v3 implemented Fregly's springs

nframes = size(SpringPos,1);

Rhsprings = 1:numSpringsBody(1);
Rtsprings = numSpringsBody(1)+1:numSpringsBody(1)+numSpringsBody(2);
Lhsprings = numSpringsBody(1)+numSpringsBody(2)+1:numSpringsBody(1)+numSpringsBody(2)+numSpringsBody(3);
Ltsprings = numSpringsBody(1)+numSpringsBody(2)+numSpringsBody(3)+1:numSpringsBody(1)+numSpringsBody(2)+numSpringsBody(3)+numSpringsBody(4);

% set spring positions for calculating moments about electrical centers
xposvals = SpringPos(:,1:3:end);
yposvals = SpringPos(:,2:3:end);
zposvals = SpringPos(:,3:3:end);

% set vertical position of springs (modified later by code)
% yPens = -SpringPos(:,2:3:end);

% set spring velocities
xvel = SpringVels(:,1:3:end)+BeltSpeed/tscale;
normvel = -SpringVels(:,2:3:end);
zvel = SpringVels(:,3:3:end);

% determine slip velocity of springs
slipvel = (xvel.^2+zvel.^2).^(1/2);

% determine normalized velocities of springs
slipOffset = 1e-4;
xvel = xvel./(slipvel+slipOffset);
zvel = zvel./(slipvel+slipOffset);

% % Continuous function that models in/out of contact condition
% c = 1000;
% yPens = yPens+log(exp(-c*yPens)+1)/c;
% yPens(isinf(yPens)) = 0;

klow = 1e-1/tscale^2;
h = 1e-3;
c = 5e-4;
ymax = 1e-2;

v = ones(nframes,1)*((Kval+klow)./(Kval-klow));
s = ones(nframes,1)*((Kval-klow)/2);
constant = -s.*(v.*ymax-c*log(cosh((ymax+h)/c)));
        
Fsprings = -s.*(v.*yposvals-c*log(cosh((yposvals+h)/c)))-constant;

Fsprings(isnan(Fsprings)) = min(min(Fsprings));
Fsprings(isinf(Fsprings)) = min(min(Fsprings));

% Add in slight slope for non-contact time frames in order to (maybe) help optimizer
% yPens = yPens+.00001*(yposvals);

% Calculate vGRF (Normal Force)
% Fy = (ones(nframes,1)*Kval).*yPens.*(1+(ones(nframes,1)*Cval).*(normvel));
Fy = Fsprings.*(1+(ones(nframes,1)*Cval).*(normvel));
RFyHF = sum(Fy(:,Rhsprings),2);
RFyT = sum(Fy(:,Rtsprings),2);
LFyHF = sum(Fy(:,Lhsprings),2);
LFyT = sum(Fy(:,Ltsprings),2);

% New friction model
%dynamic friction
mu = mu_d*tanh(slipvel/latchvel); 

%static friction
% b = mu_s-mu_d;
% c = latchvel;
% d = latchvel;
% spos = b*exp(-((slipvel-c).^2)/(2*d^2));
% sneg = -b*exp(-((slipvel+c).^2)/(2*d^2));
% mu = mu+spos+sneg;

% damping friction
damping = mu_v*slipvel/latchvel;
mu = mu+damping;

% Calculate horizontal forces
Fvals_individual = Fy.*mu;
% Resolve into components and sum forces for all springs
Fxvals_individual = -Fvals_individual.*xvel;
Fzvals_individual = -Fvals_individual.*zvel;

RFxHF = sum(Fxvals_individual(:,Rhsprings),2);
RFxT = sum(Fxvals_individual(:,Rtsprings),2);
RFzHF = sum(Fzvals_individual(:,Rhsprings),2);
RFzT = sum(Fzvals_individual(:,Rtsprings),2);
LFxHF = sum(Fxvals_individual(:,Lhsprings),2);
LFxT = sum(Fxvals_individual(:,Ltsprings),2);
LFzHF = sum(Fzvals_individual(:,Lhsprings),2);
LFzT = sum(Fzvals_individual(:,Ltsprings),2);

ECR = permute(ECR, [1 3 2]);
ECL = permute(ECL, [1 3 2]);

ECR = repmat(ECR,[1,sum(numSpringsBody([1 2])),1]);
ECL = repmat(ECL,[1,sum(numSpringsBody([3 4])),1]);

% Using the position vectors, calculate the moment contributions from each
% element about the electrical center
positionvec = cat(3, xposvals, yposvals, zposvals);
forcevec = cat(3, Fxvals_individual, Fy, Fzvals_individual); 
Moments = cross(positionvec-[ECR ECL],forcevec, 3);
xMoments = Moments(:,:,1);
yMoments = Moments(:,:,2);
zMoments = Moments(:,:,3);

% Sum the moments for all springs
RMxHF = sum(xMoments(:,Rhsprings),2);
RMxT = sum(xMoments(:,Rtsprings),2);
RMyHF = sum(yMoments(:,Rhsprings),2);
RMyT = sum(yMoments(:,Rtsprings),2);
RMzHF = sum(zMoments(:,Rhsprings),2);
RMzT = sum(zMoments(:,Rtsprings),2);

LMxHF = sum(xMoments(:,Lhsprings),2);
LMxT = sum(xMoments(:,Ltsprings),2);
LMyHF = sum(yMoments(:,Lhsprings),2);
LMyT = sum(yMoments(:,Ltsprings),2);
LMzHF = sum(zMoments(:,Lhsprings),2);
LMzT = sum(zMoments(:,Ltsprings),2);

% Set output
FspringsR = [RFxHF RFyHF RFzHF RFxT RFyT RFzT RMxHF RMyHF RMzHF RMxT RMyT RMzT];
FspringsL = [LFxHF LFyHF LFzHF LFxT LFyT LFzT LMxHF LMyHF LMzHF LMxT LMyT LMzT];
