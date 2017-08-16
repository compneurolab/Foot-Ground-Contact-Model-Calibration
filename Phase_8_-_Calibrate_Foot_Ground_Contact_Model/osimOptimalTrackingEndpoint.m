%-------------------------------------------%
% BEGIN: function brachistochroneEndpoint.m %
%-------------------------------------------%
function output = osimOptimalTrackingEndpoint(input)

% % Variables at Start and Terminus of Phase 1
t0{1} = input.phase(1).initialtime;
tf{1} = input.phase(1).finaltime;
x0{1} = input.phase(1).initialstate;
xf{1} = input.phase(1).finalstate;

integral = input.phase(1).integral;

% Objective
q = sum(integral([1:93]))/10000+sum(integral([94:99]))/10000+sum(integral([100:105]))/10000 ...
    +sum(integral([106:end]))/10000;%+SpringPen;%+outputContinuous(1).integrand*(tf{1}-t0{1})/100;%+(10^4)*(tf{1}-1);%+input.phase(2).integral;
output.objective = q;
%-------------------------------------------%
% END: function brachistochroneEndpoint.m   %
%-------------------------------------------%

