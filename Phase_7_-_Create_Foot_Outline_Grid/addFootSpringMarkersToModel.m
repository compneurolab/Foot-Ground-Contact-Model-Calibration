function addFootSpringMarkersToModel(varargin)
% Adds Markers to model based on the spring locations in order to visualize
% spring placement on the foot in the OpenSim GUI.


%import OpenSim Matlab API functions
import org.opensim.modeling.*
import java.lang.*

inModel = 'Patient4_optModel5_GPOPS.osim'
outModel='ModelWithSprings.osim';

load SpringLocsL.mat
SpringsHL = SpringsH;
SpringsTL = SpringsT;
load SpringLocsR.mat
SpringsHR = SpringsH;
SpringsTR = SpringsT;

%setup OpenSim model
model=Model(inModel);

Body = model.getBodySet.get(27);
for i = 1:size(SpringsHL,1)
    Name = String(['Left Heel Marker' num2str(i)]);
    marker=model.getMarkerSet.addMarker(Name, SpringsHL(i,:), Body);
end

Body = model.getBodySet.get(38);
for i = 1:size(SpringsTL,1)
    Name = String(['Left Toe Marker' num2str(i)]);
    marker=model.getMarkerSet.addMarker(Name, SpringsTL(i,:), Body);
end

Body = model.getBodySet.get(28);
for i = 1:size(SpringsHR,1)
    Name = String(['Right Heel Marker' num2str(i)]);
    marker=model.getMarkerSet.addMarker(Name, SpringsHR(i,:), Body);
end
%keyboard
Body = model.getBodySet.get(39);
for i = 1:size(SpringsTR,1)
    Name = String(['Right Toe Marker' num2str(i)]);
    marker=model.getMarkerSet.addMarker(Name, SpringsTR(i,:), Body);
end

%output, i.e., plate, .osim model
model.setName(outModel);%set name
model.print(outModel);%print