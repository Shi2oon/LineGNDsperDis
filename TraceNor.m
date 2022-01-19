% to check for active dislocations system for the known slip plane
clc;clear % b and t lie in the plane if the cross product is zero
% change the number
load([pwd '\82\82_Line_Data.mat'])
tracePlane = [-1 1 2]; 
ang  = @(uy,vy) rad2deg(atan2(norm(cross(uy,vy)),dot(uy,vy)));
% if angle is 180 then they are parallel (screw dislocations)
% if angle is 90 then they are perpendicular (edge dislocations)
% for dislcoation to lie in a slip plabe both t and b should be 90 
% to the slip plane

for iV = 1:length(LN.b)
    if LN.GND_per(iV)~=0
    C(iV,1) = iV;
    C(iV,2) = ang(tracePlane,[LN.b(iV,1) LN.b(iV,2) LN.b(iV,3)]);
    C(iV,3) = ang(tracePlane,[LN.t(iV,1) LN.t(iV,2) LN.t(iV,3)]);
    C(iV,6) = ang([LN.t(iV,1) LN.t(iV,2) LN.t(iV,3)],[LN.b(iV,1) LN.b(iV,2) LN.b(iV,3)]);
    C(iV,4) = LN.ang(iV);
    C(iV,5) = LN.GND_per(iV);
    if C(iV,2) ~= 90 || C(iV,3) ~= 90
        OK(iV)=iV;
    end
    else
    OK(iV)=iV;      C(iV,5) = LN.GND_per(iV);
    end
end
OK(OK==0)=[];       C(OK,:)=[];