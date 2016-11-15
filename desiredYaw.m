% Path Co-ordinates for Dynamic Path Following
%
% Input:-curvlinear reference path (ds,th)
%        -cartesian coordinates rx.ry
%       -traction speed vx
%       -current position x,y
%
% Output: Heading, Curvature, Path Length


function [ssW,ssPhi]=desiredYaw(th,s,ds,dth,k,sK,sYaw,sBeta,vx,rx,ry,xi,yi,thi,dist,dlook,Betamax)

    %nearest point
    dx=xi-rx(end);
    dy=yi-ry(end);
    dmin = sqrt(dx*dx+dy*dy); 
    amin = size(rx,2);
    %lookahead distance
    for kk=0:1:size(rx,2)-1;
        dx=xi-rx(end-kk);
        dy=yi-ry(end-kk);
        d = sqrt(dx*dx+dy*dy);
        if d<dlook && s(end-kk)>=dist
            amin=size(rx,2)-kk;
            dmin = d;
            break
        end
    end
    
    %steady state steering based on curvature response
    ssPhi = th(amin)-thi;
    %estimate steady state Beta
    ssBeta = sBeta*ssPhi;
    if abs(ssBeta) > Betamax
        ssBeta = Betamax*sign(ssBeta);
        ssPhi = ssBeta/sBeta;
    end
    %steady state yaw based on steady state steering angle corrected with
    ssW = ssPhi*sYaw*cos(ssBeta);
end