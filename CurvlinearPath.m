% Path Co-ordinates for Dynamic Path Following
% Input: Cartesian Co-ordinates for Path Planning Algorithm
% Output: Heading, Curvature, Path Length


function [th,s,ds,dth,k]=CurvlinearPath(rx,ry)
  
    %Path parameter - u
    du=0.001;              
    u=0:du:(size(rx,2)-1)*du;
    u = u/(du*size(rx,2));

    s = u*0;  %length
    th = u*0; %heading
    k = u*0;  %curvature
    w = u*0;  %yaw rate
    vy = u*0; %lateral velocity
    % rates of change
    dx=u*0; 
    dy=u*0;
    ddx=u*0;
    ddy=u*0;
    ds=u*0;

    %evaluate first derivatives and length
    for i=2:1:size(u,2)-1;
        th(i)=qtan(ry(i),rx(i),ry(i+1),rx(i+1));
        dx(i)=(rx(i+1)-rx(i));
        dy(i)=(ry(i+1)-ry(i));
        ds(i)=sqrt(dx(i)^2+dy(i)^2);
        s(i)=s(i-1)+ds(i);
    end

    %finalized end values
    dx(end)=dx(end-1);
    dy(end)=dy(end-1);
    th(1)=th(2);
    th(end)=th(end-1);
    ds(end)=sqrt(dx(end)^2+dy(end)^2);
    s(end) =s(end-1)+ds(end);
    dx(1)=dx(2);
    dy(1)=dy(2);

    %evaluate second derivatives
    for i=1:1:size(u,2)-1;
        ddy(i)=(dy(i+1)-dy(i));
        ddx(i)=(dx(i+1)-dx(i));
        dth(i)=(th(i+1)-th(i));
    end

    %finalize end values
    dth(end)=dth(end-1);
    ddy(end-1)=ddy(end-2);
    ddx(end-1)=ddx(end-2);
    ddy(end)=ddy(end-1);
    ddx(end)=ddx(end-1);

    %evaluate curvature
    k    =    (dx.*ddy-ddx.*dy)./((dx.^2+dy.^2).^(3/2));


end