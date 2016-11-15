% Lateral RHC path tracking algorithm for car-like FWS robot
% Submitted to Journal of Vibration and Control
%
% Written by Mohamed Elbanhawi [2015]
% email: mohamed.elbenhawi@rmit.edu.au
% RMIT University
%
%-------------------------------------------------------------------------%
%-Load Reference Path
%-Modify Vehicle Parameters
%-Modify appropriate dlook (PP-lookahead) based on Snider (2009).

function Main_RHC_Tracking
clc;clear;
%Load reference Path
load('xyRefStep.mat');

%simulation timestep
dt =0.01;% sec

%Vehicle Parameters
vx  =10;  % m/s    [Longitudinal Velocity]
cf =3000; % N/rad  [Front wheel coefficient]
cr =3000; % N/rad  [Rear wheel coefficient]
a1 =1.0;  % m      [Front to CG distance]
a2 =1.6;  % m      [Rear to CG distance]
L  =2.6;  % m      [Wheel Base]
Iz =1650; % Kg.m^2 [Moment of Interia]
m  =1000; % Kg     [Mass]

umax = 30*pi/180; % maximum steering angle
umin =-30*pi/180; % minimum steering angle
    
%force system coefficients
cbeta = - (cf+cr);
cphi  = cf;
    
dbeta = -(a1*cf-a2*cr);
dphi  = a1*cf;

%velocity dependent force system coefficients
cyaw = -a1*cf/vx+a2*cr/vx; 
dyaw = -a1*a1*cf/vx-a2*a2*cr/vx;

%steady state responses
%yaw rate
sYaw   = (cphi*dbeta-cbeta*dphi)/(dyaw*cbeta-cyaw*dbeta+m*vx*dbeta);
%Curvature
sK     = sYaw/vx;
%side slip
sBeta   = (dphi*(cyaw-m*vx)-dyaw*cphi)/(dyaw*cbeta-cyaw*dbeta+m*vx*dbeta);
%maximum side slip angle
Betamax = 10-7*vx*vx/(40*40);

%Lateral Control Model: time invariant model fixed longitudinal velocity
Ac =[-(cf+cr)/(m*vx),(-a1*cf+a2*cr)/(m*vx*vx)-1;(-a1*cf+a2*cr)/Iz,-(a1*a1*cf+a2*a2*cr)/(Iz*vx)];
Bc =[cf/(m*vx);a1*cf/Iz];
Cc =[0,1];
Dc = 0;

%discretize model
[Ap,Bp,Cp,~]=c2dm(Ac,Bc,Cc,Dc,dt);

%Initial conditions 
%MPC Formulations
xm = [0;0]; 
ym = 0;
um = 0;
%Kinematic state
xi = xm(1);
yi = xm(2);
thi = qtan(ry(1),rx(1),ry(2),rx(2));

%Control and prediction horizon Nc<=Np
Nc =15;
Np =50;
    
%Closed Loop tuning parameter (proportional to U and inversely proportional to error signal)
rw = 100;

%MPC gains
[Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e,BarRs]=mpcgain(Ap,Bp,Cp,Nc,Np);

I = eye(Nc);
BarR = rw*I;

%Initial augmented state
n=size(B_e,1);
X_k(:,1)=zeros(n,1);

%Linear Inequality Constrained Matrices
[Ky,Kmpc]=closedLoopGains(Phi_Phi,Phi_F,Phi_R,rw,Nc);

%------------Receeding Horizon Implementation-----------------------------%

%Vehicle State
xhat=[];
yhat=[];
that=[];
Beta=[];
xhat=[xhat,xi];
yhat=[yhat,yi];
that=[that,thi];
Beta=[Beta,0];
dist=0;

%Controller Parameters
phi(1)=0; % steering angle
r_k(1)=0; % reference yaw
t(1)=0;   % time 
HError(1)=0;


% PP lookahead distance!! needs to be modified
if vx>=20;
    dlook=10;
elseif vx>10;
    dlook =10; 
elseif vx<5 && vx>2;
    dlook=(a1+a2)*2;
else
    dlook=2;
end

dlook=5;

%index
kk=1;

%convert from cartesian to curvlinear co-ordinates
[th,s,ds,dth,k]=CurvlinearPath(rx,ry);

%distance-to-go
dtg = Inf;
dist = 0;

while (dtg>0.3)
    
    %Reference Yaw set point
    [r_k(kk),HError(kk+1)]=desiredYaw(th,s,ds,dth,k,sK,sYaw,sBeta,vx,rx,ry,xhat(end),yhat(end),that(end),dist,dlook,Betamax);
    
    %RHC 
    %Optimisation Matrices
    dU=closedLoopDeltaU(Ky,Kmpc,X_k(:,kk),r_k(kk));

    %implement only first control set
    um = um+dU(1);
    
    %Soft Constraints on steering Input
    um = min(um,umax);
    um = max(um,umin);

    %store control signal history
    phi(kk) = um;

    %update state space model using control
    xm_old = xm;
    xm=Ap*xm+Bp*um;
    ym=Cp*xm;

    %update augemented state space
    X_k(:,kk+1) = [xm-xm_old;ym];
            
    %side slip angle
    Beta(kk+1) = xm(1);
    
    %position update
    xhat=[xhat,xhat(kk)+vx*cos(that(kk)+Beta(kk))*dt];
    yhat=[yhat,yhat(kk)+vx*sin(that(kk)+Beta(kk))*dt];
    that=[that,that(kk)+vx*tan(phi(kk))*dt/L]; 
    
    %update distance-to-go and distance travelled
    dtg = sqrt((xhat(end)-rx(end))^2+(yhat(end)-yhat(end))^2);
 
    dist = dist+ sqrt((xhat(kk+1)-xhat(kk))^2+(yhat(kk+1)-yhat(kk))^2); %travelled distance
        
    %update step
    kk=kk+1;
    
    %time update
    t(kk) = t(kk-1)+dt;
    
end

%Path
figure(1);
%Show Results
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultAxesFontSize', 16)
subplot(3,2,[1,2]);plot(xhat,yhat,'color',[34,147,216]/255,'LineWidth',2);hold on
subplot(3,2,[1,2]);plot(rx,ry,'--','color',[165,165,165]/255,'LineWidth',2);axis('equal');legend('Actual','Reference');xlabel('x');ylabel('y'); hold on;
%yaw rate
subplot(3,2,[3,4]);plot(t(1:end-1),r_k,'--','color',[165,165,165]/255,'LineWidth',2);xlabel('Time [sec]');hold on 
subplot(3,2,[3,4]);plot(t,X_k(3,:),'color',[34,147,216]/255,'LineWidth',2);ylabel('Yaw Rate [deg/s]');legend('Actual','Reference');hold on %yaw rate
%Steering angle control
subplot(3,2,5);plot(t(1:end-1),umax*ones(size(phi))*180/pi,'--','color',[165,165,165]/255,'LineWidth',2);hold on;
subplot(3,2,5);plot(t(1:end-1),umin*ones(size(phi))*180/pi,'--','color',[165,165,165]/255,'LineWidth',2);hold on;
subplot(3,2,5);plot(t(1:end-1),phi*180/pi,'color',[216,41,46]/255,'LineWidth',2);ylabel('Steering Angle [deg]'); xlabel('Time [sec]');legend('Maximum','Minimum');hold on; %steering angle
xlim([0,t(end)]);
%Velocity
subplot(3,2,6);plot(t,vx*ones(size(t)),'--','color',[216,41,46]/255,'LineWidth',2); hold on;ylabel('Velocity [m/s]'); xlabel('Time [sec]');
xlim([0,t(end)]);

figure(2);
%Show Results
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultAxesFontSize', 16)
%Side slip
plot(t,Betamax*ones(size(t)),'--','color',[165,165,165]/255,'LineWidth',2); hold on;
plot(t,-Betamax*ones(size(t)),'--','color',[165,165,165]/255,'LineWidth',2); hold on;
plot(t,Beta*180/pi,'color',[216,41,46]/255,'LineWidth',2);ylabel('Side Slipe Angle [deg]');xlabel('Time [sec]');hold on; %slip angle
end


%calculated closed loop control for RHC of time invariant system
%given closed loop gain

function dU=closedLoopDeltaU (Ky,Kmpc,X_k,r_k)

%closed loop receeding horizon control for step k

dU = Ky*r_k-Kmpc*X_k; 

end


%calculate closed loop gains for one step of RHC of time invariant system
%given MPC gains

function [Ky,Kmpc]=closedLoopGains (Phi_Phi,Phi_F,Phi_R,rw,Nc)

%Tuning parameter diagonal matrix
I = eye(Nc);
BarR = rw*I;

%Set point change gain
Ky = (inv(Phi_Phi+BarR))*Phi_R;
Ky = Ky(1,1);

%State feedback controller gain
Kmpc = (inv(Phi_Phi+BarR))*Phi_F;
Kmpc = Kmpc(1,:);


end

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

function [Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e,BarRs]=mpcgain(Ap,Bp,Cp,Nc,Np)

%Augemented model (Ae,Be,Ce) from state space (Ap,Bp,Cp) model
[m1,n1]=size(Cp);
[n1,n_in]=size(Bp);
A_e=eye(n1+m1,n1+m1);
A_e(1:n1,1:n1)=Ap;
A_e(n1+1:n1+m1,1:n1)=Cp*Ap;
B_e=zeros(n1+m1,n_in);
B_e(1:n1,:)=Bp;
B_e(n1+1:n1+m1,:)=Cp*Bp;
C_e=zeros(m1,n1+m1);
C_e(:,n1+1:n1+m1)=eye(m1,m1);

%Np is prediction horizon
%Nc is control horizon
%Nc <= Np

n=n1+m1;
h(1,:)=C_e;
F(1,:)=C_e*A_e;

for kk=2:Np
    h(kk,:)=h(kk-1,:)*A_e;
    F(kk,:)= F(kk-1,:)*A_e;
end
v=h*B_e; % first column of Phi
Phi=zeros(Np,Nc); %declare the dimension of Phi
Phi(:,1)=v; % first column of Phi

for i=2:Nc
    Phi(:,i)=[zeros(i-1,1);v(1:Np-i+1,1)]; %Toeplitz matrix
end

BarRs=ones(Np,1);
Phi_Phi= Phi'*Phi;
Phi_F= Phi'*F;
Phi_R=Phi'*BarRs;

end

function theta=qtan(y1,x1,y2,x2)
    
theta=atan((y2-y1)/(x2-x1));
    
%90 degrees
if (x2==x1)&&(y2>y1);
    theta = pi/2;
   
%180 degrees
elseif (x2<x1)&&(y2==y1);
    theta = pi;
    
%270 degrees
elseif (x2==x1)&&(y2<y1);
    theta =-pi/2;  
    
%Zero
elseif (x2>x1)&&(y2==y1);
    theta =0;
    
    %First Quadrant
elseif (x2>x1)&&(y2>y1);
    theta =abs(theta);
    
%Second Quadrant
elseif (x2<x1)&&(y2>y1);
    theta = pi - abs(theta);
    
%Third Quadrant
elseif (x2<x1)&&(y2<y1);
    theta = -pi + abs(theta);
    
%Fourth Quadrant
elseif (x2>x1)&&(y2<y1);
    theta = -abs(theta); 
    
end

end

