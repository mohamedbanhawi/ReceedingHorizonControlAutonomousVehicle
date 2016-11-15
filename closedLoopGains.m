
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