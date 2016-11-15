
%calculated closed loop control for RHC of time invariant system
%given closed loop gain

function dU=closedLoopDeltaU (Ky,Kmpc,X_k,r_k)

%closed loop receeding horizon control for step k

dU = Ky*r_k-Kmpc*X_k; 

end