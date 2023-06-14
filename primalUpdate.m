function x = primalUpdate(x,stepSize,cost,lambda,Aconst,lb,ub,bConst,primalReg)
%Specifically for linear cost

%     x = x - stepSize*(cost'+(lambda'*Aconst))';
    %Strict convexity penalty term
%     Twopsi = 2*10^-2;
    %Derivative of max(0,(Ax-b)^2);
    %Adding penalty term to force it to be strictly convex for convergence
    %as outlined in 2010 here Stability of primal–dual gradient dynamics and applications to
%     %network optimization which extends methods from Arrow Hurwitz
%     penaltyTerm = Twopsi*(Aconst*x-bConst);
%     penaltyTerm(penaltyTerm<=0)=0;
%     x = x - stepSize*(cost'+((lambda+penaltyTerm)'*Aconst))';


%     regVal =10^-3;
    regVal = primalReg;

%     %Tik Reg for primal variable
    x = x - stepSize*(cost'+(lambda'*Aconst)+regVal*x')'; 
%     
    x(x<=lb)=lb;
    x(x>=ub)=ub;
    
end

