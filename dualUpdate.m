function lam = dualUpdate(lam,xCur,Aconst,bConst,stepSize, lb,ub,dualReg)
%Specifically for affine constraints
%     lamPrev = lam;
%     gradLambda =(Aconst*xCur - bConst);
    %With Regularizer
%     regVal =10^-1;
    delta = dualReg;
    
    gradLambda =(Aconst*xCur - bConst)-delta*lam;
    stepSize = delta/(delta ^ (2) + 1);
    lam = lam + stepSize*gradLambda;
    lam(lam<=lb)=lb;
    lam(lam>=ub)=ub;
end

