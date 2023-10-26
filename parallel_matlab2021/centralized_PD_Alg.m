function [xCur,lambdaCur,xList,lamList] = centralized_PD_Alg(xInit,lamInit,cost,Aconst,bConst,stepSizes,maxIterations,tolerance,lb,ub,lbDual,ubDual,xOpt,primalReg,dualReg)
xCur = xInit;
lambdaCur = lamInit;
curIt = 0;
stepPrimal = stepSizes(1);
stepDual = stepSizes(2);
xList = [];
lamList = [];
while curIt < maxIterations %&& norm(xCur-xPrev)<=tolerance % curIt>1)
   
    xPrev = xCur;
    lambdaPrev = lambdaCur;
    xCur = primalUpdate(xCur,stepPrimal,cost,lambdaCur,Aconst,lb,ub,bConst,primalReg);
    lambdaCur = dualUpdate(lambdaCur,xCur,Aconst,bConst,stepDual,lbDual,ubDual,dualReg);
    %Uncomment to run until feasible
     if norm(lambdaCur-lambdaPrev)<=tolerance && norm(xCur-xPrev)<=tolerance
        break;
     end

    %Uncomment to run until feasible
%     if all(Aconst*round(xCur)<=bConst)
%         break;
%      end
    
    %If the optimal solution is known
%     xList = [xList;sum(abs(round(xCur)-xOpt))]
    %Rounding to compare current integer solutions 
    xList = [xList;norm(xOpt-round(xCur))];
%     xList = [xList;norm(xOpt-round(xCur))];
    lamList = [lamList;sum(lambdaCur)];
    curIt = curIt +1;
end
xCur = round(xCur);
end

