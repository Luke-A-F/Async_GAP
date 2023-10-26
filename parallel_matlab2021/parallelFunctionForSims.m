function [x_async,x_async_comparedtoOptimalOverTime] = parallelFunctionForSims(commRateIndex,c,A_assignment,Aineq_local,beq,bineq,EIPS,regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal)


if commRateIndex == 1
    commRate = 1;
elseif commRateIndex == 2
    commRate = 0.75;
elseif commRateIndex == 3
    commRate = 0.5;
elseif commRateIndex == 4
    commRate = 0.1;
end
[x_async, es_async,numDualUpdates_async,convDist_async,constr_async,x_async_comparedtoOptimalOverTime] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,compute_rate,commRate,plotVals,xIntegerOptimal);

end

