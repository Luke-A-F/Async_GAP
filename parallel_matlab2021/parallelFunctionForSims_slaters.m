function [x_async,x_async_comparedtoOptimalOverTime] = parallelFunctionForSims_slaters(commRateIndex,c,A_assignment,Aineq_local,beq,bineq,EIPS,regCorrectionTest,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xIntegerOptimal)

commRate = 0.75;
compute_rate = 1;
if commRateIndex == 1
    slater_test = 0;
elseif commRateIndex == 2
    slater_test = 1;
elseif commRateIndex == 3
    slater_test = 2;
elseif commRateIndex == 4
    slater_test = 3;
end
[x_async, es_async,numDualUpdates_async,convDist_async,constr_async,x_async_comparedtoOptimalOverTime] = Async_PD(c,[A_assignment;Aineq_local],[floor(beq);floor(bineq)]+EIPS-regCorrectionTest+slater_test,stepSize,maxIterations,tolerance,lbPrimal,ubPrimal,lbDual,ubDual,primalReg,dualReg,scBlock,primalNum,dualNum,compute_rate,commRate,plotVals,xIntegerOptimal);

end