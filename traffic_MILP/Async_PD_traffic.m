function [ZFinal,es,numDualUpdates,convDist,constr,ZOverTime] = Async_PD_traffic(cost,A,b,stepSizes,maxIterations,tolerance,lb_bin,ub_bin,lb_con,ub_con,lbDual,ubDual,alpha,delta,scBlock,primalNum,dualNum,commRate,compute_rate,plotVals,xStar,binIndex,Z_init,b_exit)
%Distributed Async Saddle Point Optimization
numDualUpdates = 0;
gamma = stepSizes(1);
beta = stepSizes(2);

Zs = Z_init; %Initial Primal  Variables
LamS = zeros(size(b,1),1);%Initial Dual Variables

Zp = [];
ZOverTime =[];
for j = 1:primalNum
    Zp = [Zp Zs]; %Primal Matrix initialization
end
Zd = [];
for j= 1:dualNum
    Zd = [Zd Zs];
end


dCount = zeros(primalNum,dualNum);%Initialize dCount vector
computeCount = zeros(primalNum,1);

k = 1;
convDist = tolerance+1;
constr = 100;


%Primal Block list
zBlocks = round(linspace(1,length(Zs),primalNum+1));%Plus one is to get the start and end indices for primal number of blocks
% %Dual Block list
lamBlocks = round(linspace(1,length(LamS),dualNum+1)); %Plus one is to get the start and end indices for dual number of blocks


while k<=maxIterations %convDist(k)> tolerance&&  %&& ~all(A*Zs<=b)%Stop when feasibility is reached
    prevIter = Zs;
    lamPrev = LamS;
%       for p = 1:primalNum
%           z = Zp(:,p);
%             %Lower bound block
%             aBl = zBlocks(p);
%             %upper bound block
%             bBl = zBlocks(p+1);
%             %Scalar Block Update
%             if scBlock
%                 for i=aBl:bBl
%                     %Primal Update
%                     gradZ = cost(i,:)+LamS(:,1)'*A(:,i)+alpha*z(i,1);
%                     pUpdate = z(i,1)-gamma*gradZ; %Primal gradient descent, projected onto Z
%                     if pUpdate >= ub
%                         pUpdate = ub;
%                     elseif pUpdate <= lb
%                             pUpdate = lb;
%                     end
%                     Zp(i,p) = pUpdate;
% 
%                     Zs(i,1) = Zp(i,p);%Same as doing diag(Zp)
%                 end
%             
%             else
%                 %Vector Block Update Primal Update
%                 gradZ = cost(aBl:bBl,:)+(LamS(:,1)'*A(:,aBl:bBl))'+alpha*z(aBl:bBl,1);
%                 pUpdate = Zs(aBl:bBl,1)-gamma*gradZ; %Primal gradient descent, projected onto Z
%       
% 
% 
%                 pUpdate(pUpdate>=ub) = ub;
%                 pUpdate(pUpdate<=lb) = lb;
%                           
% 
%                 Zp(aBl:bBl,p) = pUpdate;
% 
%                 Zs(aBl:bBl,1) = Zp(aBl:bBl,p);
%             end   
%        end

        %Compute Primal Updates
        [Zs,Zp,val_primalUpdate]= infrequent__primal_computation_MILP_traffic(Zp,Zs,primalNum,zBlocks,cost,LamS,A,alpha,gamma,ub_bin,lb_bin,ub_con,lb_con,scBlock,compute_rate,binIndex);


    
    %Communicate Primal Updates
    [Zp,Zd,dval] = infrequentComms_MILP(Zp,Zd,primalNum,dualNum,zBlocks,commRate);
  

    if k ==1
        dCount = dval;
        computeCount = val_primalUpdate;
    else
        
        dCount = dCount + dval;
        computeCount = computeCount + val_primalUpdate;
    end
    dCount_nz= sum(dCount~=0,1);
    computeCount_nz= sum(computeCount~=0,1);
    
     for d = 1:dualNum
          lamNew = [];
             if dCount_nz(1,d) >=primalNum && computeCount_nz>=primalNum
                 numDualUpdates = numDualUpdates + 1;
                %Lower bound block
                aBl = lamBlocks(d);
                %upper bound block
                bBl = lamBlocks(d+1);
                if scBlock
                    %Scalar Block Update
                    for j= aBl:bBl
                        %Dual Update
                        gradLam = A(j,:)*Zd(:,d)-b(j,1)-delta*LamS(j,:);
                        dualUpdate = LamS(j,1)+beta*gradLam; %Dual gradient ascent, projected onto M
                        if dualUpdate >= ubDual
                            dualUpdate = ubDual;
                        elseif dualUpdate <= lbDual
                                dualUpdate = lbDual;
                        end
                        lamNew = [lamNew; dualUpdate];
                    end
                    LamS(aBl:bBl) = lamNew;
                 else
                     %Vector Based Dual Update
                     gradLam = A(aBl:bBl,:)*Zd(:,d)-b(aBl:bBl,1)-delta*LamS(aBl:bBl,:);
                     dualUpdate = LamS(aBl:bBl,1)+beta*gradLam; %Dual gradient ascent, projected onto M
                     dualUpdate(dualUpdate>= ubDual) = ubDual;
                     dualUpdate(dualUpdate<=lbDual) = lbDual;
                     LamS(aBl:bBl) = dualUpdate;
                end
                
                dCount = zeros(primalNum,dualNum);
             end
     end   



  
    
   
    curIter = Zs;
    lamCur = LamS;
    k = k+1;
    convDist = [convDist;norm(prevIter-curIter,2)];
%     constr(k) = sum(A*Zs-b+rho2);
    % constr(k) = norm(A*Zs-b,1);
    %Rounding to compare current integer solutions
    %For PQP formulation
    % Ztmp = Zs;
    % Ztmp(binIndex+1:end)=round(Zs(binIndex+1:end));
    Ztmp = Zs;
    Ztmp(1:binIndex)=round(Zs(1:binIndex));
     % constr(k) = norm(A*Ztmp-b,inf);
     % constr(k) = length(A*Zs<=b_exit)-sum(A*Zs<=b_exit);
      % constr(k) = length(A*Ztmp<=b_exit)-sum(A*Ztmp<=b_exit);
     constr(k) = length(A*Zs<=b_exit)-sum(A*Zs<=b_exit);
     es(k) = abs(cost'*xStar-cost'*Ztmp)/abs(cost'*xStar); %Normalized system error, distance to Sol divided by norm of Sol
    ZOverTime = [ZOverTime;norm(Zs-xStar)];
    % ZOverTime = [ZOverTime;norm(round(Zs)-xStar)];
    
     %End if feasible
    if (length(A*Ztmp<=b_exit)-sum(A*Ztmp<=b_exit))== 0
        break;
    end
    
     
end

%Only round binary part
ZFinal = Zs;
ZFinal(1:binIndex)=round(ZFinal(1:binIndex));
%PQP example
%  ZFinal = Zs;
% ZFinal(binIndex+1:end)=round(ZFinal(binIndex+1:end));



% ZFinal = round(Zs);
end

