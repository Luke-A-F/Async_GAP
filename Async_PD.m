function [ZFinal,es,numDualUpdates,convDist,constr,ZOverTime] = Async_PD(cost,A,b,stepSizes,maxIterations,tolerance,lb,ub,lbDual,ubDual,alpha,delta,scBlock,primalNum,dualNum,commRate,plotVals,xStar)
%Distributed Async Saddle Point Optimization
numDualUpdates = 0;
gamma = stepSizes(1);
beta = stepSizes(2);

Zs = ones(size(cost,1),1); %Initial Primal  Variables
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
      for p = 1:primalNum
          z = Zp(:,p);
            %Lower bound block
            aBl = zBlocks(p);
            %upper bound block
            bBl = zBlocks(p+1);
            %Scalar Block Update
            if scBlock
                for i=aBl:bBl
                    %Primal Update
                    gradZ = cost(i,:)+LamS(:,1)'*A(:,i)+alpha*z(i,1);
                    pUpdate = z(i,1)-gamma*gradZ; %Primal gradient descent, projected onto Z
                    if pUpdate >= ub
                        pUpdate = ub;
                    elseif pUpdate <= lb
                            pUpdate = lb;
                    end
                    Zp(i,p) = pUpdate;

                    Zs(i,1) = Zp(i,p);%Same as doing diag(Zp)
                end
            
            else
                %Vector Block Update Primal Update
                gradZ = cost(aBl:bBl,:)+(LamS(:,1)'*A(:,aBl:bBl))'+alpha*z(aBl:bBl,1);
                pUpdate = Zs(aBl:bBl,1)-gamma*gradZ; %Primal gradient descent, projected onto Z
      


                pUpdate(pUpdate>=ub) = ub;
                pUpdate(pUpdate<=lb) = lb;
                          

                Zp(aBl:bBl,p) = pUpdate;

                Zs(aBl:bBl,1) = Zp(aBl:bBl,p);
            end   
       end

    
    %Communicate Primal Updates
    [Zp,Zd,dval] = infrequentComms_MILP(Zp,Zd,primalNum,dualNum,zBlocks,commRate);
  

    if k ==1
        dCount = dval;
    else
        
        dCount = dCount + dval;
    end
    dCount_nz= sum(dCount~=0,1);

     for d = 1:dualNum
          lamNew = [];
             if dCount_nz(1,d) >=primalNum 
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



  
    es(k) = cost'*Zs; %Normalized system error, distance to Sol divided by norm of Sol
   
    curIter = Zs;
    lamCur = LamS;
    k = k+1;
    convDist = [convDist;norm(prevIter-curIter,2)];
%     constr(k) = sum(A*Zs-b+rho2);
    constr(k) = sum(A*Zs-b);
    %Rounding to compare current integer solutions
%     ZOverTime = [ZOverTime;norm(round(Zs)-xStar)];
    ZOverTime = [ZOverTime;norm(round(Zs)-xStar)];
    %Uncomment to run until feasible
     if norm(lamCur-lamPrev)<=tolerance && norm(curIter-prevIter)<=tolerance
        break;
     end
    %Uncomment to run until feasible
%     if all(A*round(Zs)<=b)
%         break;
%      end
    
     
end


% % % LamS
% if plotVals
%     figure()
% 
%     hold on
%     semilogy(es(2:end)','-','LineWidth',2)
%     title('Cost')
%     xlabel('Iteration Number','FontWeight','Bold')
%     ylabel('System Cost','FontWeight','Bold')
%     hold off
% 
%     figure ()
%     % % 
%     hold on
%     semilogy(convDist(2:end)','-','LineWidth',2)
%     title('Distance Between Iterates Convergence Comparison')
%     xlabel('Iteration Number','FontWeight','Bold')
%     ylabel('Distance Between Iterates','FontWeight','Bold')
%     hold off
% 
% 
% 
% 
%     figure()
%     % % 
%     hold on
%     semilogy(constr(2:end)','-','LineWidth',2)
%     title('Constraint Satisfaction Relative to Iteration Number')
%     xlabel('Iteration Number','FontWeight','Bold')
%     ylabel('Constraints','FontWeight','Bold')
%     hold off
%     
% %     if numDualUpdates ~= k
% %         figure()
% %         % % 
% %         hold on
% %         conversion = round(k/numDualUpdates);
% %         constrDuals = constr(2:conversion:end);
% %         semilogy(constrDuals','-','LineWidth',2)
% %         title('Constraint Satisfaction Relative to Dual Updates')
% %         xlabel('Dual Updates','FontWeight','Bold')
% %         ylabel('Constraints','FontWeight','Bold')
% %         hold off
% % 
% % 
% 
%         figure()
%         % % 
%         hold on
%         conversion = round(k/numDualUpdates);
%         distDuals =convDist(2:conversion:end);
%         semilogy(distDuals','-','LineWidth',2)
%         title('Distance Between Iterates Relative to Dual Updates')
%         xlabel('Dual Updates','FontWeight','Bold')
%         ylabel('Distance Between Iterates','FontWeight','Bold')
%         hold off
%     end
    
% end

ZFinal = round(Zs);
end

