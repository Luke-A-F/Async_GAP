function [Zs,Zp,val] = infrequent__primal_computation_MILP(Zp,Zs,primalNum,zBlocks,cost,LamS,A,alpha,gamma,ub,lb,scBlock,compute_rate)


    b_upper = rand(primalNum,1);
    val = zeros(primalNum,1);
    
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
                    if b_upper(p)<=compute_rate && compute_rate ~=0
                        b_upper(p) = 1
                        gradZ = cost(i,:)+LamS(:,1)'*A(:,i)+alpha*z(i,1);
                        pUpdate = z(i,1)-gamma*gradZ; %Primal gradient descent, projected onto Z
                        if pUpdate >= ub
                            pUpdate = ub;
                        elseif pUpdate <= lb
                                pUpdate = lb;
                        end
                        Zp(i,p) = pUpdate;

                        Zs(i,1) = Zp(i,p);%Same as doing diag(Zp)
                        
                       
                         %Indexing for when enough computations are done
                         %for dual update
                         val(i)=1;
                       
                    end
                end
            
            else
                %Vector Block Update Primal Update
                if b_upper(p)<=compute_rate && compute_rate ~=0
                    gradZ = cost(aBl:bBl,:)+(LamS(:,1)'*A(:,aBl:bBl))'+alpha*z(aBl:bBl,1);
                    pUpdate = Zs(aBl:bBl,1)-gamma*gradZ; %Primal gradient descent, projected onto Z



                    pUpdate(pUpdate>=ub) = ub;
                    pUpdate(pUpdate<=lb) = lb;


                    Zp(aBl:bBl,p) = pUpdate;

                    Zs(aBl:bBl,1) = Zp(aBl:bBl,p);
                     %Indexing for when enough computations are done
                     %for dual update
                     val(p)=1;
                end
            end   
       end
    
         


end


