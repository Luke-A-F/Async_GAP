function [Zp_new,Zd_new,dval] = infrequentComms_MILP(Zp,Zd,primalNum,dualNum,zBlocks,commRate)


    b_upper = rand(primalNum,primalNum+dualNum);
    Z_new = [Zp,Zd];
    Z_newnew = [Zp,Zd];
    dval = zeros(primalNum,dualNum);
    
    for i = 1:primalNum
            %Lower bound block
            aBl = zBlocks(i);
            %upper bound block
            bBl = zBlocks(i+1);
            for j=1:primalNum+dualNum
%             for j=1:primalNum
                if i~=j
                    if b_upper(i,j)<= commRate && commRate ~=0
                        b_upper(i,j)=1;
                        [~,maxLoc] = max(Zp(aBl:bBl,i));
                        newVec = zeros((bBl-aBl)+1,1);
                        newVec(maxLoc) = 1;
                        Z_new(aBl:bBl,j)= Zp(aBl:bBl,i);
                        Z_newnew(aBl:bBl,j) = newVec;
                        
                        temp = (j-primalNum);
                        if temp == 0
                            temp = 1;
                        end
                        if j>= primalNum %&& temp <= dualNum
                            
                            dval(i,temp)=1;
                        end
        

                end
            end
    end
    Zp_new = Z_new(:,1:primalNum);
    Zd_new = Z_new(:,primalNum+1:primalNum+dualNum);

    
         


end


