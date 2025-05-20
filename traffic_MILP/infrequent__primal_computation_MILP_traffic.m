function [Zs,Zp,val] = infrequent__primal_computation_MILP_erg(Zp,Zs,primalNum,zBlocks,cost,LamS,A,alpha,gamma,ub_bin,lb_bin,ub_con,lb_con,scBlock,compute_rate,binIndex)


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


                    % if aBl<binIndex && bBl <binIndex%Binary
                    %     pUpdate(pUpdate>=ub_bin) = ub_bin;
                    %     pUpdate(pUpdate<=lb_bin) = lb_bin;
                    % elseif aBl<binIndex && bBl >binIndex%Mix of Continous and binary variables
                    %     for index = aBl:bBl
                    %         if index<binIndex
                    %             p_index = (index-aBl)+1;
                    %             if pUpdate(p_index)>=ub_bin
                    %                 pUpdate(p_index) = ub_bin;
                    %             elseif pUpdate(p_index)<=lb_bin
                    %                 pUpdate(p_index) = lb_bin;
                    %             end
                    %         else
                    %         end
                    % 
                    %     end
                    % else%Continuous 
                    %     pUpdate(pUpdate>=ub_con) = ub_con;
                    %     pUpdate(pUpdate<=lb_con) = lb_con;
                    % end


                     if aBl<binIndex && bBl <binIndex%Binary
                         
                            pUpdate(pUpdate>=ub_bin) = ub_bin;
                            pUpdate(pUpdate<=lb_bin) = lb_bin;
                         
                    elseif aBl<binIndex && bBl >binIndex+1%Mix of Continous and binary variables
                        for index = aBl:bBl
                            if index<=binIndex
                                p_index = (index-aBl)+1;
                                if pUpdate(p_index)>=ub_bin
                                    pUpdate(p_index) = ub_bin;
                                elseif pUpdate(p_index)<=lb_bin
                                    pUpdate(p_index) = lb_bin;
                                end
                            else
                                p_index = (index-aBl)+1;
                                
                                % if pUpdate(p_index)>=ub_con(index)
                                %     pUpdate(p_index) = ub_con(index);
                                % elseif pUpdate(p_index)<=lb_con(index)
                                %     pUpdate(p_index) = lb_con(index);
                                % end

                                 if pUpdate(p_index)>=ub_con
                                    pUpdate(p_index) = ub_con;
                                elseif pUpdate(p_index)<=lb_con
                                    pUpdate(p_index) = lb_con;
                                end
                            end

                        end
                    else%Continuous
                        % for ind_con = aBl:bBl
                        %     p_index = (ind_con-aBl)+1;
                        %     if pUpdate(p_index)>=ub_con(ind_con)
                        %         pUpdate(p_index)=ub_con(ind_con);
                        %     end
                        %     if pUpdate(p_index)<=lb_con(ind_con)
                        %         pUpdate(p_index)=lb_con(ind_con);
                        % 
                        %     end
                            %If the continuous bounds are scalar
                            
                            pUpdate(pUpdate>=ub_con) = ub_con;
                            pUpdate(pUpdate<=lb_con) = lb_con;
                            % %Test projection. TODO refactor if it works
                            % if aBl>8 && bBl<8+9
                            %     for it=aBl:bBl
                            %         p_index = (it-aBl)+1;
                            %         %Lazy lower bound projection to test
                            %         % if pUpdate(p_index)<=-21.5
                            %         %     pUpdate(p_index) = -21.5;
                            %         % end
                            %     % pUpdate(pUpdate>=ub_con) = ub_con;
                            %     % pUpdate(pUpdate<=lb_con) = lb_con;
                            %     end
                            % end
                        % end
                    end
                    


                    Zp(aBl:bBl,p) = pUpdate;

                    Zs(aBl:bBl,1) = Zp(aBl:bBl,p);
                     %Indexing for when enough computations are done
                     %for dual update
                     val(p)=1;
                end
            end   
       end
    
         


end


