function [conmap,velop,recal] = judge_and_change(velop, ...
    q_new,q_old,conmap,F)
global N dt stack stack_bound thumb
%% add constraint
recal = 0;
for i = 1:N
    xpos = q_new(2*i-1);
    ypos = q_new(2*i);
    line = stack(1)*ypos^2+stack(2)*ypos+stack(3);
    if conmap(i,1) == 0
        if ypos < thumb
            fprintf('constraint added at node %d caused by thumb \n',i)
            conmap(i,1) =1;
            conmap(i,2) =1;
            velop(2*i-1,:) = (q_new(2*i-1)-q_old(2*i-1))/dt;
            velop(2*i,:) = (q_new(2*i)-q_old(2*i))/dt;
            recal = 1;
        elseif (ypos > stack_bound(1) && ypos < stack_bound(2) && xpos < line)
            fprintf('constraint added at node %d, caused by stack\n',i)
            conmap(i,1) =1;
            conmap(i,2) =0;
            velop(2*i-1,:) = (q_new(2*i-1)-q_old(2*i-1))/dt;
            velop(2*i,:) = (q_new(2*i)-q_old(2*i))/dt;
            recal = 1;
        end
    end
    % calculate the unit normal vector at contact point
end

%% remove constraint
for i = 1:N
    xpos = q_new(2*i-1);
    ypos = q_new(2*i);
    line = stack(1)*ypos^2+stack(2)*ypos+stack(3);
    if conmap(i,1) == 1
        frec = F(2*i-1:2*i);
        if conmap(i,2) == 1 % thumb cons 
            fnormal = dot(frec,[0,1]);
            if (fnormal <=0 && ypos >= thumb - 2e-6)
                fprintf('constraint deleted at node %d \n',i)
                conmap(i,1) = 0;
                recal = 1;
            end
        elseif conmap(i,2) == 0 % stack cons
            fnormal = dot(frec,nstack(ypos,xpos));
            if (fnormal <=0 && xpos >= line - 2e-5) || (ypos < stack_bound(1) && ypos > stack_bound(2))
                fprintf('constraint deleted at node %d \n',i)
                conmap(i,1) = 0;
                recal = 1;
            end
        end
    end

end


    
