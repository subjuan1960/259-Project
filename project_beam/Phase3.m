clc; clear; close all;
global free_index qx_comb qy_comb
Mat_prop;

q = zeros(2*N,1); 
time1 = 1.5;
free_index=1:2*N;
Phase2 = load('phase2.mat');
q_new = Phase2.q_new;
q_old = Phase2.q_old;
qd_old = Phase2.qd_old;
M = Phase2.M;
qd_new = qd_old;

qlist = [];
qdlist = [];
tlist = [];

for i = 1:(time1/dt)
    fprintf('time is %.4f s \n',i*dt)
    t= i*dt;
    err = 10;
    while err >1e-5
        e_dis = compute_external_force(i,q_new,q_old,qd_old,dt,time1);
        [F,J] = ComputeF_J(N,q_new,q_old,M,e_dis,dl,E,A,I);
        F = F+ M*(((q_new-q_old)/dt-qd_old)/dt);
        Jini = M*(1/dt)*(1/dt);
        Jdamp = dampratio/dt*eye(size(J));
        J = J + Jini+Jdamp;
        deltaX = J(free_index,free_index) \ F(free_index);
        q_new(free_index) = q_new(free_index) - deltaX;
        err = sum(abs(F(free_index)));
    end
    qd_new = (q_new - q_old) / dt;
    qlist = [qlist,q_new];
    qdlist = [qdlist,qd_new];
    tlist = [tlist,t];
    q_old = q_new;
    qd_old = qd_new;
end
%% Draw
qx = zeros(N,time1/dt);
qy = zeros(N,time1/dt);

for k =1:N
    qx(k,:) = qlist(2*k-1,:);
    qy(k,:) = qlist(2*k,:);
end

qx_comb = [qx_comb, qx];
qy_comb = [qy_comb, qy];

% figure(1)
% 
% for k = 1:i
%     plot(qx(:,k),qy(:,k),'ro-');
%     axis([-0.1,0.3,-0.1,0.3])
%     
%     %frame = getframe(gcf);
%     %writeVideo(v,frame);
%     drawnow
% end