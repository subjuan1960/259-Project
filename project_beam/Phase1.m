clc;clear;close all;
global free_index qx_comb qy_comb
Mat_prop;
% N = 50;
% l = 0.14605; % length
% dl = l/(N-1);
% r0 =  2.9074e-04/2; % radius of the rod
% rou = 436; % density kg/m3
% E = 1e9; % Young's Modulus
% I = pi/4*(r0^4);
% A = pi*r0^2;
% rou_al = 436;
% m  = pi*(r0^2)*l*rou_al/(N-1);
% dt = 0.001;
q=zeros(2*N,1);
cood_ini = linspace(0,l,N);
for i = 1:N
    xk = [0,cood_ini(i)];
    q(2*i-1:2*i) = xk';
end
M = eye(N*2)*m;
time1 = 1.5;
free_index=2:2*(N-1);
% dampratio = 0.00001;
qlist = [];
qdlist= [];
tlist = [];
q_new = q;
q_old = q;
qd_old = zeros(2*N,1);
%% Phase1
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
%%
save('phase1.mat','q_new','q_old','qd_old','M')
%%
qx = zeros(N,time1/dt);
qy = zeros(N,time1/dt);
for k =1:N
    qx(k,:) = qlist(2*k-1,:);
    qy(k,:) = qlist(2*k,:);
end
qx_comb = [];
qy_comb = [];
qx_comb = [qx_comb, qx];
qy_comb = [qy_comb, qy];
%v = VideoWriter('Output0.avi');
%open(v);


% figure(1)
% %set(gcf,'outerposition',get(0,'screensize'));
% for k = 1:i
%     plot(qx_comb(:,k),qy_comb(:,k),'ro-');
%     axis([-0.1,0.1,-0.1,0.2])
%     
%     %frame = getframe(gcf);
%     %writeVideo(v,frame);
%     drawnow
% end