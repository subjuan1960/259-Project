Phase1
Phase2
Phase3

time = 1.5+10+1.5;
dt = 0.001;
i = time/dt;
%%
figure(1)
for k = 1:i
    plot(qx_comb(:,k),qy_comb(:,k),'ro-');
    axis([-0.1,0.3,-0.1,0.7])
    
    %frame = getframe(gcf);
    %writeVideo(v,frame);
    drawnow
end