clear;clc

%% par 
pars = struct();
pars.lam   = 0.15;    pars.alpha = 1;  pars.beta = 0.01;
pars.c     = 0.1;     pars.a     = 0.7;  pars.b    = 0.8;
pars.gamma = 0.1;     pars.k    = 0.1;  pars.G     = 0.0;
%  λ u_PC = A cos(ω t) + B cos(N ω t)
pars.A = 1.0; pars.B = 0.0; pars.omega = 1;  
pars.N = 0;  
%% 
T =8.4:0.048:18; A =0.1:0.0021:0.52;
%G_min=0.0001; G_max=0.1;
%G = G_min * ((G_max/G_min).^(1/39)) .^(0:39);

%%
tic
for jj = 1:length(A)
    pars.A=A(jj);
    parfor j = 1:length(T)
        Pars(j)=pars;
        Pars(j).omega = 2*pi/T(j);
        
        [~,~,~,meanE(j,jj)] = time_trial(Pars(j));
        [mLCE(j,jj)] = Lyapunov_trial(Pars(j))
    end
    toc
end

save 31102AT.mat
%% 
figure(101); 
heatmap(meanE)
figure(102); 
heatmap(mLCE)
mmLCE=mLCE;
mmLCE(mLCE<=0)=0; mmLCE(mLCE>0)=1;
figure(103); 
heatmap(mmLCE)

figure(100); 
subplot(2,1,1); plot(T, meanE,'LineWidth',1.5); grid on; axis([-inf inf,-0.1 2]);
ylabel('E(t) = ||[x_1-x_2, y_1-y_2]||');

% figure(1); 
% subplot(2,1,1); plot(t, sita,'LineWidth',1.5); grid on; axis([-inf inf,-0.1 1]);
% ylabel('E(t) = ||[x_1-x_2, y_1-y_2]||');
% 
% figure(2); 
% subplot(3,1,1); plot(t, x(:,1),'LineWidth',1.5);  axis([-inf inf,-2 2]); grid on; ylabel('x1(t)');
% subplot(3,1,2); plot(t, x(:,3),'LineWidth',1.2); axis([-inf inf,-2 2]); grid on; ylabel('x2(t)');
% subplot(3,1,3); plot(t, x(:,1), t, x(:,3),'LineWidth',1.2); axis([-inf inf,-2 2]); grid on; xlabel('t'); ylabel('x2+x2');



