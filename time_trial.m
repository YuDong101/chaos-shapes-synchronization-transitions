function [t,output_x,sita,meanE] = time_trial(pars)

dt = 1e-2;                    
T_init = 500;  T_cal = 2000;  
x0 = [0.2 0.1 0.02 0.01 0.1];

T_plot = T_init+T_cal;
N = round(T_plot/dt);
x = x0; output_x=0;
t = dt:dt:T_plot;
sita = zeros(N,1);
%TTT=NaN;
for ii = 1:N
    [x, ~] = rk4_step_with_mid(x, t(ii), dt, @(x,t,pp) f_model(x,t,pp), pars);
    % output_x(ii,:)=x;
    sita(ii) = hypot(x(1)-x(3), x(2)-x(4)); % sqrt((x1-x2)^2+(y1-y2)^2)
        %if isnan(TTT)&& t(ii)>T_init && sita(ii)>1
       % TTT=t(ii);
       
   % end
end

meanE=mean(sita(round(T_init/dt):round(T_plot/dt)));
sita=0;

function dx = f_model(x, t, pars)
% x = [x1; y1; x2; y2; z]
    x1 = x(1); y1 = x(2); x2 = x(3); y2 = x(4); z = x(5);

  
    lam   = pars.lam;        
    alpha = pars.alpha;
    beta  = pars.beta;
    c     = pars.c;
    a     = pars.a;
    b     = pars.b;
    gamma = pars.gamma;
    k     = pars.k;
    G     = pars.G;
    % =====   λ u_PC = A cos(ω t) + B cos(N ω t) =====
    forc = pars.A*cos(pars.omega*t) + pars.B*cos(pars.N*pars.omega*t);

    dx1 = forc +(1-k)*(x1*(1 - lam) - (1/3)*x1^3 - y1) +(2*k-1)*alpha*sin(z) +(2*k-1)*beta*(x1 - x2)+ k*(x2*(1-lam) - (1/3)*x2^3 - y2);
    dy1 = c*(x1 + a - b*y1);
    dx2 = forc +(1-k)*( x2*(1 - lam) - (1/3)*x2^3 - y2) -(2*k-1) *alpha*sin(z) -(2*k-1)* beta*(x1 - x2) +k*(x1*(1 - lam) - (1/3)*x1^3 - y1);
    dy2 = c*(x2 + a - b*y2);
    dz  = gamma*(x1 - x2)+G;

    dx = [dx1; dy1; dx2; dy2; dz];
end

function J = Df_model(x, t, pars) %#ok<INUSD>
    x1 = x(1); y1 = x(2); x2 = x(3); y2 = x(4); z = x(5);
    lam = pars.lam; alpha = pars.alpha; beta = pars.beta;
    c   = pars.c;   b     = pars.b;     gamma = pars.gamma; k = pars.k;
    G     = pars.G;

    J = zeros(5,5);
    %  df1/d*
    J(1,1) = (1-k)*((1-lam) - x1^2) + (2*k - 1)*beta;           
    J(1,2) = -(1-k);                                          
    J(1,3) = k*((1-lam) - x2^2) - (2*k - 1)*beta;              
    J(1,4) = -k;                                               
    J(1,5) = (2*k - 1)*alpha*cos(z);                          

    %  df2/d*
    J(2,1) =  c;
    J(2,2) = -c*b;

    %  df3/d*
    J(3,1) =  beta-k*(2*beta-(1-lam)+x1^2);
    J(3,2) =  -k;
    J(3,3) = (1 - lam) - x2^2 - beta-k*((1-lam)-x2^2-2*beta);
    J(3,4) = -1+k;
    J(3,5) =  (1-2*k)*alpha*cos(z);

    %  df4/d*
    J(4,3) =  c;
    J(4,4) = -c*b;

    %  df5/d*
    J(5,1) =  gamma;
    J(5,3) = -gamma;
end

function DH = DH_zero(x, t, pars) %#ok<INUSD>
    DH = zeros(5,5);
end

function mLCE = msf_mLCE(x0, f, Df, DH, dt, T_init, T_cal, gamma, varargin)
% 
% 

2% =
% x0 
% f  
% Df: @(x,t,...) 
% DH: @(x,t,...) 
% dt, T_init, T_cal
% gamma
% varargin

    x = x0(:);
    n = numel(x);
    t = 0.0;

     
    W = randn(n,1);
    W = W / norm(W);

    n_forward = max(1, round(T_init / dt));
    n_compute = max(1, round(T_cal  / dt));

   
    for k = 1:n_forward
       
        [x_full, ~] = rk4_step_with_mid(x, t, dt, f, varargin{:});
        x = x_full;  t = t + dt;

        
        [x_full2, x_half] = rk4_step_with_mid(x, t, dt, f, varargin{:}); %#ok<ASGLU>
        J0    = Df(x,        t,          varargin{:}) - gamma * DH(x,        t,          varargin{:});
        Jhalf = Df(x_half,   t+0.5*dt,   varargin{:}) - gamma * DH(x_half,   t+0.5*dt,   varargin{:});
        J1    = Df(x_full,   t+dt,       varargin{:}) - gamma * DH(x_full,   t+dt,       varargin{:});

        k1 = J0    * W;
        k2 = Jhalf * (W + 0.5*dt*k1);
        k3 = Jhalf * (W + 0.5*dt*k2);
        k4 = J1    * (W + dt*k3);
        W  = W + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

      
        W  = W / norm(W);
    end

   
    sum_log = 0.0;
    for k = 1:n_compute
        
        [x_full, ~] = rk4_step_with_mid(x, t, dt, f, varargin{:});
        x = x_full;  t = t + dt;

        
        [x_full2, x_half] = rk4_step_with_mid(x, t, dt, f, varargin{:}); %#ok<ASGLU>
        J0    = Df(x,        t,          varargin{:}) - gamma * DH(x,        t,          varargin{:});
        Jhalf = Df(x_half,   t+0.5*dt,   varargin{:}) - gamma * DH(x_half,   t+0.5*dt,   varargin{:});
        J1    = Df(x_full,   t+dt,       varargin{:}) - gamma * DH(x_full,   t+dt,       varargin{:});

        k1 = J0    * W;
        k2 = Jhalf * (W + 0.5*dt*k1);
        k3 = Jhalf * (W + 0.5*dt*k2);
        k4 = J1    * (W + dt*k3);
        W  = W + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

        g = norm(W);
        W = W / g;

        sum_log = sum_log + log(g);
    end

    mLCE = sum_log / (n_compute * dt);
end

function LCE = msf_LCE(x0, f, Df, DH, dt, T_init, T_cal, gamma, varargin)


    x = x0(:);
    n = numel(x);
    t = 0.0;

   
    Q = orth(randn(n));  
    acc = zeros(n,1);

    n_forward = max(1, round(T_init / dt));
    n_compute = max(1, round(T_cal  / dt));

   
    for k = 1:n_forward
        
        [x_full, x_half] = rk4_step_with_mid(x, t, dt, f, varargin{:});

        
        J0    = Df(x,       t,          varargin{:}) - gamma * DH(x,       t,          varargin{:});
        Jhalf = Df(x_half,  t+0.5*dt,   varargin{:}) - gamma * DH(x_half,  t+0.5*dt,   varargin{:});
        J1    = Df(x_full,  t+dt,       varargin{:}) - gamma * DH(x_full,  t+dt,       varargin{:});

        
        k1 = J0    * Q;
        k2 = Jhalf * (Q + 0.5*dt*k1);
        k3 = Jhalf * (Q + 0.5*dt*k2);
        k4 = J1    * (Q + dt*k3);
        Q  = Q + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

        
        [Q, ~] = qr(Q, 0);

       
        x = x_full;  t = t + dt;
    end

  
    for k = 1:n_compute
        
        [x_full, x_half] = rk4_step_with_mid(x, t, dt, f, varargin{:});

        
        J0    = Df(x,       t,          varargin{:}) - gamma * DH(x,       t,          varargin{:});
        Jhalf = Df(x_half,  t+0.5*dt,   varargin{:}) - gamma * DH(x_half,  t+0.5*dt,   varargin{:});
        J1    = Df(x_full,  t+dt,       varargin{:}) - gamma * DH(x_full,  t+dt,       varargin{:});

        
        k1 = J0    * Q;
        k2 = Jhalf * (Q + 0.5*dt*k1);
        k3 = Jhalf * (Q + 0.5*dt*k2);
        k4 = J1    * (Q + dt*k3);
        Q  = Q + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

        
        [Q, R] = qr(Q, 0);
        acc = acc + log(abs(abs(diag(R)) + realmin));

        
        x = x_full;  t = t + dt;
    end

    LCE = acc / (n_compute * dt);
end


function mLCE_list = msf_scan(gamma_list, x0, f, Df, DH, dt, T_init, T_cal, varargin)

    m = numel(gamma_list);
    mLCE_list = zeros(m,1);
    for i = 1:m 
        gamma = gamma_list(i);
        mLCE_list(i) = msf_mLCE(x0, f, Df, DH, dt, T_init, T_cal, gamma, varargin{:});
    end
end

function [x_full, x_half] = rk4_step_with_mid(x, t, dt, f, varargin)


    
    x = x(:);

    % --- k1 ---
    k1 = f(x, t, varargin{:}); 
    k1 = k1(:);
    assert(numel(k1)==numel(x), 'f(x,t,...) must return a vector with length numel(x).');

    x1 = x + 0.5*dt*k1;

    % --- k2 ---
    k2 = f(x1, t+0.5*dt, varargin{:}); 
    k2 = k2(:);
    assert(numel(k2)==numel(x));

    x2 = x + 0.5*dt*k2;

    % --- k3 ---
    k3 = f(x2, t+0.5*dt, varargin{:}); 
    k3 = k3(:);
    assert(numel(k3)==numel(x));

    x3 = x + dt*k3;

    % --- k4 ---
    k4 = f(x3, t+dt, varargin{:}); 
    k4 = k4(:);
    assert(numel(k4)==numel(x));

  
    x_full = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x_half = x + 0.5*dt*k2;  
end



end