%%
clear variables

global A Ab B Bb C Cb Kb Lob Hb Bfb Cfb Pb Nab K D L
global N % number of followers
global n m r
global Mb Lvsb T U
global DYN DYN_AUG

%%
% System matrix
A = [-0.05 1 0 0; 
     -1 0 0 0;
     0 0 0 3 ; 
     0 0 -3 0];
B = [1   0;
     1   0; 
     0   1;
     -1   -1];
C = [1 0 0 0;
     0 0 1 0];
% 该设定下，第二个actuator可完全失效，sensor不可完全失效。

N = 5;
n = 4;
m = 2;
r = 2;

Ab = kron(eye(N),A);
Bb = kron(eye(N),B);
Cb = kron(eye(N),C);

% controlable and observable
rank(ctrb(A,B))
rank(obsv(A,C))

% fault for test
Bf = B*diag([1 1]);
Cf = diag([1 1])*C;
% controlable and observable
rank(ctrb(A,Bf))
rank(obsv(A,Cf))

% 
Bfb = kron(eye(N), Bf);
Cfb = kron(eye(N), Cf);

% topology
D = diag([1 0 0 0 0]);
W = [
0 1 1 0 1;
1 0 1 1 0;
1 1 0 1 0;
0 1 1 0 1;
1 0 0 1 0
];
L = -W + diag(sum(W, 2));
H = L+D;
delta_min = min(eig(H));
delta_max = max(eig(H));
Hb = kron(H, eye(n));

%%
Q = sdpvar(n);
M1 = [A*Q+Q*A'-2*delta_min*(B*B') Q; Q' -1/delta_min*eye(n)];
constraint = [M1<=0*eye(2*n)   , Q >= eye(n)];
options = sdpsettings('solver','sdpt3','verbose',0);
optimize(constraint,[],options);
check(constraint)
K = B'*inv(value(Q));

Lo = -place(A',C',[-2,-3,-2.5,-4])';
eig(A+Lo*C)

% validate
max(eig([kron(eye(N),A)-kron(H,B*K) -kron(H,B*K); zeros(N*n) kron(eye(N),A+Lo*C)]))

% 得到 K 和 Lo
Kb = kron(eye(N), K);
Lob = kron(eye(N), Lo);

%%
Nab = zeros(m*N);
Pb = zeros(r*N);
for i = 1:N
    Bfi = Bfb(n*i-n+1:n*i, m*i-m+1:m*i);
    Nai = pinv(Bfi)*B;
    Nab(m*i-m+1:m*i, m*i-m+1:m*i) = Nai;

    Cfi = Cfb(r*i-r+1:r*i, n*i-n+1:n*i);
    Pi = C*pinv(Cfi);
    Pb(r*i-r+1:r*i, r*i-r+1:r*i) = Pi;
end

[flag, Mb, Lvsb, T, U] = solve_hete_fault(Bfb, Cfb);
% [flag, Mb, Lvsb] = solve_hete_fault_no_etm(Bfb, Cfb);
if flag == 1
    disp("Solved!")
end

DYN_AUG = get_fault_dynamic();
% DYN = get_fault_dynamic_no_etm();
% max(real(eig(DYN)))
%%


%% Simulation

rng(0)
% 运算向量： deltaf, e, s, xdelta, xbar, eps(trigger)
deltaf_initial = rand(N*n, 1);
e_initial = rand(N*n, 1);
% s_initial = rand(N*n, 1);
s_initial = zeros(N*n, 1);
% xdelta_initial = rand(N*n, 1);
xdelta_initial = zeros(N*n, 1);
eps_initial = zeros(N*n,1);
x_initial = [deltaf_initial; e_initial; s_initial; xdelta_initial; eps_initial];

% trigger times 
times_trigger = zeros(N,1);
instant_trigger_history = cell(5,1);

tstart = 0;
tfinal = 50;
steps = 50000;
t_history = linspace(tstart, tfinal, steps);
dt = t_history(2)-t_history(1);
x_history = zeros(5*N*n, steps);

x = x_initial;
x_history(:,1) = x;

fault_config_act_initial = ones(N,2); % 每行是每个agent的actuator的衰减系数 1 为正常 0为丢失
fault_config_sens_initial = ones(N,2); % 每行是每个agent的sensor的衰减系数 1正常 0丢失

fault_config_act_curr = fault_config_act_initial;
fault_config_sens_curr = fault_config_sens_initial;
fault_config_act_last = fault_config_act_curr;  % 记录上一时刻的fault设置
fault_config_sens_last = fault_config_sens_curr; % 记录上一时刻的fault设置

eta_initial = 100; % detm auxiliary variable
lam = 1;
the = 1;
eta = eta_initial*ones(N,1);
eta_history = zeros(N, steps);
eta_history(:,1) = eta;

for is = 2:steps

    t = t_history(is);

    if (t>5)
        % 设置一个fault，比如 agent 2 的actuator fault
        fault_config_act_curr = [1 1; 0.5 0; 1 1; 1 1; 1 1];

    end

    if (t > 10)
        fault_config_act_curr = [1 1; 0.5 0; 1 1; 1 0.2; 1 1];
        fault_config_sens_curr = [0.3 1; 1 1; 1 1; 1 1; 1 1];
    end

    if (t > 20)
        fault_config_act_curr = [1 1; 1 1; 1 1; 1 1; 1 1];
        fault_config_sens_curr = [1 1; 1 1; 1 1; 0.8 0.3; 0.8 1];
    end

    if (t > 35)
        fault_config_act_curr = [0.5 0; 1 1; 1 1; 1 1; 1 0];
        fault_config_sens_curr = [1 1; 1 1; 1 1; 0.8 0.3; 0.8 1];
    end

    % if (t > 35)
    %     fault_config_act_curr = [1 1; 1 1; 1 1; 1 1; 1 1];
    %     fault_config_sens_curr = [1 1; 1 1; 1 1; 1 1; 1 1];
    % end





    if ( sum(fault_config_act_curr ~= fault_config_act_last,'all')  || sum(fault_config_sens_curr ~= fault_config_sens_last,'all') )
        disp('Fault detected')
        % 构造Bfb,Cfb矩阵
        [Bfb, Cfb] = construct_fault_matrix(fault_config_act_curr, fault_config_sens_curr, B, C);
        % 计算新的VA，VS参数，新的event参数
        [flag, Mb, Lvsb, T, U] = solve_hete_fault(Bfb, Cfb);
        % 计算新的仿真用DYN矩阵
        DYN_AUG = get_fault_dynamic();
        % 记录fault设置
        fault_config_act_last = fault_config_act_curr;  % 记录上一时刻的fault设置
        fault_config_sens_last = fault_config_sens_curr; % 记录上一时刻的fault设置

        % 设置初值
        if (t>5)
            % fault_config_act_curr = [1 1; 0.5 0; 1 1; 1 1; 1 1];
            x(3*N*n+1 + 1*n : 3*N*n+1 + 1*n + n-1) = rand(n,1);  % xdelta
        end

        if (t > 10)
            % fault_config_act_curr = [1 1; 0.5 0; 1 1; 1 0.2; 1 1];
            % fault_config_sens_curr = [0.3 1; 1 1; 1 1; 1 1; 1 1];
            x(3*N*n+1 + 3*n : 3*N*n+1 + 3*n + n-1) = rand(n,1);  % xdelta
            x(2*N*n+1 + 0*n : 2*N*n+1 + 0*n + n-1) = rand(n,1);  % s
        end

        if (t > 20)
            % fault_config_act_curr = [1 1; 1 1; 1 1; 1 1; 1 1];
            % fault_config_sens_curr = [1 1; 1 1; 1 1; 0.8 0.3; 0.8 1];
            x(3*N*n+1: 4*N*n ) = zeros(n*N,1);  % xdelta
            x(2*N*n+1 + 0*n : 2*N*n+1 + 0*n + n-1) = zeros(n,1);  % s
            x(2*N*n+1 + 3*n : 2*N*n+1 + 3*n + n-1) = rand(n,1);  % s
            x(2*N*n+1 + 4*n : 2*N*n+1 + 4*n + n-1) = rand(n,1);  % s
        end

        if (t > 25)
            % fault_config_act_curr = [0.5 0; 1 1; 1 1; 1 1; 1 0];
            % fault_config_sens_curr = [1 1; 1 1; 1 1; 0.8 0.3; 0.8 1];

            x(3*N*n+1 + 0*n : 3*N*n+1 + 0*n + n-1) = rand(n,1);  % xdelta
            x(3*N*n+1 + 4*n : 3*N*n+1 + 4*n + n-1) = rand(n,1);  % xdelta
        end

        % if (t > 35)
        %     % fault_config_act_curr = [1 1; 1 1; 1 1; 1 1; 1 1];
        %     % fault_config_sens_curr = [1 1; 1 1; 1 1; 1 1; 1 1];
        % 
        %     x(3*N*n+1: 4*N*n ) = zeros(n*N,1);  % xdelta
        %     x(2*N*n+1: 3*N*n ) = zeros(n*N,1);  % xdelta
        % end

    end
    
    % 后续可以加一个estimatefault需要的时间

    % detm
    deltaf = x(1:N*n);
    e = x(N*n+1:2*N*n);
    eps = x(4*N*n+1:5*N*n);
    z = -kron(H,eye(n))*(deltaf+e+eps);
    detadt = -lam*eta-eps'*kron(eye(N),U)*eps+z'*kron(eye(N),T)*z;
    eta = eta + detadt*dt;

    % Dynamic
    dxdt = DYN_AUG*x;
    x = x + dxdt*dt;

    deltaf = x(1:N*n);
    e = x(N*n+1:2*N*n);
    s = x(2*N*n+1:3*N*n);
    xdelta = x(3*N*n+1:4*N*n);
    eps = x(4*N*n+1:5*N*n);

    z = -kron(H,eye(n))*(deltaf+e+eps);

    % event
    for ia = 1:N
        eps_a = eps(ia*n-n+1:ia*n);
        z_a = z(ia*n-n+1:ia*n);
        eta_a = eta(ia);
        % theta_a = eps_a'*U*eps_a - z_a'*T*z_a;
        theta_a = the*(eps_a'*U*eps_a - z_a'*T*z_a)-eta_a;
        % theta_a = 1;

        if theta_a >= 0 % measurement error 归零
            eps(ia*n-n+1:ia*n) = zeros(n,1);
            times_trigger(ia) = times_trigger(ia)+1;
            instant_trigger_history{ia} = [instant_trigger_history{ia} t];
        end

    end


    x(4*N*n+1:5*N*n) = eps;

    % record
    x_history(:, is) = x;
    eta_history(:,is) = eta;




end

deltaf_history = x_history(1:N*n,:);
e_history = x_history(N*n+1:2*N*n, :);
s_history = x_history(2*N*n+1:3*N*n, :);
xdelta_history = x_history(3*N*n+1:4*N*n, :);
eps_history = x_history(4*N*n+1:5*N*n, :);

fprintf("Average IET: %f s\n", tfinal/(sum(times_trigger))*N);

for i = 1:N
    instant_trigger_i = instant_trigger_history{i};
    iet_i = instant_trigger_i(2:end) - instant_trigger_i(1:end-1);
    fprintf("Agent %d, min IET: %f s, max IET: %f s, Average: %f s \n", ...
        i, min(iet_i), max(iet_i), mean(iet_i));
end


%% better plot
[all_themes, all_colors] = GetColors();
% 检查，线宽，颜色，字体
%% Consensus error
figure() 
margin_up = 0.02; % 总图上边缘
margin_down = 0.12; % 总图下边缘
margin_left = 0.1; % 总图左边缘
margin_right = 0.03; % 总图右边缘
gap = 0.03; % 子图间隔
ywidth = (1-margin_up-margin_down-gap*(N-1))/N;

for i = 1:N
    subplot(N,1,i)
    plot(t_history, deltaf_history(i*n-n+1:i*n,:), 'LineWidth',1)
    set(gca,'xlim',[tstart tfinal]);
    set(gca, 'colororder', all_themes{1}(2:end,:));
    if i ~= N
        set(gca,'xticklabel',[]);
    end
    margin_sub_down = 1-margin_up-i*ywidth-(i-1)*gap;
    set(gca,'position',[margin_left margin_sub_down 1-margin_right-margin_left ywidth])
    set(gca, 'Fontname', 'Times New Roman', 'Fontsize',11)
    ylabel('$\delta^f_{'+string(i)+'}$','interpreter','latex')
    grid on
end
xlabel("Time (s)")

%% Observer error
figure() 
margin_up = 0.02; % 总图上边缘
margin_down = 0.12; % 总图下边缘
margin_left = 0.1; % 总图左边缘
margin_right = 0.03; % 总图右边缘
gap = 0.03; % 子图间隔
ywidth = (1-margin_up-margin_down-gap*(N-1))/N;
for i = 1:N
    subplot(N,1,i)
    plot(t_history, e_history(i*n-n+1:i*n,:))
    set(gca,'xlim',[tstart tfinal]);
    if i ~= N
        set(gca,'xticklabel',[]);
    end
    margin_sub_down = 1-margin_up-i*ywidth-(i-1)*gap;
    set(gca,'position',[margin_left margin_sub_down 1-margin_right-margin_left ywidth])
    set(gca, 'Fontname', 'Times New Roman', 'Fontsize',11)
    ylabel('$e_{'+string(i)+'}$','interpreter','latex')
    grid on
end
xlabel("Time (s)")

%% Virtual actuator error
figure() 
margin_up = 0.02; % 总图上边缘
margin_down = 0.12; % 总图下边缘
margin_left = 0.1; % 总图左边缘
margin_right = 0.03; % 总图右边缘
gap = 0.03; % 子图间隔
ywidth = (1-margin_up-margin_down-gap*(N-1))/N;
for i = 1:N
    subplot(N,1,i)
    plot(t_history, xdelta_history(i*n-n+1:i*n,:),'LineWidth',1)
    set(gca,'xlim',[tstart tfinal]);
    set(gca, 'colororder', all_themes{1}(2:end,:));
    if i ~= N
        set(gca,'xticklabel',[]);
    end
    if i == 3
        set(gca,'ylim',[-1 1]);
    end
    margin_sub_down = 1-margin_up-i*ywidth-(i-1)*gap;
    set(gca,'position',[margin_left margin_sub_down 1-margin_right-margin_left ywidth])
    set(gca, 'Fontname', 'Times New Roman', 'Fontsize',11)
    ylabel('$x^{\Delta}_{'+string(i)+'}$','interpreter','latex')
    grid on
end
xlabel("Time (s)")

%% Virtual Sensor error
figure() 
margin_up = 0.02; % 总图上边缘
margin_down = 0.12; % 总图下边缘
margin_left = 0.1; % 总图左边缘
margin_right = 0.03; % 总图右边缘
gap = 0.03; % 子图间隔
ywidth = (1-margin_up-margin_down-gap*(N-1))/N;
for i = 1:N
    subplot(N,1,i)
    plot(t_history, s_history(i*n-n+1:i*n,:),'LineWidth',1)
    set(gca, 'colororder', all_themes{1}(2:end,:));
    set(gca,'xlim',[tstart tfinal]);
    if i ~= N
        set(gca,'xticklabel',[]);
    end
    margin_sub_down = 1-margin_up-i*ywidth-(i-1)*gap;
    set(gca,'position',[margin_left margin_sub_down 1-margin_right-margin_left ywidth])
    set(gca, 'Fontname', 'Times New Roman', 'Fontsize',11)
    ylabel('$s_{'+string(i)+'}$','interpreter','latex')
    grid on
end
xlabel("Time (s)")

%% VA VS together
figure()
tiledlayout(5, 2, 'TileSpacing', 'tight','Padding','tight');
for i = 1 : N
    nexttile;
    
    % VA
    plot(t_history, xdelta_history(i*n-n+1:i*n,:),'LineWidth',0.57)
    ylabel('$x^{\Delta}_{'+string(i)+'}$','interpreter','latex')
    grid on
    if i ~= N
        set(gca,'xticklabel',[]);
    end
    set(gca, 'Fontname', 'Times New Roman', 'Fontsize',11)
    set(gca,'xlim',[tstart tfinal]);
    set(gca, 'colororder', all_themes{1}(2:end,:));
    if i == 3
        set(gca,'ylim',[-1 1]);
    end
    if i == N
        xlabel("Time (s)")
    end

    nexttile;
    % VS
    plot(t_history, s_history(i*n-n+1:i*n,:),'LineWidth',0.7)
    ylabel('$s_{'+string(i)+'}$','interpreter','latex')
    grid on
    if i ~= N
        set(gca,'xticklabel',[]);
    end
    set(gca, 'Fontname', 'Times New Roman', 'Fontsize',11)
    set(gca,'xlim',[tstart tfinal]);
    set(gca, 'colororder', all_themes{2}(2:end,:));
    if i == N
        xlabel("Time (s)")
    end

end





%% Event error
figure() 
margin_up = 0.02; % 总图上边缘
margin_down = 0.12; % 总图下边缘
margin_left = 0.1; % 总图左边缘
margin_right = 0.03; % 总图右边缘
gap = 0.03; % 子图间隔
ywidth = (1-margin_up-margin_down-gap*(N-1))/N;


for i = 1:N
    subplot(N,1,i)
    plot(t_history, eps_history(i*n-n+1:i*n,:))
    set(gca,'xlim',[tstart tfinal]);
    if i ~= N
        set(gca,'xticklabel',[]);
    end
    margin_sub_down = 1-margin_up-i*ywidth-(i-1)*gap;
    set(gca,'position',[margin_left margin_sub_down 1-margin_right-margin_left ywidth])
    set(gca, 'Fontname', 'Times New Roman', 'Fontsize',11)
    ylabel('$\varepsilon_{'+string(i)+'}$','interpreter','latex')
    grid on
end
xlabel("Time (s)")

%% Event error with norm
figure() 
margin_up = 0.02; % 总图上边缘
margin_down = 0.12; % 总图下边缘
margin_left = 0.1; % 总图左边缘
margin_right = 0.03; % 总图右边缘
gap = 0.03; % 子图间隔
ywidth = (1-margin_up-margin_down-gap*(N-1))/N;
for i = 1:N
    subplot(N,1,i)
    plot(t_history, vecnorm(eps_history(i*n-n+1:i*n,:)),'Color','k','LineWidth',0.7)
    set(gca,'xlim',[tstart tfinal]);
    if i ~= N
        set(gca,'xticklabel',[]);
    end
    margin_sub_down = 1-margin_up-i*ywidth-(i-1)*gap;
    set(gca,'position',[margin_left margin_sub_down 1-margin_right-margin_left ywidth])
    set(gca, 'Fontname', 'Times New Roman', 'Fontsize',11)
    ylabel('$\| \varepsilon_{'+string(i)+'} \|$','interpreter','latex')
    grid on
end
xlabel("Time (s)")



%% ausxiliary variable
figure() 
margin_up = 0.02; % 总图上边缘
margin_down = 0.12; % 总图下边缘
margin_left = 0.1; % 总图左边缘
margin_right = 0.03; % 总图右边缘
gap = 0.03; % 子图间隔
ywidth = (1-margin_up-margin_down-gap*(N-1))/N;
for i = 1:N
    subplot(N,1,i)
    plot(t_history, eta_history(i,:))
    set(gca,'xlim',[tstart tfinal]);
    if i ~= N
        set(gca,'xticklabel',[]);
    end
    margin_sub_down = 1-margin_up-i*ywidth-(i-1)*gap;
    set(gca,'position',[margin_left margin_sub_down 1-margin_right-margin_left ywidth])
    set(gca, 'Fontname', 'Times New Roman', 'Fontsize',11)
    ylabel('$\eta_{'+string(i)+'}$','interpreter','latex')
    grid on
end
xlabel("Time (s)")

%% IET
figure() 
margin_up = 0.02; % 总图上边缘
margin_down = 0.12; % 总图下边缘
margin_left = 0.1; % 总图左边缘
margin_right = 0.03; % 总图右边缘
gap = 0.03; % 子图间隔
ywidth = (1-margin_up-margin_down-gap*(N-1))/N;
for i = 1:N
    subplot(N,1,i)

    instant_trigger_i = instant_trigger_history{i};
    iet_i = instant_trigger_i(2:end) - instant_trigger_i(1:end-1);
    plot(instant_trigger_i(2:end), iet_i,'LineWidth',1)

    set(gca,'xlim',[tstart tfinal]);
    if i ~= N
        set(gca,'xticklabel',[]);
    end
    margin_sub_down = 1-margin_up-i*ywidth-(i-1)*gap;
    set(gca,'position',[margin_left margin_sub_down 1-margin_right-margin_left ywidth])
    set(gca, 'Fontname', 'Times New Roman', 'Fontsize',11)
    ylabel('IET$_{'+string(i)+'}$','interpreter','latex')
    grid on
end
xlabel("Time (s)")

%% IET one figure
figure() 
tiledlayout(1, 1, 'TileSpacing', 'tight','Padding','tight');
nexttile;
for i = 1:N
    instant_trigger_i = instant_trigger_history{i};
    iet_i = instant_trigger_i(2:end) - instant_trigger_i(1:end-1);
    plot(instant_trigger_i(2:end), iet_i,'LineWidth',1)
    hold on
end
set(gca, 'colororder', all_themes{1}(1:end,:));
set(gca, 'Fontname', 'Times New Roman', 'Fontsize',11)
ylabel('IET','interpreter','latex')
grid on
set(gca,'xlim',[tstart tfinal]);
xlabel("Time (s)")
legend('agent 1', 'agent 2', 'agent 3', 'agent 4', 'agent 5')


%%




