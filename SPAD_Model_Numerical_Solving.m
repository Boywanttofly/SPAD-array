clear
%参数定义
N_array=3600; 
Tf=112e-9;%死时间-下降140e-9
Tr=16e-9;%死时间-上升20e-9
Td=Tr+Tf;%总死时间
tf=2.5009e-08;%SPAD输出波形下降时间常数 2.526e-08 2.668e-08
tr=6.5e-09;%SPAD输出波形上升时间常数 6.5e-09
Vs=0.019327;%根据20250525noise1获取，峰值减谷值 0.0192 0.019327
p_xt=0.1647;%0.1647
SR=1250e6;%采样速率
SI=1/SR;%采样间隔
ook_EC=1e-6;%qkd码元周期
k_ook=round(ook_EC/SI);
PC_arr=[19.2 36.7 83.4 192.4 605.0 1259.7];%平均光子数，单位：1e6
addpath('F:\USTC_reasearch_topic\UOWC\Long Distance');%createFit_for_exp_data_DC.m定义在该文件夹中

%ook发送数据
%% *************************基于指数函数的模型预测***************************
Order_number_ook=15;%m序列的阶数等于9
ook_data=(idinput((2^(Order_number_ook)-1),'prbs')'+1)/2;%生成m序列
ook_data(size(ook_data,2)+1)=1;

ch1=load(['F:\USTC_reasearch_topic\UOWC\Exp_data\20250525\optical-signal_000.mat']);%optical-signal','_00',num2str(en-1),'.mat
t=ch1.time;
ch1=ch1.data;
% figure();
for i=1:size(ch1,1)
    if(i<3 || i>(size(ch1,1)-2))
        ch1_filter(i)=ch1(i);
    else
        ch1_filter(i)=mean(ch1(i-2:i+2));
    end
end

% 1位置数据提取
% cnt = 0;
% for ook_i = 1 : size(ook_data,2)
%     if(ook_data(ook_i) == 1)
%         cnt = cnt + 1;
%         ch1_1((cnt-1) * k_ook + 1 : cnt * k_ook) = ch1_filter((ook_i-1) * k_ook + 1 : ook_i * k_ook);
%     end
% end
% xlabel('t/s');
% ylabel('Voltage (V)');
% 噪声参数获取
[BinHeight,BinCenter] = createFit_for_exp_data_DC(ch1_filter);
n_loc_max = find(BinHeight == max(BinHeight));
%认为BinHeight中位数左侧的为真实噪声分布，将左侧反转获取右侧分布
nfit_y = [BinHeight(1:n_loc_max) fliplr(BinHeight(1:n_loc_max-1))];
nfit_gap = BinCenter(2)-BinCenter(1);
nfit_x = BinCenter(1) : nfit_gap : BinCenter(1) + (size(nfit_y,2)-1)*nfit_gap;
[fitresult, gof] = noise_fit(nfit_x, nfit_y);
param_name = coeffnames(fitresult); % 获取参数名称
param_value = coeffvalues(fitresult); % 获取参数值
n_mean=param_value(2);%出现过-3.7E-4
n_var=param_value(3)^2/2;%底噪方差=c1^2/2
a=tr*Vs*(exp(-Tr/tr)-1+Tr/tr);
b=tf*Vs*(1-exp(-Tf/tf));
thr1 = 1e-15;%低概率忽略门限值
for en = 1:6%exp_number
    % 接收信号读取
    ch1=load(['F:\USTC_reasearch_topic\UOWC\Exp_data\20250525\optical-signal','_00',num2str(en-1),'.mat']);%optical-signal','_00',num2str(en-1),'.mat
    t=ch1.time;
    ch1=ch1.data;
    % figure();
    for i=1:size(ch1,1)
        if(i<3 || i>(size(ch1,1)-2))
            ch1_filter(i)=ch1(i);
        else
            ch1_filter(i)=mean(ch1(i-2:i+2));
        end
    end
    % 1位置数据提取
    % cnt = 0;
    % for ook_i = 1 : size(ook_data,2)
    %     if(ook_data(ook_i) == 1)
    %         cnt = cnt + 1;
    %         ch1_1((cnt-1) * k_ook + 1 : cnt * k_ook) = ch1_filter((ook_i-1) * k_ook + 1 : ook_i * k_ook);
    %     end
    % end
    lambda_est(en)=1/((a+b)/(mean(ch1_filter)-n_mean)/(1-p_xt)-Td/N_array);
    % xlabel('t/s');
    % ylabel('Voltage (V)');
    [BinHeight,BinCenter] = createFit_for_exp_data_DC(ch1_filter);
    clear V_range P_d P_V
    % axes('position',[0.55,0.55,0.3,0.3]);%关键在这句！所画的小图
% end

%积分运算
%20250525各组数据估计光强
%统计串扰概率的数据总计数率为2.1e6，信号（分布反推计数）：[18.91	36.05	81.49	187.47	588.52	1224.25]e6
%脉冲堆叠区计数：20.26	38.54	86.90	199.30	617.11	1260.41
%直接卷积法
Photon_count_arrive=lambda_est(en)/N_array;%光子数需通过模型均值-入射光强估计法获取/(1/(1-p_crosstalk))
Photon_count_arrive=PC_arr(en)*1e6/N_array;
%预先设置好输出电压为0的概率
P_V(1)=1/(Photon_count_arrive*Td+1);
%最大输出电压是像素数乘以单像素最大输出电压
v_step = BinCenter(2)-BinCenter(1);
V_range = 0 : v_step : max(ch1_filter);
j_up=N_array;%定义串扰数上界
for i=2:size(V_range,2)
    P_V(i)=0;
    v=V_range(i);
    %考虑串扰导致的幅度变化
    %如果v是单像素最高幅度的j_down倍，则v至少是由j_down上取整个串扰产生
    j_down=ceil(v/Vs);
    for j=j_down:j_up
        f(i) = 0;
        %上升沿的P(V)
        if(v <= j * Vs * (1 - exp(-Tr / tr)))
            f(i) = f(i) + (1-p_xt)*p_xt^(j-1)*P_V(1)*Photon_count_arrive*tr/(j*Vs-v);
            %sigmoid函数
            % f(i) = f(i) + (1-p_xt)*p_xt^(j-1)*P_V(1)*Photon_count_arrive*tr2*(1/v+1/(j*0.0207-v));
        end
        %下降沿的P(V)
        if((v >= j * Vs * exp(-Tf / tf))&&(v <= j * Vs))
            f(i) = f(i) + (1-p_xt)*p_xt^(j-1)*Photon_count_arrive*tf*P_V(1)/v;
        end
        %将极小的增量略去，可以大幅降低运算时间
        if(f(i) * v_step > thr1)
            P_V(i)=P_V(i)+f(i)*v_step;
        end
    end
end
%将极小的概率略去，可以大幅降低后续卷积运算时间
V_range = V_range(P_V > thr1);
P_V = P_V(P_V > thr1);
P_V = P_V/sum(P_V);%概率归一化
P_V_N=P_V;
%V_range_N记录筛选后电压范围
V_range_N = V_range;
for z=1:N_array-1
    P_V_N=conv(P_V_N,P_V);
    V_range_N = min(V_range_N) + min(V_range) : v_step : min(V_range_N) + min(V_range) + (size(P_V_N,2) - 1) * v_step;
    V_range_N = V_range_N(P_V_N > thr1);
    P_V_N = P_V_N(P_V_N > thr1);
    P_V_N = P_V_N / sum(P_V_N);%概率归一化
end
% V = 0 : v_step : (size(P_V_N,2)-1)*v_step;
% plot(V,P_V_N / v_step);
% hold on
%噪声分布计算
%在噪声为主导的波形中，找到概率最大值以获取噪声分布中位数

%需要将噪声取值设置为BinCenter减去v_step的整数倍，以保证叠加信号后取值可以准确取到BinCenter的值
nlow_len = round((min(BinCenter) + 0.0105) / v_step);
nup_len = round((min(BinCenter) - 0.0105) / v_step);
V_noise = min(BinCenter) - nlow_len * v_step : v_step : min(BinCenter) - nup_len * v_step;

P_V_noise = 1/sqrt(2*pi*n_var).*exp(-(V_noise-n_mean).^2/2/n_var)*v_step;

P_V_out = conv(P_V_N,P_V_noise);
V = min(V_range_N) + min(V_noise) : v_step : (size(P_V_out,2) - 1) * v_step + min(V_range_N) + min(V_noise);
% 将概率极小的位置的去除
% 处理模型预测数据
% V = V(P_V_out > thr1);
P_V_out(P_V_out < thr1) = thr1;
% P_V_out = P_V_out / sum(P_V_out);%概率归一化
% 处理原始数据（可能出现0值，因此必须滤除）
% BinCenter = BinCenter(BinHeight > thr1);
%将自变量小数点后四位之后的误差消除,1e-4为实验所采集的电压值的常见量级
V = round(V/1e-4)*1e-4;
BinCenter = round(BinCenter/1e-4)*1e-4;

% BinHeight = BinHeight / sum(BinHeight);%概率归一化
% 计算两个分布之间的KL距离
% 由于两个分布自变量非严格对齐，需确定二者交集
% 确定下界:r开头代表实际分布，a开头代表模型预测分布

%addBH、addV必须初始化，否则变长后即使下一次实际长度变短也依然保留长长度
addBH = [];
addV = [];
if(min(V) < min(BinCenter))
    %区间最小值为V最小值
    rlow = find(BinCenter == min(BinCenter));
    [~,alow] = min(abs(V - min(BinCenter)));
    addlen = round((min(BinCenter) - min(V)) / v_step);%计算从V到BinCenter需要多少步
    addBH(1 : addlen) = 0;
    BinCenter = [min(V) : v_step : min(V) + (addlen-1) * v_step BinCenter];
    BinHeight = [addBH BinHeight];
else
    %区间最小值为BinCenter最小值
    [~,rlow] = min(abs(BinCenter - min(V)));
    alow = find(V == min(V));
    addlen = round((min(V) - min(BinCenter)) / v_step);
    addV(1 : addlen) = 0;
    V = [min(BinCenter) : v_step : min(BinCenter) + (addlen-1) * v_step V];
    P_V_out = [addV P_V_out];
end
% 确定上界
addBH = [];
addV = [];
if(max(V) > max(BinCenter))
    %区间最大值为V最大值
    rup = find(BinCenter == max(BinCenter));
    [~,aup] = min(abs(V - max(BinCenter)));
    addlen = round((max(V) - max(BinCenter)) / v_step);
    addBH(1 : addlen) = 0;
    BinCenter = [BinCenter(1 : end - 1) max(BinCenter) : v_step : max(BinCenter) + addlen * v_step];
    BinHeight = [BinHeight addBH];
else
    %区间最大值为BinHeight最大值
    [~,rup] = min(abs((BinCenter - max(V))));
    aup = find(V == max(V));
    addlen = round((max(BinCenter) - max(V)) / v_step);
    addV(1 : addlen) = 0;
    V = [V(1 : end - 1) max(V) : v_step : max(V) + addlen * v_step];
    P_V_out = [P_V_out addV];
end
BinHeight(BinHeight * v_step < thr1) = thr1 / v_step;
P_V_out(P_V_out < thr1) = thr1;%这是概率
p_V_out = P_V_out / v_step;%获取pdf
er_len = min([size(P_V_out,2) size(BinHeight,2)]);
% 储存分布均方差计算结果。如果需要计算不考虑串扰的分布差，需将p_xt置0再重新运行。
er_kl(en) = 0;
for kl_i = 1 : er_len
    er_kl(en) = er_kl(en) + BinHeight(kl_i) * v_step * log(BinHeight(kl_i) * v_step / P_V_out(kl_i));%是概率所以用大写((BinHeight(kl_i) - P_V_out(kl_i)/v_step))^2
end
% 计算分布残差
% pdf_resi_fit = p_V_out - BinHeight;
% er_re(en) = er_re(en) / er_len;
%绘制模型预测分布图
loc1 = figure;
plot(BinCenter(1:end),BinHeight(1:end),'LineWidth',3,'Color',[0,0,0.8]);
hold on
plot(V,p_V_out,'LineWidth',3,'Color',[0.8,0,0]);
% 绘制残差图
% loc2 = figure;
% plot(BinCenter(1:end),pdf_resi);
% hold on
% plotformat(loc2)

%% ************************基于原始雪崩波形的斜率计算************************
%提取原始波形的斜率倒数数据（form 工位电脑matlab处理）
Photon_count_arrive=PC_arr(en)*1e6/N_array;
V_range = 0 : v_step : max(ch1_filter);
v_r = importdata('v_f.mat');%电压值
d_t = importdata('d_ft.mat');%对应概率密度
% Photon_count_arrive = 20.26e6 / N_array;
% Photon_count_arrive = Photon_count_arrive / 1.1;
d_t = d_t*Photon_count_arrive/(Photon_count_arrive*Td+1);%实际概率密度是反函数导数*Photon_count_arrive/(Photon_count_arrive*Td+1)
% locdt = figure;
% plot(v_r,d_t,'LineWidth',2,'Color',[0.8,0,0]);
% xlabel('Voltage (V)');
% ylabel('p_V(v)');
% plotformat(locdt)
Vs = 0.019327;%10点滤波0.018984 五点滤波0.019327 三点滤波0.019397
%以下归一化操作确保pdf与对应v的积分之和为1
% step_r = v_r(1:end-1) - v_r(2:end);%概率密度需要乘以区间长度以获取概率，方便归一化，直接用概率密度不好归一化
% d_t = d_t(2:end);%保证区间与v_r一致
% d_t = d_t.*step_r;%此时d_t具有概率含义
% d_t = [P_V(1), d_t];%上述计算均为非0处概率，需补充0处概率。
% % %将极小的概率略去，可以大幅降低后续卷积运算时间
% % v_r = v_r(d_t > thr1);
% % d_t = d_t(d_t > thr1);
% d_t = d_t/sum(d_t);%概率归一化
% d_t(2:end) = d_t(2:end)./step_r;%得到调整后的概率密度值
% for i=1:size(d_t,2)
%     if(i<3 || i>(size(d_t,2)-2))
%         d_t_filter(i)=d_t(i);
%     else
%         d_t_filter(i)=mean(d_t(i-2:i+2));
%     end
% end
% v_r = v_r(3:end);
% d_t = d_t(3:end);
%预先设置好输出电压为0的概率
P_d(1)=1/(Photon_count_arrive*Td+1);%所有基于实际曲线导数的概率/密度用d表示
for i=2:size(V_range,2)
    P_d(i)=0;
    v=V_range(i);
    %考虑串扰导致的幅度变化
    %如果v是单像素最高幅度的j_down倍，则v至少是由j_down上取整个串扰产生
    j_down=ceil(v/Vs);
    for j=j_down:j_up
        %计算每一种串扰增量的概率前必须将f(i)清0
        delta(i) = 0;
        %寻找最接近v/(j)的v_r作为概率密度值
        [~,loc] = min(abs(v_r-v/j));
        %注意概率密度值要除以j
        delta(i) = (1-p_xt)*p_xt^(j-1)*d_t(loc) / j;
        %将极小的增量略去，可以大幅降低运算时间
        if(delta(i) * v_step > thr1)
            P_d(i) = P_d(i) + delta(i) * v_step;
        end
    end
end
%将极小的概率略去，可以大幅降低后续卷积运算时间
V_range = V_range(P_d > thr1);
P_d = P_d(P_d > thr1);
P_d = P_d/sum(P_d);%概率归一化
P_d_N=P_d;
%V_range_N记录筛选后电压范围
V_range_N = V_range;
for z=1:N_array-1
    P_d_N=conv(P_d_N,P_d);
    V_range_N = min(V_range_N) + min(V_range) : v_step : min(V_range_N) + min(V_range) + (size(P_d_N,2) - 1) * v_step;
    V_range_N = V_range_N(P_d_N > thr1);
    P_d_N = P_d_N(P_d_N > thr1);
    P_d_N = P_d_N / sum(P_d_N);%概率归一化
end
P_d_out = conv(P_d_N,P_V_noise);
V = min(V_range_N) + min(V_noise) : v_step : (size(P_d_out,2) - 1) * v_step + min(V_range_N) + min(V_noise);
P_d_out = P_d_out / v_step;%获取pdf
figure(loc1)
plot(V,P_d_out,'LineWidth',3,'Color',[0,0.8,0]);
legend('Exp','Ana1','Ana2');
xlabel('V_a (V)');
ylabel('PDF');
plotformat(loc1)
V = round(V/1e-4)*1e-4;
BinCenter = round(BinCenter/1e-4)*1e-4;

% BinHeight = BinHeight / sum(BinHeight);%概率归一化
% 计算两个分布之间的KL距离
% 由于两个分布自变量非严格对齐，需确定二者交集
% 确定下界:r开头代表实际分布，a开头代表模型预测分布

%addBH、addV必须初始化，否则变长后即使下一次实际长度变短也依然保留长长度
addBH = [];
addV = [];
if(min(V) < min(BinCenter))
    %区间最小值为V最小值
    rlow = find(BinCenter == min(BinCenter));
    [~,alow] = min(abs(V - min(BinCenter)));
    addlen = round((min(BinCenter) - min(V)) / v_step);%计算从V到BinCenter需要多少步
    addBH(1 : addlen) = 0;
    BinCenter = [min(V) : v_step : min(V) + (addlen-1) * v_step BinCenter];
    BinHeight = [addBH BinHeight];
else
    %区间最小值为BinCenter最小值
    [~,rlow] = min(abs(BinCenter - min(V)));
    alow = find(V == min(V));
    addlen = round((min(V) - min(BinCenter)) / v_step);
    addV(1 : addlen) = 0;
    V = [min(BinCenter) : v_step : min(BinCenter) + (addlen-1) * v_step V];
    P_d_out = [addV P_d_out];
end
% 确定上界
addBH = [];
addV = [];
if(max(V) > max(BinCenter))
    %区间最大值为V最大值
    rup = find(BinCenter == max(BinCenter));
    [~,aup] = min(abs(V - max(BinCenter)));
    addlen = round((max(V) - max(BinCenter)) / v_step);
    addBH(1 : addlen) = 0;
    BinCenter = [BinCenter(1 : end - 1) max(BinCenter) : v_step : max(BinCenter) + addlen * v_step];
    BinHeight = [BinHeight addBH];
else
    %区间最大值为BinHeight最大值
    [~,rup] = min(abs((BinCenter - max(V))));
    aup = find(V == max(V));
    addlen = round((max(BinCenter) - max(V)) / v_step);
    addV(1 : addlen) = 0;
    V = [V(1 : end - 1) max(V) : v_step : max(V) + addlen * v_step];
    P_d_out = [P_d_out addV];
end
BinHeight(BinHeight * v_step < thr1) = thr1 / v_step;
P_d_out(P_d_out * v_step < thr1) = thr1 / v_step;%这是概率密度
er_len = min([size(P_d_out,2) size(BinHeight,2)]);
% 储存分布均方差计算结果。如果需要计算不考虑串扰的分布差，需将p_xt置0再重新运行。
er_kl_wave(en) = 0;
for kl_i = 1 : er_len
    er_kl_wave(en) = er_kl_wave(en) + BinHeight(kl_i) * v_step * log(BinHeight(kl_i) * v_step / (P_d_out(kl_i) * v_step));%((BinHeight(kl_i) - P_V_out(kl_i)/v_step))^2
end
% 计算分布残差
pdf_resi_wave = P_d_out - BinHeight;
% 绘制残差图
% loc2 = figure;
% plot(BinCenter(1:end),pdf_resi);
% hold on
% plotformat(loc2)
end
floc = figure(7);
plot(PC_arr*1e6,er_kl);
hold on
plot(PC_arr*1e6,er_kl_wave);
legend('with-crosstalk-ana1','with-crosstalk-ana2','without-crosstalk-ana1','without-crosstalk-ana2');
xlabel('Average photon number /s');
ylabel('KL');
plotformat(floc)
% loc3 = figure(7);
% plot(lambda_est,er_kl_wave,'LineWidth',2)
% xlabel('Average photon number /s');
% ylabel('kl');
% hold on
% plotformat(loc3)
% 
% rmpath('D:\USTC-PhD\Reasearch\UOWC\Exp\matlab\ook');