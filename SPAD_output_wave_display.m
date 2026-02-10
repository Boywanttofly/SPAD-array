clearvars -except st_Vpp
clc
ch1=load('F:\USTC_reasearch_topic\UOWC\Exp_data\20250525\optical-noise_001.mat');
t=ch1.time;
ch1=ch1.data;
SI = 8e-10;
for i=1:size(ch1,1)
    if(i<3 || i>(size(ch1,1)-2))
        ch1_filter(i)=ch1(i);
    else
        ch1_filter(i)=mean(ch1(i-2:i+2));
    end
end
ch1_noise = ch1_filter(ch1_filter < 2e-3);
st_Vpp=importdata("st_vpp_20250525_n.mat");%需通过exp_ook_data_process.m的光子计数判决模块统计st_Vpp
st=40;
ta=159;%159为基准
x(1:st+ta+1)=0:8e-10:(st+ta)*8*1e-10;
y(1:st+ta+1)=0;
vpp=0.019327;
% figure();
cnt=0;
v0=0;
% noise_mean=3.52e-04;
% noise_var=(6.674e-04/sqrt(2))^2;
for i=1:size(st_Vpp,2)
    if(abs(ch1_filter(st_Vpp(i))-vpp)<1e-3)
        err(1:st+ta+1)=0;
        for j=1:st+ta+1
            %希望提取什么样的波形，就选一个参考波形作为基准，选择与其误差小于一定阈值的波形
            err(j)=(ch1_filter(st_Vpp(i)-st+j-1)-ch1_filter(st_Vpp(11)-st+j-1))^2;
        end
        rmse=mean(err);
        if(max(err)<1e-5)
            cnt=cnt+1;
            y=y+ch1_filter(st_Vpp(i)-st:st_Vpp(i)+ta);
            y_sum(cnt)=sum(ch1_filter(st_Vpp(i)-st:st_Vpp(i)+ta));
            vpp_st(cnt) = max(ch1_filter(st_Vpp(i)-st:st_Vpp(i)+ta));
            % plot(x,ch1_filter(st_Vpp(i)-st:st_Vpp(i)+ta));
            % hold on
        end
    end
end
y=y/cnt;
% v0 = max(y);
v0 = 0.019327;
% s_std=sum(y*SI)-noise_mean*size(y,2)*SI;
% v0=0.019327;%0.0196
% 确定上升沿参数
tr=6.5e-09;%5.139e-9
% 确定下降沿参数
tf=2.5009e-08;%28.45e-9
%删去脉冲前部的水平线，使起始位置函数值为0以便于拟合
x_fit = x(20:end);
y_fit = y(20:end);
% x_fit = x_fit - x_fit(1);
% y_fit = y_fit - y_fit(1);
x_fit_r = x_fit(1:21);%上升沿部分
y_fit_r = y_fit(1:21);
x_fit_f = x_fit(22:162);%下降沿部分
y_fit_f = y_fit(22:162);

%下降沿参数拟合：y取log用直线拟合,斜率得到tf，截距得到vpp
% y_fit_f = log(y_fit_f/vpp);
% x_fit_f = x_fit_f - x_fit_f(1);
% [fitresult, gof] = parafit(x_fit_f, y_fit_f);
% tf = 1/fitresult.p1;
% % 
% % %上升沿参数拟合：
% x_fit_r = x_fit(1:21)-x_fit(1);%上升沿部分
% y_fit_r = y_fit(1:21)-y_fit(1);%上升沿部分
% x_fit_r = x_fit_r(18:21);
% y_fit_r = y_fit_r(18:21);
% y_fit_r = abs(log(1-y_fit_r/vpp));
% [fitresult, gof] = parafit(x_fit_r, y_fit_r);
% tr = 1/fitresult.p1;

figure();
x_fit = x_fit - x_fit(1);
plot(x_fit,y_fit,'Color',[0,0,0.8],'LineWidth',2);
hold on
peak = find(y_fit==max(y_fit));
% y_ana(1:size(x,2))=noise_mean+sqrt(noise_var)*randn(1,size(x,2));
for j=1:size(x_fit,2)
    if(j<=peak)
        y_ana(j)=v0*(1-exp((-x_fit(j))/tr));
        %三次函数模型
        % y_ana(j)=x_fit(j)*(532727/64*(x_fit(j)/1e-9-16)^2+1076010);
    elseif(j>=peak)
        y_ana(j)=v0*exp((-x_fit(j)+peak*SI)/tf);
    end
end
% y_ana = [zeros(1,19) y_ana];
% y_ana = [y_ana  zeros(1,60)];
plot(x_fit,y_ana);
xlabel('t /s');
ylabel('V');
legend('Exp','Model');
%% *************************反函数导数计算**********************************
%% 基于上升沿导数
%由于dy/dx=1/(dx/dy),有dU^(-1)(v)/d(v) = 1/ (dU(t)/dt)|U=v,U为原始雪崩波形
%首先计算原始波形导数
% figure;
% %上升沿导数
% y_dr = abs((y_fit_r(2:end)-y_fit_r(1:end-1))./(x_fit_r(2:end)-x_fit_r(1:end-1)));%先求导数
% y_dr = 1 ./ y_dr;%再求倒数得到反函数导数
% %下降沿导数,由于下降沿导数较平缓，相邻间隔导数区分不明显，需拉大间隔
% y_df = abs((y_fit_f(1:end-2)-y_fit_f(3:end))./(x_fit_f(1:end-2)-x_fit_f(3:end)));%先求导数
% y_df = 1 ./ y_df;%再求倒数得到反函数导数
% % plot(x_der,y_der,'Color',[0,0,0.8])
% hold on
% % 将电压值与导数对应
% v_r = y_fit_r(1:end-1);%上升沿的电压值
% v_f = y_fit_f(2:end-1);%下降沿的电压值
% %将相同电压值的概率密度相加
% for v_i = 1 : size(v_r,2)
%     [~,loc] = min(abs(v_r(v_i)-v_f));
%     d_t(v_i) = y_dr(v_i) + y_df(loc);%d_t为合并概率密度
% end
% plot(v_r,d_t,'Color',[0.8,0,0]);
% save d_rt d_t
% save v_r v_r

%% 基于下降沿导数
%由于dy/dx=1/(dx/dy),有dU^(-1)(v)/d(v) = 1/ (dU(t)/dt)|U=v,U为原始雪崩波形
%首先计算原始波形导数
figure;
% 求整体导数
y_d = (y_fit(3:end)-y_fit(1:end-2))./(x_fit(3:end)-x_fit(1:end-2));
x_d = x_fit(2:end-1);
plot(x_d,y_d);
loc = find(y_fit==max(y_fit));%y_fit最大值处是区分上升沿和下降沿的点
%上升沿导数
%间隔2的取法
% y_dr = y_d(y_d>0);%上升沿与下降沿的分界线就是y_d由正变负的位置
% y_dr = 1 ./ y_dr;
% v_r = y_fit(find(y_d>0)+1);%y_d表示y_fit从第二个点开始的导数
% y_dr = y_dr(1:end-1);
% v_r = v_r(1:end-1);
%间隔1的取法
y_dr = abs((y_fit_r(2:end)-y_fit_r(1:end-1))./(x_fit_r(2:end)-x_fit_r(1:end-1)));%先求导数
y_dr = 1 ./ y_dr;%再求倒数得到反函数导数
v_r = y_fit_r(1:end-1);
%下降沿导数
%间隔2采样点的取法
% y_df = y_d(y_d<0);
% y_df = abs(1 ./ y_df);
% v_f = y_fit(find(y_d<0)+1);
%间隔1采样点的取法
y_df = abs((y_fit_f(1:end-2)-y_fit_f(3:end))./(x_fit_f(1:end-2)-x_fit_f(3:end)));%先求导数
y_df = 1 ./ y_df;%再求倒数得到反函数导数
v_f = y_fit_f(2:end-1);
% % plot(x_der,y_der,'Color',[0,0,0.8])
% hold on
% % 将电压值与导数对应
% v_r = y_fit_r(1:end-1);%上升沿的电压值
% v_f = y_fit_f(2:end);%下降沿的电压值
%将相同电压值的概率密度相加
for v_i = 1 : size(v_f,2)
    [~,loc] = min(abs(v_r-v_f(v_i)));
    d_t(v_i) = y_df(v_i) + y_dr(loc);%d_t为合并概率密度
end
locdt = figure;
plot(v_f,d_t,'LineWidth',2,'Color',[0.8,0,0]);
xlabel('Voltage (V)');
ylabel('p_V(v)');
plotformat(locdt)
save d_ft d_t
save v_f v_f