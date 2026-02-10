clear;
%% **************************原始信号生成**************************
Order_number_ook=16;%m序列的阶数等于9
ook_data=(idinput((2^(Order_number_ook)-1),'prbs')'+1)/2;%生成m序列
ook_data(size(ook_data,2)+1)=1;

SR=125e7 / 1;%采样速率
SI=1/SR;%采样间隔
ook_EC=1e-6;%qkd码元周期
k_ook=round(ook_EC/SI);

Td_rise=16e-9;%上升时间
Td_fall=96e-9;%猝灭死时间
Td_all=Td_fall+Td_rise;
tao_rise=6.2e-9;%SPAD输出波形下降时间常数
tao_fall=25.3e-9;%SPAD输出波形下降时间常数
N_array=3600;
noise_mean=-0.00054596;%出现过-3.7E-4
noise_var=4.46085e-07;%底噪方差

%% *****************************单文件处理*********************************
% err_count=0;%误码统计
% % 接收信号读取
ch1=load('F:\USTC_reasearch_topic\UOWC\Exp_data\20250525\optical-noise_001.mat');
t=ch1.time;
ch1=ch1.data;

%%
% ch1=V_sample';
% ch1=V_sample;
figure();
plot(t,ch1);
hold on
xlabel('t/s');
ylabel('Voltage (V)');
% 均值滤波
for i=1:size(ch1,1)
    if(i<6 || i>(size(ch1,1)-5))
        ch1_filter(i)=ch1(i);
    else
        ch1_filter(i)=mean(ch1(i-5:i+5));
    end
end
plot(t,ch1_filter);
vol_thr=0.01;%脉冲门限值
pw_std=33e-9;%单脉冲半波全宽基本值
Vpp_std=0.0193;%单脉冲幅度基本值
gap_cnt=0;
ook_ENum=round((max(t)-min(t))/ook_EC);
k_ook=round(ook_EC/SI);
PC(1:ook_ENum)=0;
%%记录脉冲峰值
cnt=0;
for i=1:ook_ENum
    st=0;
    tail=0;
    for j=1:k_ook
        if(j<k_ook && ch1_filter((i-1)*k_ook+j) >= vol_thr)
            if(st==0)
                st=j;
                tail=j;
            else
                tail=j;
            end
        else
            if(st~=0)
                x=ch1_filter((i-1)*k_ook+st:(i-1)*k_ook+tail);
                %极值检测
                for z=2:size(x,2)-1
                    if((x(z)>=x(z-1))&&(x(z)>=x(z+1)))
                        PC(i)=PC(i)+round(x(z)/Vpp_std);
                        cnt=cnt+1;
                        vpp_original(cnt)=x(z);
                        st_Vpp(cnt)=round(x(z)/Vpp_std);
                    end
                end
                Vpp=max(ch1_filter((i-1)*k_ook+st:(i-1)*k_ook+tail));
                % t_dur(cnt)=(tail-st+1);
                % st_Vpp(cnt)=(i-1)*k_ook+st+find(x==max(x),1)-1;
                % Vpp_cnt(cnt)=round(Vpp/Vpp_std);
                %*round((tail-st+1)*SI/pw_std)
                st=0;
                tail=0;
            end
        end
    end
end
%统计串扰幅值分布
m=max(st_Vpp);
for i=1:m
    y(i)=size(find(st_Vpp==i),2);
end
y=y/sum(y);
x=1:m;
p_crosstalk=mean(y(2:end)/y(1:end-1));

loc=figure;
plot(x,y,'LineWidth',2,'Color',[0,0,0.8]);
hold on
plot(x,y(1)*p_crosstalk.^(x-1),'LineWidth',2,'Color',[0.8,0,0]);
xlabel('Vpp/Vpp_{std}');
ylabel('PMF');
legend('Exp data','Logarithmic Function Fit')
addpath("F:\USTC_reasearch_topic\UOWC\Long Distance");
plotformat(loc)
rmpath("F:\USTC_reasearch_topic\UOWC\Long Distance");