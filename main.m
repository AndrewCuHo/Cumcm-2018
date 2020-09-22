clc;
clear all;
close all;
format short g ;

%for h=1:28     %第一组迭代
% for h=1:26    %第二组迭代
glo_cnt = 0;
List_store = [1, 2, 3, 4, 5, 6, 7, 8];
P = zeros(1, 500); %确保初始所有工序为0

Final_CNC_1 = [];
Final_liao_1 = [1];%第一个不算
extra_diejia = [0]; %只有一个！！！！！只算最终优

Final_CNC_2 = [];
Final_liao_2 = [1];%第一个不算

Flag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%先编码



for h=1:25
N=8;%CNC个数
glo_cnt = glo_cnt + 1;
%方面存储，两次时间后大约清除一次,和两次时间一致

    
   
M=500;%预设种群规模 
D=7;%迭代次数
DC=D;%迭代计数器
pc=0.65;%交叉概率
pm=0.24;%0.05;%变异概率
%wait_cnc = 25;%CNC等待的时间,第一组
%wait_cnc = 30;%CNC等待的时间,第二组
wait_cnc = 25;%CNC等待的时间,第三组

% alpha = 400;%第一道工序时间  1组
% beta = 378;%第二道工序时间   1组

% alpha = 280;%第一道工序时间  2组
% beta = 500;%第二道工序时间   2组

alpha = 455;%第一道工序时间  3组
beta = 182;%第二道工序时间   3组

%CD 是CityDistance缩写，表示城市间距离, 代价矩阵,第一组


% CD = [0, 28, 48, 48, 61, 61, 74, 74;
%         31, 0, 51, 51, 64, 64, 77, 77;
%         48, 48, 0, 28, 48, 48, 61, 61;
%         48, 48, 28, 0, 48, 48, 61, 61;
%         51, 51, 31, 51, 0, 51, 64, 64;
%         64, 64, 51, 51, 31, 0, 51, 51;
%         74, 74, 61, 61, 48, 48, 0, 28;
%         77, 77, 64, 64, 51, 51, 31, 0
%         ]; % 第一组信息

% CD = [0, 30, 53, 53, 71, 71, 89, 89;
%         35, 0, 58, 58, 76, 76, 94, 94;
%         53, 53, 0, 30, 53, 53, 71, 71;
%         58, 58, 35, 0, 58, 58, 76, 76;
%         71, 71, 53, 53, 0, 30, 53, 53;
%         76, 76, 58, 58, 35, 0, 58, 58;
%         89, 89, 71, 71, 53, 30, 0, 30;
%         94, 94, 76, 76, 58, 58, 35, 0
%         ];%第二组信息
%  
    CD = [0, 27, 45, 45, 59, 59, 73, 73;
        32, 0, 50, 50, 64, 64, 78, 78;
        45, 45, 0, 27, 45, 45, 59, 59;
        50, 50, 32, 0, 50, 50, 64, 64;
        59, 59, 45, 45, 0, 27, 45, 45;
        64, 64, 50, 50, 32, 0, 50, 50;
        73, 73, 59, 59, 45, 45, 0, 27;
        78, 78, 64, 64, 50, 50, 32, 0
        ];%第三组信息


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
min_store = [1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%初始化种群Farm，采用随机生成方法
Temp2=zeros(M,1);%1列用于存放定义的距离
Temp1=zeros(M,N);%N列用于存放路径

for index=1:M %循环M次
Temp1(index,:)=randperm(N);%生成一个1到N随机序列
end
% 种群规模M为行数  城市的个数N+1为列数 ,最后一列放适应度
Farm=[Temp1,Temp2];%矩阵的合并
List = ones(size(Farm));%LOgistc解决是否清洗

% if h > 1
%     List_t = [];
%     for k=1:size(min_store, 1)
%     List_t = [List_t, List_store(min_store(k))];%存储上一周期的变化，及能量最低级
%     end
%     List = [List(:, 1:N), List(:, N+1) .* List_t];
% end


%先随机编码,双随机编码
CD = coding(CD);

while 1    %循环体1  开始标识符
% %debug
% FarmSize=(size(Farm));
% R=FarmSize(1);
% fprintf('\n第%.f次迭代开始，种群的个数%.f\n',D-DC+1,R);
% %debug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%评价种群，计算每个个体的适应度
[Farm, juli, fenbu, P_last, FIR, SEC, x_min] = new_AdaptFun(Farm,CD, min_store, P, alpha, beta);%调用适应度函数AdaptFun，生成适应度,加了tig

P = P_last;

%分布为旺假设概率密度

% if DC==0  %如果迭代次数达到了预设值
%     %输出结果 ：最优路径、最优路径的代价
%     %beta_1 = 560 ;%常数偏移量，第一次
%     %beta_1 = 580 ;%常数偏移量，第二次
%     beta_1 = 545 ;%常数偏移量，第三次
%     result=min(Farm(:,N+1));%在所有行的第N+1列中寻找最小值
%     result_fil = result + beta_1;
%     fprintf('\n最优路径的代价:\n %.1f',result_fil);
%     
%     fprintf('\n\n\n最优路径:\n');
%    [r,c]=find(Farm(:,N+1)==result);%寻找存放最小值的对应行号
%    for index=1:N  %打印路径
%    fprintf('%.0f  ',Farm(r(1),index));
%    %r可能出现多行情况（有多个解），故r(1)
%           
%    end
%    fprintf('\n');
%    
%    %第一次最优路径的代价:4802.0
%    %最优路径:8  6  5  3  4  2  1  7
%    %第二次：5005.0
%    %最优路径:4  1  2  8  7  6  5  3
%    %第二次：4696.0
%    %最优路径:1  4  3  8  7  6  5  2  
%     
%    
%  
%     break;%那么跳出循环体1
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%选择

Farm=Selection(Farm,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%设delta
t_2 = zeros(500, 8);
delta_t  = 1 ./ (1 + zscore(Farm(:, N+1)) .* zscore(Farm(:, N+1)));
delta = cat(2,t_2, delta_t);

min_t = find(Farm == min(Farm(:,N+1)));%用到

List = delta .* wait_cnc .* 210;    %超参
Farm = List + Farm;

result_t1 = min(Farm(:,N+1))        %输出迭代最优值

extra_diejia = [extra_diejia, result_t1 + extra_diejia(end)];

Final_liao_1 = [Final_liao_1, extra_diejia(end)];
Final_liao_2 = [Final_liao_2, extra_diejia(end)];


fprintf('\n\n\n工序一CNC编号: \n')

[r,c]=find(Farm(:,N+1)==result_t1);%寻找存放最小值的对应行号
   for index=1:N  %打印路径
   fprintf('%.0f ',Farm(r(1),index));
   Final_CNC_1 = [Final_CNC_1, Farm(r(1),index)];
   end
   
fprintf('\n\n\n工序二CNC编号: \n')
   for index=1:N  %打印路径
        try
            fprintf('%.0f ',Farm(r(2),index));
            Final_CNC_2 = [Final_CNC_2, Farm(r(2),index)];
        catch
            fprintf('%.0f ',Farm(r(1),index));
            Final_CNC_2 = [Final_CNC_2, Farm(r(1),index)];
        end
   end
   
   
if DC==0  %如果迭代次数达到了预设值
    %输出结果 ：最优路径、最优路径的代价
    %beta_1 = 560 ;%常数偏移量，第一次
    %beta_1 = 580 ;%常数偏移量，第二次
    %beta_1 = 545 ;%常数偏移量，第三次
    result=min(Farm(:,N+1));%在所有行的第N+1列中寻找最小值
    result_fil = result;% + beta_1;
    fprintf('\n最优路径的代价:\n %.1f',result_fil);
    
    extra_diejia = [extra_diejia, result_fil + extra_diejia(end)];
     
    Final_liao_1 = [Final_liao_1, extra_diejia(end)];
    Final_liao_2 = [Final_liao_2, extra_diejia(end)];
    
    fprintf('\n\n\n最优路径:\n');
   [r,c]=find(Farm(:,N+1)==result);%寻找存放最小值的对应行号
   for index=1:N  %打印路径
   fprintf('%.0f  ',Farm(r(1),index));
   
   Final_CNC_1 = [Final_CNC_1, Farm(r(1),index)];
   
   Final_CNC_2 = [Final_CNC_2, Farm(r(1),index)];
   
   %r可能出现多行情况（有多个解），故r(1)
          
   end
   fprintf('\n');
   
   %第一次最优路径的代价:4802.0
   %最优路径:8  6  5  3  4  2  1  7
   %第二次：5005.0
   %最优路径:4  1  2  8  7  6  5  3
   %第二次：4696.0
   %最优路径:1  4  3  8  7  6  5  2  
    
   
 
    break;%那么跳出循环体1
end


% %debug
% FarmSize=(size(Farm));
% R=FarmSize(1);
% fprintf('\n选择操作结束，种群的个数%.f\n',R);
% 
% 
%     result=min(Farm(:,N+1));%在所有行的第N+1列中寻找最小值
%     fprintf('\n最优路径的代价:\n %.1f',result);
% %debug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%交叉
Farm=Crossover(Farm,pc);
% %debug
% FarmSize=(size(Farm));
% R=FarmSize(1);
% fprintf('\n交叉操作结束，种群的个数%.f\n',R);
% %debug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%变异
Farm=Mutation(Farm,pm);
% %debug
% FarmSize=(size(Farm));
% R=FarmSize(1);
% fprintf('\n第%.f次迭代结束，种群的个数%.f\n',D-DC+1,R);
% %debug

DC=DC-1; %迭代计数器
%重新跳转到“评价种群”这一步

for s=1:7
    
 if (Farm(x_min(s),9) > alpha)&&(P(x_min(s)) == 0)
            P(x_min(s)) = 1;
            Farm(x_min(s),9)= Farm(x_min(s),9) + alpha;
            FIR = [FIR, Farm(x_min(s),9)./11];%%%%%%%%c超参数
            if(size(FIR, 2) > 8)
                FIR = FIR(1,1:8)
            end
 end
    
     if (Farm(x_min(s),9) > beta)&&(P(x_min(s)) == 1)
            P(x_min(s)) = 0;
            Farm(x_min(s),9)= Farm(x_min(s),9) + beta;
            SEC = [SEC, Farm(x_min(s),9)./11];%%%%%%%%%超参数
            if(size(SEC, 2) > 8)
                SEC = SEC(1,1:8)
            end
         
     end
end

     

fprintf('\n最优路径的代价: %.1f',juli);            %双时不管

fprintf('\n工序一最优路径的代价: %.1f',FIR);

for i = 1:size(FIR, 2)
    Final_liao_1 = [Final_liao_1, FIR(i) +  extra_diejia(end)];
end

fprintf('\n工序二最优路径的代价: %.1f',SEC);

for i = 1:size(SEC, 2)
    Final_liao_2 = [Final_liao_2, SEC(i) +  extra_diejia(end)];
end

min_store = min_t(1,1);%存储最小值
List_store = min_t;

if (Final_liao_1(end)>28800)||(Final_liao_1(end)>28800)
    Flag = 1;
    break;

end    %循环体1  结束标识符

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%msgbox('运行结束','YYM');
if rem(glo_cnt, 2) == 0
   
    
    List_store = [];
   %P =  Caltime(DC,250);            %最大数量为250，两次为一个断点，确保所有工序已经+1
                                        %   比如数字801表示工件8的工序01
end


end

  if Flag == 1
    break;
end
    
end

xlswrite('C:\Users\ZCH\Desktop\2018-B-Chinese\D_3_CNC1_xiaohao.xlsx',Final_liao_1');
xlswrite('C:\Users\ZCH\Desktop\2018-B-Chinese\D_3_CNC1_track.xlsx',Final_CNC_1');
xlswrite('C:\Users\ZCH\Desktop\2018-B-Chinese\D_3_CNC2_xiaohao.xlsx',Final_liao_2');
xlswrite('C:\Users\ZCH\Desktop\2018-B-Chinese\D_3_CNC2_track.xlsx',Final_CNC_2');
