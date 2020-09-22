function [Farm, juli, fenbu, P_last, FIR, SEC, x_min, juli_1, x_min_1] = new_AdaptFun(Farm,CD, min_store, P, alpha1, beita)  %P为计算的调度
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
%依照排队论给的适应度函数
%计算Farm现在的行数R ，列数C
FarmSize=(size(Farm));
R=FarmSize(1);
C=FarmSize(2);
juli = [];
juli_1 = [];
FIR=[];
SEC = [];
x_min = [];
x_min_1 = [];
[m, n] = find(Farm, Farm(min_store));
m = m(1, 1);
n = n(1, end);
%计算适应度时间，规划放在最后一列
%第一组数据

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%第一组参数
% miu_1 = 378;
% 
% t5_all = 1/61 + 1/64 + 1/48 + 1/51;
% t6_all = 1/61 + 1/64 + 1/48 + 1/51;
% t7_all = 1/74 + 1/77 + 1/61 + 1/64;
% t8_all = 1/74 + 1/77 + 1/61 + 1/64;

% fenzi5 =1 / ((400) * 4 + miu_1 +61 + 64 +48 +51);
% fenzi6 = 1 / (400 * 4 + miu_1 + 61 + 64 + 48 + 51);
% fenzi7 = 1 / (400 * 4 + miu_1 + 74 + 77 +61 +64);
% fenzi8 = fenzi7;

%T = ((61 + 64 + 48 + 51) * 2 + (74 + 77 + 61 + 64) * 2) / 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% miu_2 = 500;
% t5_all = 1/71 + 1/76 + 1/53 + 1/58;
% t6_all = 1/71 + 1/76 + 1/53 + 1/58;
% t7_all = 1/89 + 1/94 + 1/71 + 1/76;
% t8_all = 1/89 + 1/94 + 1/71 + 1/76;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

miu_3 = 182;

t5_all = 1/59 + 1/64 + 1/45 + 1/50;
t6_all = 1/59 + 1/64 + 1/45 + 1/50;
t7_all = 1/73 + 1/78 + 1/59 + 1/64;
t8_all = 1/73 + 1/78 + 1/59 + 1/64;


fenmu = t5_all + t6_all + t7_all + t8_all;

fenzi5 =1 / ((455) * 4 + miu_3 +59 + 64 +45 +50);
fenzi6 = fenzi5;
fenzi7 = 1 / (455 * 4 + miu_3 + 73 + 78 +59 +64);
fenzi8 = fenzi7;

w5 = fenzi5 / fenmu;
w6 = fenzi6 / fenmu;
w7 = fenzi7 / fenmu;
w8 = fenzi8 / fenmu;

T = ((59 + 64 + 45 + 50) * 2 + (73 + 78 + 59 + 64) * 2) / 4;

lamda5 = 1 / (T * w5); 
lamda6 = 1 / (T * w6);
lamda7 = 1 / (T * w7);
lamda8 = 1 / (T * w8);

%扰动系数
alpha = 0.00009;
%第一组
%t_1 = alpha*abs(((1/lamda5 - 1/miu_1)^2 + (1/lamda6 - 1/miu_1)^2 + (1/lamda7 - 1/miu_1)^2 + (1/lamda8 - 1/miu_1)^2)- T * (lamda5 + lamda6 + lamda7 + lamda8));
%第二组
%t_1 = alpha*abs(((1/lamda5 - 1/miu_2)^2 + (1/lamda6 - 1/miu_2)^2 + (1/lamda7 - 1/miu_2)^2 + (1/lamda8 - 1/miu_2)^2)- T * (lamda5 + lamda6 + lamda7 + lamda8));

t_1 = alpha*abs(((1/lamda5 - 1/miu_3)^2 + (1/lamda6 - 1/miu_3)^2 + (1/lamda7 - 1/miu_3)^2 + (1/lamda8 - 1/miu_3)^2)- T * (lamda5 + lamda6 + lamda7 + lamda8));

fenbu = abs(((1/lamda5 - 1/miu_3)^2 + (1/lamda6 - 1/miu_3)^2 + (1/lamda7 - 1/miu_3)^2 + (1/lamda8 - 1/miu_3)^2));

%CD(:, 1:C) = CD(:, 1:C) .* fenbu;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%解码，计算工序已完成
%计算时间和工序大头喽

P_last = P;

P = P';


for i=1:R
    distance=0;
    %第一次
    %EW =  560;
    for j=1:C-1-1 
        distance=(distance+ CD(Farm(i,j),Farm(i,j+1)));
        if (m==i)
            juli = [juli, distance];
            x_min = [x_min, i];    
        end
        
        if (n==i)
            juli_1 = [juli_1, distance];
            x_min_1 = [x_min_1, i];    
        end
        
            %juli = [juli, ]
    end
    distance=(distance+CD(Farm(i,C-1),Farm(i,1)));
    Farm(i,C)= distance ;  %+ EW;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %调整超参数
%     if (distance > alpha1)&&(P(i) == 0)
%             P(i) = 1;
%             Farm(i,C)= Farm(i,C) + alpha1;
%             FIR = [FIR, distance];
%     end
%     
%      if (distance > beita)&&(P(i) == 1)
%             P(i) = 0;
%             Farm(i,C)= Farm(i,C) + beita;
%             SEC = [SEC, distance];
%          
%      end
    
end
   


end

