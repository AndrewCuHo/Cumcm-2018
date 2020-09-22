function [ Chrom ] = coding( CD)
%CODING 此处显示有关此函数的摘要
%   此处显示详细说明

FarmSize=(size(CD));
R=FarmSize(1);
C=FarmSize(2);

W_num = R * C;%工序总个数

Num = zeros(1, 8); %某周期花费时间

for i=1:8
    Num(i) = C; %8个一周其工序个数
end
%个体编码
%这里我们考虑隔一个为第一/第二工序；
%8）为周期时间迭代，（8为工序
%随机生成一个周期工序:
Chrom = zeros(R, C*2); 
for i=1: R
    t_num = Num;
    for j=1 : C
        %此处指第二段工序,          0,1交叉为工序
        val = unidrnd(C);%
        val_t = [zeros(1, C/2), ones(1, C/2)]';
        
        val = val(randperm(numel(val)));
        val_t = val_t(randperm(size(val_t,1)),:);
         %对前时间矩阵正常随机，  
        %Chrom(i, j) = val;
        %产生后半段工序
        Chrom(i, j + C) = val_t(j);
    end
 
end
   Chrom = [CD, Chrom];
        
        
        
        
        
    

end

