function [ Chrom ] = coding( CD)
%CODING �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

FarmSize=(size(CD));
R=FarmSize(1);
C=FarmSize(2);

W_num = R * C;%�����ܸ���

Num = zeros(1, 8); %ĳ���ڻ���ʱ��

for i=1:8
    Num(i) = C; %8��һ���乤�����
end
%�������
%�������ǿ��Ǹ�һ��Ϊ��һ/�ڶ�����
%8��Ϊ����ʱ���������8Ϊ����
%�������һ�����ڹ���:
Chrom = zeros(R, C*2); 
for i=1: R
    t_num = Num;
    for j=1 : C
        %�˴�ָ�ڶ��ι���,          0,1����Ϊ����
        val = unidrnd(C);%
        val_t = [zeros(1, C/2), ones(1, C/2)]';
        
        val = val(randperm(numel(val)));
        val_t = val_t(randperm(size(val_t,1)),:);
         %��ǰʱ��������������  
        %Chrom(i, j) = val;
        %�������ι���
        Chrom(i, j + C) = val_t(j);
    end
 
end
   Chrom = [CD, Chrom];
        
        
        
        
        
    

end

