function P=Caltime(DC,PNumber)

% 功能说明：          根据基因S,计算调度工序P,时间计算在func里，此处解码调整工序
% 输入参数：
%        S           为基因  
%        PNumber     为最大加工个数 
%        alpha       alpha为每组第一道工序所耗时间   
% 输出参数: 
%        P           为输出的调度工序 
%                    比如数字801表示工件8的工序01
[R, C] = size(DC);

WNumber=length(DC);%工序总个数
WNumber=WNumber/2;

%取工序基因，取基因的一半
DC = DC(1:WNumber, end);

%初始化
temp=zeros(1,PNumber);
P=zeros(1,WNumber);

%解码生成调度工序
for i=1: WNumber 
    
   %工序加+1
  temp(DC(i))=temp(DC(i))+1; 
  if temp(i)> 1
      
      temp(i) = 1;
  end
  if temp(i) <=0 
     
      temp(i)=0;
  end
  
  P(i)=DC(i)*100+temp(DC(i));
  
end


