function P=Caltime(DC,PNumber)

% ����˵����          ���ݻ���S,������ȹ���P,ʱ�������func��˴������������
% ���������
%        S           Ϊ����  
%        PNumber     Ϊ���ӹ����� 
%        alpha       alphaΪÿ���һ����������ʱ��   
% �������: 
%        P           Ϊ����ĵ��ȹ��� 
%                    ��������801��ʾ����8�Ĺ���01
[R, C] = size(DC);

WNumber=length(DC);%�����ܸ���
WNumber=WNumber/2;

%ȡ�������ȡ�����һ��
DC = DC(1:WNumber, end);

%��ʼ��
temp=zeros(1,PNumber);
P=zeros(1,WNumber);

%�������ɵ��ȹ���
for i=1: WNumber 
    
   %�����+1
  temp(DC(i))=temp(DC(i))+1; 
  if temp(i)> 1
      
      temp(i) = 1;
  end
  if temp(i) <=0 
     
      temp(i)=0;
  end
  
  P(i)=DC(i)*100+temp(DC(i));
  
end


