clc;
clear all;
close all;
format short g ;

%for h=1:28     %��һ�����
% for h=1:26    %�ڶ������
glo_cnt = 0;
List_store = [1, 2, 3, 4, 5, 6, 7, 8];
P = zeros(1, 500); %ȷ����ʼ���й���Ϊ0

Final_CNC_1 = [];
Final_liao_1 = [1];%��һ������
extra_diejia = [0]; %ֻ��һ������������ֻ��������

Final_CNC_2 = [];
Final_liao_2 = [1];%��һ������

Flag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%�ȱ���



for h=1:25
N=8;%CNC����
glo_cnt = glo_cnt + 1;
%����洢������ʱ����Լ���һ��,������ʱ��һ��

    
   
M=500;%Ԥ����Ⱥ��ģ 
D=7;%��������
DC=D;%����������
pc=0.65;%�������
pm=0.24;%0.05;%�������
%wait_cnc = 25;%CNC�ȴ���ʱ��,��һ��
%wait_cnc = 30;%CNC�ȴ���ʱ��,�ڶ���
wait_cnc = 25;%CNC�ȴ���ʱ��,������

% alpha = 400;%��һ������ʱ��  1��
% beta = 378;%�ڶ�������ʱ��   1��

% alpha = 280;%��һ������ʱ��  2��
% beta = 500;%�ڶ�������ʱ��   2��

alpha = 455;%��һ������ʱ��  3��
beta = 182;%�ڶ�������ʱ��   3��

%CD ��CityDistance��д����ʾ���м����, ���۾���,��һ��


% CD = [0, 28, 48, 48, 61, 61, 74, 74;
%         31, 0, 51, 51, 64, 64, 77, 77;
%         48, 48, 0, 28, 48, 48, 61, 61;
%         48, 48, 28, 0, 48, 48, 61, 61;
%         51, 51, 31, 51, 0, 51, 64, 64;
%         64, 64, 51, 51, 31, 0, 51, 51;
%         74, 74, 61, 61, 48, 48, 0, 28;
%         77, 77, 64, 64, 51, 51, 31, 0
%         ]; % ��һ����Ϣ

% CD = [0, 30, 53, 53, 71, 71, 89, 89;
%         35, 0, 58, 58, 76, 76, 94, 94;
%         53, 53, 0, 30, 53, 53, 71, 71;
%         58, 58, 35, 0, 58, 58, 76, 76;
%         71, 71, 53, 53, 0, 30, 53, 53;
%         76, 76, 58, 58, 35, 0, 58, 58;
%         89, 89, 71, 71, 53, 30, 0, 30;
%         94, 94, 76, 76, 58, 58, 35, 0
%         ];%�ڶ�����Ϣ
%  
    CD = [0, 27, 45, 45, 59, 59, 73, 73;
        32, 0, 50, 50, 64, 64, 78, 78;
        45, 45, 0, 27, 45, 45, 59, 59;
        50, 50, 32, 0, 50, 50, 64, 64;
        59, 59, 45, 45, 0, 27, 45, 45;
        64, 64, 50, 50, 32, 0, 50, 50;
        73, 73, 59, 59, 45, 45, 0, 27;
        78, 78, 64, 64, 50, 50, 32, 0
        ];%��������Ϣ


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
min_store = [1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ʼ����ȺFarm������������ɷ���
Temp2=zeros(M,1);%1�����ڴ�Ŷ���ľ���
Temp1=zeros(M,N);%N�����ڴ��·��

for index=1:M %ѭ��M��
Temp1(index,:)=randperm(N);%����һ��1��N�������
end
% ��Ⱥ��ģMΪ����  ���еĸ���N+1Ϊ���� ,���һ�з���Ӧ��
Farm=[Temp1,Temp2];%����ĺϲ�
List = ones(size(Farm));%LOgistc����Ƿ���ϴ

% if h > 1
%     List_t = [];
%     for k=1:size(min_store, 1)
%     List_t = [List_t, List_store(min_store(k))];%�洢��һ���ڵı仯����������ͼ�
%     end
%     List = [List(:, 1:N), List(:, N+1) .* List_t];
% end


%���������,˫�������
CD = coding(CD);

while 1    %ѭ����1  ��ʼ��ʶ��
% %debug
% FarmSize=(size(Farm));
% R=FarmSize(1);
% fprintf('\n��%.f�ε�����ʼ����Ⱥ�ĸ���%.f\n',D-DC+1,R);
% %debug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������Ⱥ������ÿ���������Ӧ��
[Farm, juli, fenbu, P_last, FIR, SEC, x_min] = new_AdaptFun(Farm,CD, min_store, P, alpha, beta);%������Ӧ�Ⱥ���AdaptFun��������Ӧ��,����tig

P = P_last;

%�ֲ�Ϊ����������ܶ�

% if DC==0  %������������ﵽ��Ԥ��ֵ
%     %������ ������·��������·���Ĵ���
%     %beta_1 = 560 ;%����ƫ��������һ��
%     %beta_1 = 580 ;%����ƫ�������ڶ���
%     beta_1 = 545 ;%����ƫ������������
%     result=min(Farm(:,N+1));%�������еĵ�N+1����Ѱ����Сֵ
%     result_fil = result + beta_1;
%     fprintf('\n����·���Ĵ���:\n %.1f',result_fil);
%     
%     fprintf('\n\n\n����·��:\n');
%    [r,c]=find(Farm(:,N+1)==result);%Ѱ�Ҵ����Сֵ�Ķ�Ӧ�к�
%    for index=1:N  %��ӡ·��
%    fprintf('%.0f  ',Farm(r(1),index));
%    %r���ܳ��ֶ���������ж���⣩����r(1)
%           
%    end
%    fprintf('\n');
%    
%    %��һ������·���Ĵ���:4802.0
%    %����·��:8  6  5  3  4  2  1  7
%    %�ڶ��Σ�5005.0
%    %����·��:4  1  2  8  7  6  5  3
%    %�ڶ��Σ�4696.0
%    %����·��:1  4  3  8  7  6  5  2  
%     
%    
%  
%     break;%��ô����ѭ����1
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ѡ��

Farm=Selection(Farm,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��delta
t_2 = zeros(500, 8);
delta_t  = 1 ./ (1 + zscore(Farm(:, N+1)) .* zscore(Farm(:, N+1)));
delta = cat(2,t_2, delta_t);

min_t = find(Farm == min(Farm(:,N+1)));%�õ�

List = delta .* wait_cnc .* 210;    %����
Farm = List + Farm;

result_t1 = min(Farm(:,N+1))        %�����������ֵ

extra_diejia = [extra_diejia, result_t1 + extra_diejia(end)];

Final_liao_1 = [Final_liao_1, extra_diejia(end)];
Final_liao_2 = [Final_liao_2, extra_diejia(end)];


fprintf('\n\n\n����һCNC���: \n')

[r,c]=find(Farm(:,N+1)==result_t1);%Ѱ�Ҵ����Сֵ�Ķ�Ӧ�к�
   for index=1:N  %��ӡ·��
   fprintf('%.0f ',Farm(r(1),index));
   Final_CNC_1 = [Final_CNC_1, Farm(r(1),index)];
   end
   
fprintf('\n\n\n�����CNC���: \n')
   for index=1:N  %��ӡ·��
        try
            fprintf('%.0f ',Farm(r(2),index));
            Final_CNC_2 = [Final_CNC_2, Farm(r(2),index)];
        catch
            fprintf('%.0f ',Farm(r(1),index));
            Final_CNC_2 = [Final_CNC_2, Farm(r(1),index)];
        end
   end
   
   
if DC==0  %������������ﵽ��Ԥ��ֵ
    %������ ������·��������·���Ĵ���
    %beta_1 = 560 ;%����ƫ��������һ��
    %beta_1 = 580 ;%����ƫ�������ڶ���
    %beta_1 = 545 ;%����ƫ������������
    result=min(Farm(:,N+1));%�������еĵ�N+1����Ѱ����Сֵ
    result_fil = result;% + beta_1;
    fprintf('\n����·���Ĵ���:\n %.1f',result_fil);
    
    extra_diejia = [extra_diejia, result_fil + extra_diejia(end)];
     
    Final_liao_1 = [Final_liao_1, extra_diejia(end)];
    Final_liao_2 = [Final_liao_2, extra_diejia(end)];
    
    fprintf('\n\n\n����·��:\n');
   [r,c]=find(Farm(:,N+1)==result);%Ѱ�Ҵ����Сֵ�Ķ�Ӧ�к�
   for index=1:N  %��ӡ·��
   fprintf('%.0f  ',Farm(r(1),index));
   
   Final_CNC_1 = [Final_CNC_1, Farm(r(1),index)];
   
   Final_CNC_2 = [Final_CNC_2, Farm(r(1),index)];
   
   %r���ܳ��ֶ���������ж���⣩����r(1)
          
   end
   fprintf('\n');
   
   %��һ������·���Ĵ���:4802.0
   %����·��:8  6  5  3  4  2  1  7
   %�ڶ��Σ�5005.0
   %����·��:4  1  2  8  7  6  5  3
   %�ڶ��Σ�4696.0
   %����·��:1  4  3  8  7  6  5  2  
    
   
 
    break;%��ô����ѭ����1
end


% %debug
% FarmSize=(size(Farm));
% R=FarmSize(1);
% fprintf('\nѡ�������������Ⱥ�ĸ���%.f\n',R);
% 
% 
%     result=min(Farm(:,N+1));%�������еĵ�N+1����Ѱ����Сֵ
%     fprintf('\n����·���Ĵ���:\n %.1f',result);
% %debug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����
Farm=Crossover(Farm,pc);
% %debug
% FarmSize=(size(Farm));
% R=FarmSize(1);
% fprintf('\n���������������Ⱥ�ĸ���%.f\n',R);
% %debug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����
Farm=Mutation(Farm,pm);
% %debug
% FarmSize=(size(Farm));
% R=FarmSize(1);
% fprintf('\n��%.f�ε�����������Ⱥ�ĸ���%.f\n',D-DC+1,R);
% %debug

DC=DC-1; %����������
%������ת����������Ⱥ����һ��

for s=1:7
    
 if (Farm(x_min(s),9) > alpha)&&(P(x_min(s)) == 0)
            P(x_min(s)) = 1;
            Farm(x_min(s),9)= Farm(x_min(s),9) + alpha;
            FIR = [FIR, Farm(x_min(s),9)./11];%%%%%%%%c������
            if(size(FIR, 2) > 8)
                FIR = FIR(1,1:8)
            end
 end
    
     if (Farm(x_min(s),9) > beta)&&(P(x_min(s)) == 1)
            P(x_min(s)) = 0;
            Farm(x_min(s),9)= Farm(x_min(s),9) + beta;
            SEC = [SEC, Farm(x_min(s),9)./11];%%%%%%%%%������
            if(size(SEC, 2) > 8)
                SEC = SEC(1,1:8)
            end
         
     end
end

     

fprintf('\n����·���Ĵ���: %.1f',juli);            %˫ʱ����

fprintf('\n����һ����·���Ĵ���: %.1f',FIR);

for i = 1:size(FIR, 2)
    Final_liao_1 = [Final_liao_1, FIR(i) +  extra_diejia(end)];
end

fprintf('\n���������·���Ĵ���: %.1f',SEC);

for i = 1:size(SEC, 2)
    Final_liao_2 = [Final_liao_2, SEC(i) +  extra_diejia(end)];
end

min_store = min_t(1,1);%�洢��Сֵ
List_store = min_t;

if (Final_liao_1(end)>28800)||(Final_liao_1(end)>28800)
    Flag = 1;
    break;

end    %ѭ����1  ������ʶ��

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%msgbox('���н���','YYM');
if rem(glo_cnt, 2) == 0
   
    
    List_store = [];
   %P =  Caltime(DC,250);            %�������Ϊ250������Ϊһ���ϵ㣬ȷ�����й����Ѿ�+1
                                        %   ��������801��ʾ����8�Ĺ���01
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
