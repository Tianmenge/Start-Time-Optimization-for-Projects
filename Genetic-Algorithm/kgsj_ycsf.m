% �Ŵ��㷨
% ��Լ���Ż�ת��Ϊ��Լ���Ż��������
% ����ͬС���Ŀ�ʼʱ��S��Ϊһ��Ⱦɫ�壬���г�ʼ�����룬����֮��õ���Ӧ��xֵ��yֵ��Ȼ���ж��Ƿ�����Լ������
tic;    % ��¼����ʱ��
t1 = clock;

T = 156;   % ����С����������ʱ��
I = 36;    % С����Ŀ
J = 195;   % ���������ܣ�
code = 8;  % ��8λ����������s���б���
ss = 120;  % sȡֵ��Χ������
pop_size = 36;              % ��Ⱥ��С
chromo_size = I*code;       % Ⱦɫ���С
generation_size = 10;     % ��������
cross_rate = 0.6;          % �������
mutate_rate = 0.01;        % �������
elitism = 2;
fitness_avg = zeros(generation_size,1);  % ÿһ������СĿ�꺯��ֵ
best_fitness = Inf;    % ����Ŀ�꺯��ֵ
best_generation = 0;   % ����Ŀ�꺯��ֵ���ֵĴ���
best_individual = zeros(1,chromo_size);    % ���Ÿ���
best_s = zeros(1,I);  % ���Ÿ����Ӧ�Ŀ���ʱ��

% ��ȡ����
data = csvread('data.csv',1,1);
M = data(:,1);
A = data(:,2);
demand = csvread('demand.csv',1);
P = demand(:,2);
B = demand(:,3);
G = demand(:,4);

gene = zeros(pop_size,chromo_size);
s = zeros(pop_size,I);
pop_X = cell(1,pop_size);
pop_Y = cell(1,pop_size);

% ���ȶԳ�ʼ���ĵ�1��������Ӧ�ȡ�����ȣ�Ȼ����е���
% ���ν���ѡ�񡢽��桢����Ȳ�������������Ӧ�ȡ������

%----------------------------------��ʼ��-----------------------------------
% ��ʼ����Ⱥ�����ͣ��Բ�ͬС���Ŀ���ʱ�����8λ�����������б���

for k=1:pop_size
    for i=1:chromo_size
        gene(k,i) = round(rand());
    end
end
   
%-----------------------------------����------------------------------------
% �Գ�ʼ���Ļ����ͽ��н��룬�õ���Ӧ�Ĳ�ͬС���Ŀ���ʱ���Լ���Ӧ��x��y

for k=1:pop_size
    for i=1:I
        b = gene(k,(code*i-code+1):(code*i));
        a = bin2dec(num2str(b));
        s(k,i) = 1 + round(a*(ss-1)/(2^code-1));
    end
end

% ����s�õ���Ӧ��Y��X��Y��һ��156*36��0-1����y(ji)��ʾ�ڵ�j��С��i�Ƿ��ڽ���ͬʱ����ͬһ���ڽ�С����������15��
% X��һ��195*36��0-1����x(ji)��ʾ�ڵ�j��С��i�Ƿ��Ѿ����,�ڵ�T��֮��XȡֵΪ1
for k=1:pop_size
    temp_X = zeros(J,I);
    temp_Y = zeros(J,I);
    ini_pos = s(k,:);  % ini_pos��ȡֵ������[0,112]����֤������С����������ʱ��Ϊ��156��
    for i=1:I
        duration = ini_pos(i)+18+M(i);
        temp_X(duration:J,i) = 1;
        temp_Y(ini_pos(i):duration,i) = 1;
    end
    pop_X{k} = temp_X;
    pop_Y{k} = temp_Y;
end

% % {
%--------------------------------������Ӧ��----------------------------------
%������Ⱥ������Ӧ�ȣ��Բ�ͬ���Ż�Ŀ�꣬�˴���Ҫ��д
%pop_size: ��Ⱥ��С
%chromo_size: Ⱦɫ�峤��

fitness_value=zeros(1,pop_size);
for k=1:pop_size
    temp_X = pop_X{k};
    temp_Y = 15-sum(pop_Y{k},2);
    
    % ����С�����ͼ�����Ӧ��ʵ�ʽ������
    P_real = zeros(size(P));
    B_real = zeros(size(B));
    G_real = zeros(size(G));
    for j=1:J
        for i=1:4
            P_real(j) = P_real(j) + temp_X(j,i)*A(i);
        end
        for i=5:12
            B_real(j) = B_real(j) + temp_X(j,i)*A(i);
        end
        for i=13:36
            G_real(j) = G_real(j) + temp_X(j,i)*A(i);
        end
        fitness_value(k) = fitness_value(k) + 10^8*((min(0,B_real(j)-B(j)))^2 + (min(0,P_real(j)-P(j)))^2 + (min(0,G_real(j)-G(j)))^2 + (min(0,temp_Y(j)))^2);
    end
end

%----------------------�Ը��尴����Ӧ�ȴ�С��������---------------------------
%�Ը��尴��Ӧ�ȴ�С�������򣬲��ұ�����Ѹ���
%pop_size: ��Ⱥ��С
%chromo_size: Ⱦɫ�峤��

[fitness_value,index] = sort(fitness_value);
temp_X = cell(size(pop_X));
temp_Y = cell(size(pop_Y));
temp_s = zeros(size(s));
temp_gene = zeros(size(gene));
for k=1:pop_size
    temp_X{k} = pop_X{index(k)};
    temp_Y{k} = pop_Y{index(k)};
    temp_s(k,:) = s(index(k),:);
    temp_gene(k,:) = gene(k,:);
end
pop_X = temp_X;
pop_Y = temp_Y;
s = temp_s;
gene = temp_gene;

% ����Ҫ��fitness_value��С��pop�����Զ�fitness_value��С����������򣬵õ�rank����rank����scale�����õ�scale_value
fitness_table = zeros(1,pop_size);  % �ۻ�
rank = 1:1:pop_size;
temp = 1./(rank.^0.5);
scale_value = temp.*sum(temp)/pop_size;

for i=1:pop_size
    if i==1
        fitness_table(i) = fitness_table(i) + scale_value(i);    
    else
        fitness_table(i) = fitness_table(i-1) + scale_value(i);
    end
end
fitness_avg(1) = sum(fitness_value)/pop_size;

if fitness_value(1) < best_fitness
    fitness_value(1)
    best_fitness = fitness_value(1);
    best_generation = 1;
    best_s = s(1,:);
end
%}

%========================�������� generation_size===========================
for generation=2:generation_size
    generation

%---------------------------------ѡ�����----------------------------------
% �����޳�������Լ����������Ⱥ����ʣ�µ���Ⱥ������Ӧ�Ƚ�����������µ�36����Ⱥ
temp_X = {};
temp_Y = {};
temp_s = [];
temp_gene = [];
count = 0;
temp_value = [];
% �޳�������Լ����������Ⱥ
for k=1:pop_size
    temp1 = pop_Y{k};
    temp2 = pop_X{k};
    if sum( sum(temp1,2) <= 15 ) == J
        % ����С�����ͼ�����Ӧ��ʵ�ʽ������
%         P_real = zeros(size(P));
%         B_real = zeros(size(B));
%         G_real = zeros(size(G));
%         for j=1:J
%             for i=1:4
%                 P_real(j) = P_real(j) + temp2(j,i)*A(i);
%             end
%             for i=5:12
%                 B_real(j) = B_real(j) + temp2(j,i)*A(i);
%             end
%             for i=13:36
%                 G_real(j) = G_real(j) + temp2(j,i)*A(i);
%             end   
%         end
%         if (sum(P_real>=P) == J) || (sum(B_real>=B) == J) || (sum(G_real>=G) == J)
            count = count + 1;
%             temp_value(count) = 10^8*(sum((min(0,P_real-P)).^2) + sum((min(0,B_real-B)).^2) + sum((min(0,G_real-G)).^2));
            temp_X{count} = pop_X{k};
            temp_Y{count} = pop_Y{k};
            temp_s(count,:) = s(k,:);
            temp_gene(count,:) = gene(k,:);
%         end
    end
end
count
%{
% ������Ӧ�ȴ�С��������
[temp_value,index] = sort(temp_value);
temp1_X = cell(size(temp_X));
temp1_Y = cell(size(temp_Y));
temp1_s = zeros(size(temp_s));
temp1_gene = zeros(size(temp_gene));
for k=1:count
    temp1_X{k} = temp_X{index(k)};
    temp1_Y{k} = pop_Y{index(k)};
    temp1_s(k,:) = temp_s(index(k),:);
    temp1_gene(k,:) = temp_gene(k,:);
end
temp_X = temp1_X;
temp_Y = temp1_Y;
temp_s = temp1_s;
temp_gene = temp1_gene;

fitness_avg(generation) = sum(fitness_value)/count;

if fitness_value(1) < best_fitness
    best_fitness = fitness_value(1);
    best_generation = 1;
    best_s = s(1,:);
end
%}
% ���countС��pop_size����temp_X��ѭ��д��pop_X
if count < pop_size
    for i=1:count:pop_size        
        pop_X(i:count+i-1) = temp_X;
        pop_Y(i:count+i-1) = temp_Y;
        s(i:count+i-1,:) = temp_s;
        gene(i:count+i-1,:) = temp_gene;
    end
    for i=1:mod(pop_size,count)
        pop_X{pop_size-mod(pop_size,count)+i} = temp_X{i};
        pop_Y{pop_size-mod(pop_size,count)+i} = temp_Y{i};
        s(pop_size-mod(pop_size,count)+i,:) = temp_s(i,:);
        gene(pop_size-mod(pop_size,count)+i,:) = temp_gene(i,:);
    end
elseif count == 0
    continue;
else
    pop_X = temp_X;
    pop_Y = temp_Y;
    s = temp_s;
    gene = temp_gene;
end

% ���̶�ѡ������������ŵ�ǰ��������ֱ���Ŵ�����һ��
new_X = cell(size(pop_X));
new_Y = cell(size(pop_Y));
new_s = zeros(size(s));
new_gene = zeros(size(gene));

new_X(1:elitism) = pop_X(1:elitism);
new_Y(1:elitism) = pop_Y(1:elitism);
new_s(1:elitism,:) = s(1:elitism,:);
new_gene(1:elitism,:) = gene(1:elitism,:);

for i=elitism+1:pop_size
    r = rand() * fitness_table(pop_size);
    first = 1;
    last = pop_size;  
    idx = -1;
    % ʹ�ö��ַ�ѡ�����
    while (first <= last) && (idx == -1) 
        mid = round((last+first)/2);
        if r > fitness_table(mid)
            first = mid;
        elseif r < fitness_table(mid)
            last = mid;     
        else
            idx = mid;
            break;
        end
        if mid == first
            idx = last;
            break;
        else
            idx = first;
            break;
        end
   end
    new_X{i} = pop_X{idx};
    new_Y{i} = pop_Y{idx};
    new_s(i,:) = s(idx,:);
    new_gene(i,:) = gene(idx,:);
end

pop_X = new_X;
pop_Y = new_Y;
s = new_s;
gene = new_gene;

%----------------------------------�������----------------------------------
% ������������Ⱥ���е��㽻��������������Ϊcross_rate
% ��ѡ��Ҫ�����������Ⱥ��ѡ�����⽻��㣬�����֮��Ļ���ȫ������

for i=1+elitism:2:pop_size
    if rand() < cross_rate
        cross_pos = 1 + round(rand() * (chromo_size-1));
        if cross_pos == 1
            continue;
        end
        temp = gene(i,cross_pos:chromo_size);
        gene(i,cross_pos:chromo_size) = gene(i+1,cross_pos:chromo_size);
        gene(i+1,cross_pos:chromo_size) = temp;
    end
end

%------------------------------����֮�����----------------------------------
% �Գ�ʼ���Ļ����ͽ��н��룬�õ���Ӧ�Ĳ�ͬС���Ŀ���ʱ���Լ���Ӧ��x��y

for k=1+elitism:pop_size
    for i=1:I
        b = gene(k,(code*i-code+1):(code*i));
        a = bin2dec(num2str(b));
        s(k,i) = 1 + round(a*(ss-1)/(2^code-1));
    end
end

% ����s�õ���Ӧ��Y��X��Y��һ��156*36��0-1����y(ji)��ʾ�ڵ�j��С��i�Ƿ��ڽ���ͬʱ����ͬһ���ڽ�С����������15��
% X��һ��195*36��0-1����x(ji)��ʾ�ڵ�j��С��i�Ƿ��Ѿ����,�ڵ�T��֮��XȡֵΪ1
for k=1+elitism:pop_size
    temp_X = zeros(J,I);
    temp_Y = zeros(J,I);
    ini_pos = s(k,:);  % ini_pos��ȡֵ������[0,112]����֤������С����������ʱ��Ϊ��156��
    for i=1:I
        duration = ini_pos(i)+18+M(i);
        temp_X(duration:J,i) = 1;
        temp_Y(ini_pos(i):duration,i) = 1;
    end
    pop_X{k} = temp_X;
    pop_Y{k} = temp_Y;
end

%---------------------------------�������-----------------------------------
% �������������������mutate_rate����ѡ�����е�һ������б��죬0-1����

for i=1+elitism:pop_size
    if rand < mutate_rate
        mutate_pos = 1 + round(rand() * (chromo_size-1));
        gene(i,mutate_pos) = 1 - gene(i, mutate_pos);
    end
end

%------------------------------����֮�����----------------------------------
% �Գ�ʼ���Ļ����ͽ��н��룬�õ���Ӧ�Ĳ�ͬС���Ŀ���ʱ���Լ���Ӧ��x��y

for k=1+elitism:pop_size
    for i=1:I
        b = gene(k,(code*i-code+1):(code*i));
        a = bin2dec(num2str(b));
        s(k,i) = 1 + round(a*(ss-1)/(2^code-1));
    end
end

% ����s�õ���Ӧ��Y��X��Y��һ��156*36��0-1����y(ji)��ʾ�ڵ�j��С��i�Ƿ��ڽ���ͬʱ����ͬһ���ڽ�С����������15��
% X��һ��195*36��0-1����x(ji)��ʾ�ڵ�j��С��i�Ƿ��Ѿ����,�ڵ�T��֮��XȡֵΪ1
for k=1+elitism:pop_size
    temp_X = zeros(J,I);
    temp_Y = zeros(J,I);
    ini_pos = s(k,:);  % ini_pos��ȡֵ������[0,112]����֤������С����������ʱ��Ϊ��156��
    for i=1:I
        duration = ini_pos(i)+18+M(i);
        temp_X(duration:J,i) = 1;
        temp_Y(ini_pos(i):duration,i) = 1;
    end
    pop_X{k} = temp_X;
    pop_Y{k} = temp_Y;
end

%-------------------------------������Ӧ��-----------------------------------
%������Ⱥ������Ӧ�ȣ��Բ�ͬ���Ż�Ŀ�꣬�˴���Ҫ��д
%pop_size: ��Ⱥ��С
%chromo_size: Ⱦɫ�峤��

fitness_value=zeros(1,pop_size);
for k=1:pop_size
    temp_X = pop_X{k};
    temp_Y = 15-sum(pop_Y{k},2);
    
    % ����С�����ͼ�����Ӧ��ʵ�ʽ������
    P_real = zeros(size(P));
    B_real = zeros(size(B));
    G_real = zeros(size(G));
    for j=1:J
        for i=1:4
            P_real(j) = P_real(j) + temp_X(j,i)*A(i);
        end
        for i=5:12
            B_real(j) = B_real(j) + temp_X(j,i)*A(i);
        end
        for i=13:36
            G_real(j) = G_real(j) + temp_X(j,i)*A(i);
        end
        fitness_value(k) = fitness_value(k) + 10^8*((min(0,B_real(j)-B(j)))^2 + (min(0,P_real(j)-P(j)))^2 + (min(0,G_real(j)-G(j)))^2 + (min(0,temp_Y(j)))^2);
    end
end

%---------------------�Ը��尴����Ӧ�ȴ�С��������----------------------------
%�Ը��尴��Ӧ�ȴ�С�������򣬲��ұ�����Ѹ���
%pop_size: ��Ⱥ��С
%chromo_size: Ⱦɫ�峤��

[fitness_value,index] = sort(fitness_value);
temp_X = cell(size(pop_X));
temp_Y = cell(size(pop_X));
temp_s = zeros(size(s));
temp_gene = zeros(size(gene));
for k=1:pop_size
    temp_X{k} = pop_X{index(k)};
    temp_Y{k} = pop_Y{index(k)};
    temp_s(k,:) = s(index(k),:);
    temp_gene(k,:) = gene(k,:);
end
pop_X = temp_X;
pop_Y = temp_Y;
s = temp_s;
gene = temp_gene;

% ����Ҫ��fitness_value��С��pop�����Զ�fitness_value��С����������򣬵õ�rank����rank����scale�����õ�scale_value
fitness_table = zeros(1,pop_size);  % �ۻ�
rank = 1:1:pop_size;
temp = 1./(rank.^0.5);
scale_value = temp.*sum(temp)/pop_size;

for i=1:pop_size
    if i==1
        fitness_table(i) = fitness_table(i) + scale_value(i);    
    else
        fitness_table(i) = fitness_table(i-1) + scale_value(i);
    end
end
fitness_avg(generation) = sum(fitness_value)/pop_size;
fitness_value(1)
if fitness_value(1) < best_fitness
    best_fitness = fitness_value(1);
    best_generation = generation;
    best_s = s(1,:);
end

end

%-----------------------------��ӡ�㷨��������-------------------------------
x = 1:1:generation_size;
y = fitness_avg;
figure
hold on
plot(x,y)

%---------------����ʵ�ʽ������ۼ�������Ӧ�������������---------------------
figure
hold on
plot(P,'-b','Linewidth',2)        % ��������
plot(P_real,'-r','Linewidth',2)  % ʵ�ʽ��������
title('Community Type P','FontSize',12,'FontName','Times New Roman');
legend('Demand area','Actual delivery area',2);
legend('boxon');
xlabel('Week','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
ylabel('Area','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
 % ����ͼƬ
% saveas(gcf,'Community Type P.fig');
% saveas(gcf,'Community Type P.eps','psc2');
% saveas(gcf,'Community Type P.jpg');
% close(gcf);

figure
hold on
plot(B,'-b','Linewidth',2)        % ��������
plot(B_real,'-r','Linewidth',2)  % ʵ�ʽ��������
title('Community Type B','FontSize',12,'FontName','Times New Roman');
legend('Demand area','Actual delivery area',2);
legend('boxon');
xlabel('Week','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
ylabel('Area','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
 % ����ͼƬ
% saveas(gcf,'Community Type B.fig');
% saveas(gcf,'Community Type B.eps','psc2');
% saveas(gcf,'Community Type B.jpg');
% close(gcf);

figure
hold on
plot(G,'-b','Linewidth',2)        % ��������
plot(G_real,'-r','Linewidth',2)  % ʵ�ʽ��������
title('Community Type G','FontSize',12,'FontName','Times New Roman');
legend('Demand area','Actual delivery area',2);
legend('boxon');
xlabel('Week','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
ylabel('Area','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
 % ����ͼƬ
% saveas(gcf,'Community Type G.fig');
% saveas(gcf,'Community Type G.eps','psc2');
% saveas(gcf,'Community Type G.jpg');
% close(gcf);

disp "���ſ���ʱ��"
best_s
disp "������Ӧ��"
best_fitness
disp "�õ����Ž���Ĵ���"
best_generation

tic ;%tic2
t2 = clock;
disp(['����������ʱ�䣺',num2str(etime(clock,t1))]);