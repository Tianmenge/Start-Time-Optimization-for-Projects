% ���������ʼ��Ⱥ���޳�������yԼ������Ⱥ
% ���Ƚ���ĳһ��ʵ�����С�������ڸ���֮�����깤ʱ�������С����ǰ���޳�������yԼ������Ⱥ
% Ȼ����15����Ϸ�ʽ������ĳһ��ʵ��������������ڸ���֮ǰ���깤ʱ�������С���Ƴ٣�����Ƴ�֮��ͬʱ����y��xԼ�������Ƴ�

tic;    % ��¼����ʱ��
t1 = clock;

T = 156;   % ����С����������ʱ��
I = 36;    % С����Ŀ
J = 195;   % ���������ܣ�
ss = 156;  % ���ó�ʼ��sȡֵ��Χ������156-26-18
pop_size = 500;           % ��Ⱥ��С
generation_size = 5;    % ��������
combination_size = 6;   % �������ȫ������������С�����ͣ�������С�������Ż����Ż����Ⱥ����2*2^3-1=15����ֻ����A(3,3)=6
fitness_avg = zeros(1,generation_size);  % ÿһ����ƽ��Ŀ�꺯��ֵ��Ŀ�꺯����С
best_fitness = Inf;    % ����Ŀ�꺯��ֵ
best_generation = 0;   % ����Ŀ�꺯��ֵ���ֵĴ���
best_combination = 0;  % ����Ŀ�꺯����Ӧ�����
best_s = zeros(1,I);   % ���Ÿ����Ӧ�Ŀ���ʱ��
best_X = zeros(J:I);
best_Y = zeros(J:I);

% ��ȡ����
data = csvread('data.csv',1,1);
M = data(:,1);
A = data(:,2);
demand = csvread('demand.csv',1);
P = demand(:,2);
B = demand(:,3);
G = demand(:,4);

global s;
global pop_X;
global pop_Y;
% s = zeros(pop_size,I);
% pop_X = cell(1,pop_size);
% pop_Y = cell(1,pop_size);

%----------------------------------��ʼ��-----------------------------------
% �����������ͬһ���ڽ�С����������15����pop_X��pop_Y��s����Ⱥ����Ϊpop_size=36
% pop_X��һ��1*36��Ԫ����ÿһ����һ��195*36��0-1����x(ji)��ʾ�ڵ�j��С��i�Ƿ��Ѿ����
% pop_Y��һ��1*36��Ԫ����ÿһ����һ��195*36��0-1����y(ji)��ʾ�ڵ�j��С��i�Ƿ��ڽ���ͬʱ����ͬһ���ڽ�С����������15��
% s��һ��36*36��0-1����ÿһ�б�ʾһ��36��С���Ŀ���ʱ�䣬����36��
count = 0;
while(1)
    temp_X = zeros(J,I);
    temp_Y = zeros(J,I);
    ini_pos = round(rand(1,36)*(ss-1)+1);  % ini_pos��ȡֵ���䷶Χ��[1,112]����֤����С����������ʱ��Ϊ��156��
    for i = 1:I
        duration = ini_pos(i)+18+M(i);
        temp_X(duration:J,i) = 1;
        temp_Y(ini_pos(i):duration,i) = 1;
    end
    if sum( sum(temp_Y,2) <= 15 ) == J
        count = count + 1;
        pop_X{count} = temp_X;
        pop_Y{count} = temp_Y;
        s(count,:) = ini_pos;
    end
    if count == pop_size
        break;
    end 
end

%--------------------------------������Ӧ��----------------------------------
% ������Ⱥ������Ӧ�ȣ��Բ�ͬ���Ż�Ŀ�꣬�˴���Ҫ��д

fitness_value=zeros(1,pop_size);
for k=1:pop_size
    temp_X = pop_X{k};
    % ����С�����ͼ�����Ӧ��ʵ�ʽ������
    P_real = zeros(size(P));
    B_real = zeros(size(B));
    G_real = zeros(size(G));
    for i=1:4
        P_real = P_real + temp_X(:,i)*A(i);
    end
    for i=5:12
        B_real = B_real + temp_X(:,i)*A(i);
    end
    for i=13:36
        G_real = G_real + temp_X(:,i)*A(i);
    end
    fitness_value(k) = sum((P_real-P).^2) + sum((B_real-B).^2) + sum((G_real-G).^2) + 10^8*(sum((min(0,P_real-P)).^2) + sum((min(0,B_real-B)).^2) + sum((min(0,G_real-G)).^2));
end

%------------------------------������Ѹ���----------------------------------
% �Ը��尴��Ӧ����С��Ϊ��Ѹ���

[temp,index] = sort(fitness_value);
fitness_avg(1) = sum(fitness_value)/pop_size;
fitness_value(1)
if fitness_value(index(1)) < best_fitness
    best_fitness = fitness_value(index(1));
    best_generation = 1;
    best_s = s(index(1),:);
    best_X = pop_X{index(1)};
    best_Y = pop_Y{index(1)};
end

%-------------------------------�����깤ʱ��---------------------------------
% ����ĳһ��j��ʵ�ʽ������С������ģ����ղ�ֵѰ��j��֮���깤ʱ�����j�����С��
% �����깤ʱ��x����j��Ȼ��õ���Ӧ��y��s
while(1)
for k=1:pop_size
    temp_X = pop_X{k};
    temp_Y = pop_Y{k};
    % ����С�����ͼ�����Ӧ��ʵ�ʽ������
    P_real = zeros(size(P));
    B_real = zeros(size(B));
    G_real = zeros(size(G));
    for i=1:4
        P_real = P_real + temp_X(:,i)*A(i);
    end
    for i=5:12
        B_real = B_real + temp_X(:,i)*A(i);
    end
    for i=13:36
        G_real = G_real + temp_X(:,i)*A(i);
    end
    for j=1:J
        if P_real(j) < P(j)
            j;
%             temp_i = find( A(1:4) == (P(j)-P_real(j)) );
            temp_i = 1:1:4;
            temp_s = s(k,temp_i);
            temp_x = temp_s' + M(temp_i) + 18;
            pos = find( (temp_x>j) == 1 );  % ���ղ�ֵѰ��j��֮���깤ʱ�����j�����С��
            min_x = temp_x(pos(1));
            temp = temp_i(pos(1));
            for i=1:length(pos)
                if temp_x(pos(i)) < min_x
                    min_x = temp_x(pos(i));
                    temp = temp_i(pos(i));
                end
            end
            temp_X(j:J,temp) = 1;  % ����temp��С�������ʱ����ǰ
            temp_Y(:,temp) = 0;
            temp_Y((j-M(temp)-18):j,temp) = 1;
            s(k,temp) = j-M(temp)-18;
            P_real = zeros(size(P));
            for i=1:4
                P_real = P_real + temp_X(:,i)*A(i);
            end
        end
    end
    for j=1:J
        if B_real(j) < B(j)
            j;
%             temp_i = 4 + find( A(5:12) == (B(j)-B_real(j)) );
            temp_i = 5:12;
            temp_s = s(k,temp_i);
            temp_x = temp_s' + M(temp_i) + 18;
            pos = find( (temp_x>j) == 1 );  % ���ղ�ֵѰ��j��֮���깤ʱ�����j�����С��
            min_x = temp_x(pos(1));
            temp = temp_i(pos(1));
            for i=1:length(pos)          
                if temp_x(pos(i)) < min_x
                    min_x = temp_x(pos(i));
                    temp = temp_i(pos(i));
                end
            end
            temp_X(j:J,temp) = 1;  % ����temp��С�������ʱ����ǰ
            temp_Y(:,temp) = 0;
            temp_Y((j-M(temp)-18):j,temp) = 1;
            s(k,temp) = j-M(temp)-18;
            B_real = zeros(size(B));
            for i=5:12
                B_real = B_real + temp_X(:,i)*A(i);
            end
        end
    end
    for j=1:J
        if G_real(j) < G(j)
            j;
%             temp_i = 12 + find( A(13:36) >= (G(j)-G_real(j)) );
            temp_i = 13:1:36;
            temp_s = s(k,temp_i);
            temp_x = temp_s' + M(temp_i) + 18;
            pos = find( (temp_x>j) == 1 );  % ���ղ�ֵѰ��j��֮���깤ʱ�����j�����С��
            min_x = temp_x(pos(1));
            temp = temp_i(pos(1));
            for i=1:length(pos)
                if temp_x(pos(i)) < min_x
                    min_x = temp_x(pos(i));
                    temp = temp_i(pos(i));
                end
            end
            temp_X(j:J,temp) = 1;  % ����temp��С�������ʱ����ǰ
            temp_Y(:,temp) = 0;
            temp_Y((j-M(temp)-18):j,temp) = 1;
            s(k,temp) = j-M(temp)-18;
            G_real = zeros(size(G));
            for i=13:36
                G_real = G_real + temp_X(:,i)*A(i);
            end
        end
    end
    pop_X{k} = temp_X;
    pop_Y{k} = temp_Y;
end
if (sum(P_real>=P) == J) || (sum(B_real>=B) == J) || (sum(G_real>=G) == J)
    break;
end
end

% �޳�������yԼ������Ⱥ��һ�㵱ssС��156ʱ��ʵ�ʽ�����������������
% ���Ե����깤ʱ��֮�����Ⱥ����ʼ��Ⱥ��࣬Ҳ����˵count����Ϊ0
temp_X = {};
temp_Y = {};
temp_s = [];
count = 0;
for k=1:pop_size
    temp1 = pop_Y{k};
    temp2 = pop_X{k};
    if sum( sum(temp1,2) <= 15 ) == J
        count = count + 1;
        temp_X{count} = pop_X{k};
        temp_Y{count} = pop_Y{k};
        temp_s(count,:) = s(k,:);
    end
end
s = temp_s;
pop_X = temp_X;
pop_Y = temp_Y;
count
pop_size = count;

for generation=2:generation_size
    generation
%-------------------------------��������ʱ��---------------------------------
% ����ĳһ��j��ʵ�ʽ��������������ģ����ղ�ֵѰ��j��֮ǰ�깤ʱ�����j��Զ��С��
% �����깤ʱ��x����j��Ȼ��õ���Ӧ��y�����y����Լ�������������x,s�Լ�ʵ�ʽ������
% %{
type=['P','B','G'];
for k=1:pop_size
%     [pop_X, pop_Y, s] = computPreal(k, pop_X, pop_Y, s, P, M, A, J);
%     [pop_X, pop_Y, s] = computBreal(k, pop_X, pop_Y, s, B, M, A, J);
%     [pop_X, pop_Y, s] = computGreal(k, pop_X, pop_Y, s, G, M, A, J);
%     %{
    com_X = cell(1, combination_size);
    com_Y = cell(1, combination_size);
    com_s = zeros(combination_size, I);
    combination = zeros(1, I);
%     n = cell(1,15);
%     count = 0;
    best_com_fitness = Inf;
    temp = perms(type);
    for c=1:size(temp,1)
        for i=1:3
            temp(c,i);
            if temp(c,i)=='P'
                [pop_X, pop_Y, s] = computPreal(k, pop_X, pop_Y, s, P, M, A, J);
            elseif temp1(c,i)=='B'
                [pop_X, pop_Y, s] = computBreal(k, pop_X, pop_Y, s, B, M, A, J);
            else
                [pop_X, pop_Y, s] = computGreal(k, pop_X, pop_Y, s, G, M, A, J);
            end
        end
        com_X{c} = pop_X{k};
        com_Y{c} = pop_Y{k};
        com_s(c,:) = s(k,:);
    end
    %{
    for t=1:3
        temp = combnk(type,t);
        for c=1:size(temp,1)
            temp1=perms(temp(c,:));
            for j=1:size(temp1,1)
%                 n{count+j}=temp1(j,:);
                temp1(j,:);
                for i=1:t
                    temp1(j,i);
                    if temp1(j,i)=='P'
                        [pop_X, pop_Y, s] = computPreal(k, pop_X, pop_Y, s, P, M, A, J);
                    elseif temp1(j,i)=='B'
                        [pop_X, pop_Y, s] = computBreal(k, pop_X, pop_Y, s, B, M, A, J);
                    else
                        [pop_X, pop_Y, s] = computGreal(k, pop_X, pop_Y, s, G, M, A, J);
                    end
                end
                com_X{count+j} = pop_X{k};
                com_Y{count+j} = pop_Y{k};
                com_s(count+j,:) = s(k,:);
            end
            count = count + size(temp1,1);
        end
    end
    %}
    [best_com_fitness, best_combination, best_s, best_X, best_Y] = com_fitness(combination_size, pop_X, pop_Y, s, P, B, G, A, best_com_fitness, best_combination, best_s, best_X, best_Y);
    pop_X{k} = best_X;
    pop_Y{k} = best_Y;
    s(k,:) = best_s;
    combination(k) = best_combination;
    %}
end
[fitness_avg_pop, best_fitness, best_generation, best_combination, best_s, best_X, best_Y] = fitness(pop_size, pop_X, pop_Y, s, P, B, G, A, best_fitness, generation, combination, best_generation, best_combination, best_s, best_X, best_Y);
fitness_avg(generation) = fitness_avg_pop;
end

%-----------------------------��ӡ�㷨��������-------------------------------
x = 2:1:generation_size;
y = fitness_avg(2:generation_size);
figure
hold on
plot(x,y)

% �������ſ���ʱ�䣬����С�����ͼ�����Ӧ��ʵ�ʽ������
P_real = zeros(size(P));
B_real = zeros(size(B));
G_real = zeros(size(G));
for i=1:4
    P_real = P_real + best_X(:,i)*A(i);
end
for i=5:12
    B_real = B_real + best_X(:,i)*A(i);
end
for i=13:36
    G_real = G_real + best_X(:,i)*A(i);
end

%---------------����ʵ�ʽ������ۼ�������Ӧ�������������---------------------
figure
hold on
plot(P,'-b','Linewidth',2)        % ��������
plot(P_real,'-r','Linewidth',2)   % ʵ�ʽ��������
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
plot(B_real,'-r','Linewidth',2)   % ʵ�ʽ��������
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
plot(G_real,'-r','Linewidth',2)   % ʵ�ʽ��������
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
disp "�õ����Ž�������"
best_combination

tic ;  % tic2
t2 = clock;
disp(['����������ʱ�䣺',num2str(etime(clock,t1))]);
