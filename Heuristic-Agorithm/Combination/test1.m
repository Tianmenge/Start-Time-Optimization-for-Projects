% 随机产生初始种群，剔除不满足y约束的种群
% 首先将在某一期实际面积小于需求，在该期之后且完工时间最近的小区提前，剔除不满足y约束的种群
% 然后按照15中组合方式，将在某一期实际面积大于需求，在该期之前且完工时间最近的小区推迟，如果推迟之后不同时满足y和x约束，则不推迟

tic;    % 记录运行时间
t1 = clock;

T = 156;   % 所有小区的最迟完成时间
I = 36;    % 小区数目
J = 195;   % 总期数（周）
ss = 156;  % 设置初始的s取值范围的上限156-26-18
pop_size = 500;           % 种群大小
generation_size = 5;    % 迭代次数
combination_size = 6;   % 所有组合全排列数，三种小区类型，对那种小区进行优化，优化的先后次序，2*2^3-1=15，若只考虑A(3,3)=6
fitness_avg = zeros(1,generation_size);  % 每一代的平均目标函数值，目标函数最小
best_fitness = Inf;    % 最优目标函数值
best_generation = 0;   % 最优目标函数值出现的代数
best_combination = 0;  % 最优目标函数对应的组合
best_s = zeros(1,I);   % 最优个体对应的开工时间
best_X = zeros(J:I);
best_Y = zeros(J:I);

% 读取数据
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

%----------------------------------初始化-----------------------------------
% 随机产生满足同一期在建小区数不超过15个的pop_X、pop_Y和s，种群数量为pop_size=36
% pop_X是一个1*36的元胞，每一项是一个195*36的0-1矩阵，x(ji)表示在第j期小区i是否已经完成
% pop_Y是一个1*36的元胞，每一项是一个195*36的0-1矩阵，y(ji)表示在第j期小区i是否在建，同时满足同一期在建小区数不超过15个
% s是一个36*36的0-1矩阵，每一行表示一组36个小区的开工时间，共有36组
count = 0;
while(1)
    temp_X = zeros(J,I);
    temp_Y = zeros(J,I);
    ini_pos = round(rand(1,36)*(ss-1)+1);  % ini_pos的取值区间范围是[1,112]，保证所有小区的最迟完成时间为第156期
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

%--------------------------------计算适应度----------------------------------
% 计算种群个体适应度，对不同的优化目标，此处需要改写

fitness_value=zeros(1,pop_size);
for k=1:pop_size
    temp_X = pop_X{k};
    % 按照小区类型计算相应的实际交付面积
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

%------------------------------保存最佳个体----------------------------------
% 对个体按适应度最小的为最佳个体

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

%-------------------------------调整完工时间---------------------------------
% 对于某一期j内实际建筑面积小于需求的，按照差值寻找j期之后完工时间距离j最近的小区
% 令其完工时间x等于j，然后得到相应的y和s
while(1)
for k=1:pop_size
    temp_X = pop_X{k};
    temp_Y = pop_Y{k};
    % 按照小区类型计算相应的实际交付面积
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
            pos = find( (temp_x>j) == 1 );  % 按照差值寻找j期之后完工时间距离j最近的小区
            min_x = temp_x(pos(1));
            temp = temp_i(pos(1));
            for i=1:length(pos)
                if temp_x(pos(i)) < min_x
                    min_x = temp_x(pos(i));
                    temp = temp_i(pos(i));
                end
            end
            temp_X(j:J,temp) = 1;  % 将第temp个小区的完成时间提前
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
            pos = find( (temp_x>j) == 1 );  % 按照差值寻找j期之后完工时间距离j最近的小区
            min_x = temp_x(pos(1));
            temp = temp_i(pos(1));
            for i=1:length(pos)          
                if temp_x(pos(i)) < min_x
                    min_x = temp_x(pos(i));
                    temp = temp_i(pos(i));
                end
            end
            temp_X(j:J,temp) = 1;  % 将第temp个小区的完成时间提前
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
            pos = find( (temp_x>j) == 1 );  % 按照差值寻找j期之后完工时间距离j最近的小区
            min_x = temp_x(pos(1));
            temp = temp_i(pos(1));
            for i=1:length(pos)
                if temp_x(pos(i)) < min_x
                    min_x = temp_x(pos(i));
                    temp = temp_i(pos(i));
                end
            end
            temp_X(j:J,temp) = 1;  % 将第temp个小区的完成时间提前
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

% 剔除不满足y约束的种群，一般当ss小于156时，实际建筑面积都会大于需求
% 所以调整完工时间之后的种群跟初始种群差不多，也就是说count不会为0
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
%-------------------------------调整开工时间---------------------------------
% 对于某一期j内实际建筑面积大于需求的，按照差值寻找j期之前完工时间距离j最远的小区
% 令其完工时间x等于j，然后得到相应的y，如果y满足约束条件，则更新x,s以及实际建筑面积
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

%-----------------------------打印算法迭代过程-------------------------------
x = 2:1:generation_size;
y = fitness_avg(2:generation_size);
figure
hold on
plot(x,y)

% 对于最优开工时间，按照小区类型计算相应的实际交付面积
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

%---------------绘制实际交付的累计面积与对应的需求面积曲线---------------------
figure
hold on
plot(P,'-b','Linewidth',2)        % 需求的面积
plot(P_real,'-r','Linewidth',2)   % 实际交付的面积
title('Community Type P','FontSize',12,'FontName','Times New Roman');
legend('Demand area','Actual delivery area',2);
legend('boxon');
xlabel('Week','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
ylabel('Area','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
% 保存图片
% saveas(gcf,'Community Type P.fig');
% saveas(gcf,'Community Type P.eps','psc2');
% saveas(gcf,'Community Type P.jpg');
% close(gcf);

figure
hold on
plot(B,'-b','Linewidth',2)        % 需求的面积
plot(B_real,'-r','Linewidth',2)   % 实际交付的面积
title('Community Type B','FontSize',12,'FontName','Times New Roman');
legend('Demand area','Actual delivery area',2);
legend('boxon');
xlabel('Week','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
ylabel('Area','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
% 保存图片
% saveas(gcf,'Community Type B.fig');
% saveas(gcf,'Community Type B.eps','psc2');
% saveas(gcf,'Community Type B.jpg');
% close(gcf);

figure
hold on
plot(G,'-b','Linewidth',2)        % 需求的面积
plot(G_real,'-r','Linewidth',2)   % 实际交付的面积
title('Community Type G','FontSize',12,'FontName','Times New Roman');
legend('Demand area','Actual delivery area',2);
legend('boxon');
xlabel('Week','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
ylabel('Area','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
% 保存图片
% saveas(gcf,'Community Type G.fig');
% saveas(gcf,'Community Type G.eps','psc2');
% saveas(gcf,'Community Type G.jpg');
% close(gcf);

disp "最优开工时间"
best_s
disp "最优适应度"
best_fitness
disp "得到最优结果的代数"
best_generation
disp "得到最优结果的组合"
best_combination

tic ;  % tic2
t2 = clock;
disp(['程序总运行时间：',num2str(etime(clock,t1))]);
