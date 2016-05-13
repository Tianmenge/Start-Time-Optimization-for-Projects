% 遗传算法
% 将约束优化转化为无约束优化进行求解
% 将不同小区的开始时间S作为一个染色体，进行初始化编码，解码之后得到相应的x值和y值，然后判断是否满足约束条件
tic;    % 记录运行时间
t1 = clock;

T = 156;   % 所有小区的最迟完成时间
I = 36;    % 小区数目
J = 195;   % 总期数（周）
code = 8;  % 用8位二进制数对s进行编码
ss = 120;  % s取值范围的上限
pop_size = 36;              % 种群大小
chromo_size = I*code;       % 染色体大小
generation_size = 10;     % 迭代次数
cross_rate = 0.6;          % 交叉概率
mutate_rate = 0.01;        % 变异概率
elitism = 2;
fitness_avg = zeros(generation_size,1);  % 每一代的最小目标函数值
best_fitness = Inf;    % 最优目标函数值
best_generation = 0;   % 最优目标函数值出现的代数
best_individual = zeros(1,chromo_size);    % 最优个体
best_s = zeros(1,I);  % 最优个体对应的开工时间

% 读取数据
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

% 首先对初始化的第1代计算适应度、排序等，然后进行迭代
% 依次进行选择、交叉、变异等操作，并计算适应度、排序等

%----------------------------------初始化-----------------------------------
% 初始化种群基因型，对不同小区的开工时间采用8位二进制数进行编码

for k=1:pop_size
    for i=1:chromo_size
        gene(k,i) = round(rand());
    end
end
   
%-----------------------------------解码------------------------------------
% 对初始化的基因型进行解码，得到相应的不同小区的开工时间以及相应的x和y

for k=1:pop_size
    for i=1:I
        b = gene(k,(code*i-code+1):(code*i));
        a = bin2dec(num2str(b));
        s(k,i) = 1 + round(a*(ss-1)/(2^code-1));
    end
end

% 根据s得到相应的Y和X，Y是一个156*36的0-1矩阵，y(ji)表示在第j期小区i是否在建，同时满足同一期在建小区数不超过15个
% X是一个195*36的0-1矩阵，x(ji)表示在第j期小区i是否已经完成,在第T期之后X取值为1
for k=1:pop_size
    temp_X = zeros(J,I);
    temp_Y = zeros(J,I);
    ini_pos = s(k,:);  % ini_pos的取值区间是[0,112]，保证了所有小区的最迟完成时间为第156期
    for i=1:I
        duration = ini_pos(i)+18+M(i);
        temp_X(duration:J,i) = 1;
        temp_Y(ini_pos(i):duration,i) = 1;
    end
    pop_X{k} = temp_X;
    pop_Y{k} = temp_Y;
end

% % {
%--------------------------------计算适应度----------------------------------
%计算种群个体适应度，对不同的优化目标，此处需要改写
%pop_size: 种群大小
%chromo_size: 染色体长度

fitness_value=zeros(1,pop_size);
for k=1:pop_size
    temp_X = pop_X{k};
    temp_Y = 15-sum(pop_Y{k},2);
    
    % 按照小区类型计算相应的实际交付面积
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

%----------------------对个体按照适应度大小进行排序---------------------------
%对个体按适应度大小进行排序，并且保存最佳个体
%pop_size: 种群大小
%chromo_size: 染色体长度

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

% 由于要求fitness_value最小的pop，所以对fitness_value从小到大进行排序，得到rank，对rank进行scale处理，得到scale_value
fitness_table = zeros(1,pop_size);  % 累积
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

%========================迭代次数 generation_size===========================
for generation=2:generation_size
    generation

%---------------------------------选择操作----------------------------------
% 首先剔除不满足约束条件的种群，对剩下的种群按照适应度进行排序，组成新的36个种群
temp_X = {};
temp_Y = {};
temp_s = [];
temp_gene = [];
count = 0;
temp_value = [];
% 剔除不满足约束条件的种群
for k=1:pop_size
    temp1 = pop_Y{k};
    temp2 = pop_X{k};
    if sum( sum(temp1,2) <= 15 ) == J
        % 按照小区类型计算相应的实际交付面积
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
% 按照适应度大小进行排序
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
% 如果count小于pop_size，则将temp_X等循环写入pop_X
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

% 轮盘赌选择操作，将最优的前两个个体直接遗传到下一代
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
    % 使用二分法选择个体
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

%----------------------------------交叉操作----------------------------------
% 对相邻两个种群进行单点交叉操作，交叉概率为cross_rate
% 任选需要交叉的两个种群，选择任意交叉点，交叉点之后的基因全部互换

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

%------------------------------交叉之后解码----------------------------------
% 对初始化的基因型进行解码，得到相应的不同小区的开工时间以及相应的x和y

for k=1+elitism:pop_size
    for i=1:I
        b = gene(k,(code*i-code+1):(code*i));
        a = bin2dec(num2str(b));
        s(k,i) = 1 + round(a*(ss-1)/(2^code-1));
    end
end

% 根据s得到相应的Y和X，Y是一个156*36的0-1矩阵，y(ji)表示在第j期小区i是否在建，同时满足同一期在建小区数不超过15个
% X是一个195*36的0-1矩阵，x(ji)表示在第j期小区i是否已经完成,在第T期之后X取值为1
for k=1+elitism:pop_size
    temp_X = zeros(J,I);
    temp_Y = zeros(J,I);
    ini_pos = s(k,:);  % ini_pos的取值区间是[0,112]，保证了所有小区的最迟完成时间为第156期
    for i=1:I
        duration = ini_pos(i)+18+M(i);
        temp_X(duration:J,i) = 1;
        temp_Y(ini_pos(i):duration,i) = 1;
    end
    pop_X{k} = temp_X;
    pop_Y{k} = temp_Y;
end

%---------------------------------变异操作-----------------------------------
% 单点变异操作，变异概率mutate_rate，任选基因中的一个点进行变异，0-1互换

for i=1+elitism:pop_size
    if rand < mutate_rate
        mutate_pos = 1 + round(rand() * (chromo_size-1));
        gene(i,mutate_pos) = 1 - gene(i, mutate_pos);
    end
end

%------------------------------变异之后解码----------------------------------
% 对初始化的基因型进行解码，得到相应的不同小区的开工时间以及相应的x和y

for k=1+elitism:pop_size
    for i=1:I
        b = gene(k,(code*i-code+1):(code*i));
        a = bin2dec(num2str(b));
        s(k,i) = 1 + round(a*(ss-1)/(2^code-1));
    end
end

% 根据s得到相应的Y和X，Y是一个156*36的0-1矩阵，y(ji)表示在第j期小区i是否在建，同时满足同一期在建小区数不超过15个
% X是一个195*36的0-1矩阵，x(ji)表示在第j期小区i是否已经完成,在第T期之后X取值为1
for k=1+elitism:pop_size
    temp_X = zeros(J,I);
    temp_Y = zeros(J,I);
    ini_pos = s(k,:);  % ini_pos的取值区间是[0,112]，保证了所有小区的最迟完成时间为第156期
    for i=1:I
        duration = ini_pos(i)+18+M(i);
        temp_X(duration:J,i) = 1;
        temp_Y(ini_pos(i):duration,i) = 1;
    end
    pop_X{k} = temp_X;
    pop_Y{k} = temp_Y;
end

%-------------------------------计算适应度-----------------------------------
%计算种群个体适应度，对不同的优化目标，此处需要改写
%pop_size: 种群大小
%chromo_size: 染色体长度

fitness_value=zeros(1,pop_size);
for k=1:pop_size
    temp_X = pop_X{k};
    temp_Y = 15-sum(pop_Y{k},2);
    
    % 按照小区类型计算相应的实际交付面积
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

%---------------------对个体按照适应度大小进行排序----------------------------
%对个体按适应度大小进行排序，并且保存最佳个体
%pop_size: 种群大小
%chromo_size: 染色体长度

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

% 由于要求fitness_value最小的pop，所以对fitness_value从小到大进行排序，得到rank，对rank进行scale处理，得到scale_value
fitness_table = zeros(1,pop_size);  % 累积
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

%-----------------------------打印算法迭代过程-------------------------------
x = 1:1:generation_size;
y = fitness_avg;
figure
hold on
plot(x,y)

%---------------绘制实际交付的累计面积与对应的需求面积曲线---------------------
figure
hold on
plot(P,'-b','Linewidth',2)        % 需求的面积
plot(P_real,'-r','Linewidth',2)  % 实际交付的面积
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
plot(B_real,'-r','Linewidth',2)  % 实际交付的面积
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
plot(G_real,'-r','Linewidth',2)  % 实际交付的面积
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

tic ;%tic2
t2 = clock;
disp(['程序总运行时间：',num2str(etime(clock,t1))]);