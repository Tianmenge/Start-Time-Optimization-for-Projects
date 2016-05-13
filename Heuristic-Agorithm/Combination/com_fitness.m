% 不同组合对应的适应度，同一个种群中选择适应度最小的组合
function [best_com_fitness, best_combination, best_s, best_X, best_Y] = com_fitness(combination_size, pop_X, pop_Y, s, P, B, G, A, best_com_fitness, best_combination, best_s, best_X, best_Y)
% , best_generation
fitness_value=zeros(1,combination_size);
for k=1:combination_size
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

[temp,index] = sort(fitness_value);
% fitness_avg(generation) = sum(fitness_value)/pop_size;
% best_combination;
if fitness_value(index(1)) < best_com_fitness
    best_com_fitness = fitness_value(index(1));
    best_combination = index(1);
    best_s = s(index(1),:);
    best_X = pop_X{index(1)};
    best_Y = pop_Y{index(1)};
end

clear k;
clear i;
clear temp_X;
clear fitness_value;
clear index;
clear temp;
