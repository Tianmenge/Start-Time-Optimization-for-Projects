%--------------------------------计算适应度----------------------------------
% 计算种群个体适应度，对不同的优化目标，此处需要改写
function [fitness_avg_pop, best_fitness, best_generation, best_combination, best_s, best_X, best_Y] = fitness(pop_size, pop_X, pop_Y, s, P, B, G, A, best_fitness, generation, combination, best_generation, best_combination, best_s, best_X, best_Y)

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
fitness_avg_pop = sum(fitness_value)/pop_size;
if fitness_value(index(1)) < best_fitness
    best_fitness = fitness_value(index(1));
    best_generation = generation;
    best_combination = combination(index(1));
    best_s = s(index(1),:);
    best_X = pop_X{index(1)};
    best_Y = pop_Y{index(1)};
end
