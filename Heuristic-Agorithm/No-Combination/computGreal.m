% ����G����С����ʵ�ʽ������
function [pop_X, pop_Y, s] = computGreal(k, pop_X, pop_Y, s, G, M, A, J)
    temp_X = pop_X{k};
    % ����С�����ͼ�����Ӧ��ʵ�ʽ������
    G_real = zeros(size(G));
    for i=13:36
        G_real = G_real + temp_X(:,i)*A(i);
    end
    
     for j=1:J
         if G_real(j) > G(j)
            temp_X = pop_X{k};
            temp_Y = pop_Y{k};
            pos1 = j;
            for temp_j = j:J
                if G_real(pos1) - G(pos1) ~= G_real(temp_j) - G(temp_j)
%                     temp_i = 12 + find( A(13:36) <= (G_real(j)-G(j)) );
                    temp_i = 13:1:36;
                    temp_s = s(k,temp_i);
                    temp_x = temp_s' + M(temp_i) + 18;
                    pos = find( (temp_x<temp_j) == 1 );  % ���ղ�ֵѰ��j��֮ǰ�깤ʱ�����j�����С��
                    for i=1:length(pos)
                        temp = temp_i(pos(i));
                        temp_Y(:,temp) = 0;
                        temp_Y((temp_j-M(temp)-18):temp_j,temp) = 1;
                        temp_X(:,temp) = 0;
                        temp_X(temp_j:J,temp) = 1;  % ����temp��С�������ʱ���Ƴٵ�j
                        G_real = zeros(size(G));
                        for i=13:36
                            G_real = G_real + temp_X(:,i)*A(i);
                        end
                        % �ж��Ƴ�֮���y�Ƿ�����Ҫ���Լ�ʵ�ʽ�������Ƿ��������
                        if (sum( sum(temp_Y,2) <= 15 ) == J) && (sum(G_real >= G) == J)
                            s(k,temp) = temp_j-M(temp)-18;
                            pop_X{k} = temp_X;
                            pop_Y{k} = temp_Y;
                        else
                            temp_X = pop_X{k};
                            temp_Y = pop_Y{k};
                        end
                    end
                    j = temp_j;
                    break;
                end
            end
         end
     end
%{
                    max_x = temp_x(pos(1));
                    temp = temp_i(pos(1));
                    for i=1:length(pos)
                        if temp_x(pos(i)) > max_x
                            max_x = temp_x(pos(i));
                            temp = temp_i(pos(i));
                        end
                    end
                    temp_Y(:,temp) = 0;
                    temp_Y((temp_j-M(temp)-18):temp_j,temp) = 1;
                    temp_X(:,temp) = 0;
                    temp_X(temp_j:J,temp) = 1;  % ����temp��С�������ʱ���Ƴٵ�j
                    G_real = zeros(size(G));
                    for i=13:36
                        G_real = G_real + temp_X(:,i)*A(i);
                    end
                    j = temp_j;
                    break;
                end
            end
            % �ж��Ƴ�֮���y�Ƿ�����Ҫ���Լ�ʵ�ʽ�������Ƿ��������
            if (sum( sum(temp_Y,2) <= 15 ) == J) && (sum(G_real >= G) == J)
                s(k,temp) = temp_j-M(temp)-18;
                pop_X{k} = temp_X;
                pop_Y{k} = temp_Y;
            else
                continue;
            end
         end
     end
%}
                