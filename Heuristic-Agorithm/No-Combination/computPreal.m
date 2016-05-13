% ����P����С����ʵ�ʽ������
function [pop_X, pop_Y, s] = computPreal(k, pop_X, pop_Y, s, P, M, A, J)
    temp_X = pop_X{k};
    % ����С�����ͼ�����Ӧ��ʵ�ʽ������
    P_real = zeros(size(P));
    for i=1:4
        P_real = P_real + temp_X(:,i)*A(i);
    end

    for j=1:J
        if P_real(j) > P(j)
            temp_X = pop_X{k};
            temp_Y = pop_Y{k};
            pos1 = j;
            for temp_j = j:J
                if P_real(pos1) - P(pos1) ~= P_real(temp_j) - P(temp_j)
%                     temp_i = find( A(1:4) == (P_real(pos1)-P(pos1)) );
                    temp_i = 1:1:4;
                    temp_s = s(k,temp_i);
                    temp_x = temp_s' + M(temp_i) + 18;
                    pos = find( (temp_x<temp_j) == 1 );  % ���ղ�ֵѰ��j��֮ǰ�깤ʱ�����j�����С��
                    for i=1:length(pos)
                        temp = temp_i(pos(i));
                        temp_Y(:,temp) = 0;
                        temp_Y((temp_j-M(temp)-18):temp_j,temp) = 1;
                        temp_X(:,temp) = 0;
                        temp_X(temp_j:J,temp) = 1;  % ����temp��С�������ʱ���Ƴٵ�j
                        P_real = zeros(size(P));
                        for i=13:36
                            P_real = P_real + temp_X(:,i)*A(i);
                        end
                        % �ж��Ƴ�֮���y�Ƿ�����Ҫ���Լ�ʵ�ʽ�������Ƿ��������
                        if (sum( sum(temp_Y,2) <= 15 ) == J) && (sum(P_real >= P) == J)
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
                    P_real = zeros(size(P));
                    for i=1:4
                        P_real = P_real + temp_X(:,i)*A(i);
                    end
                    j = temp_j;
                    break;
                end
            end
            % �ж��Ƴ�֮���y�Ƿ�����Ҫ���Լ�ʵ�ʽ�������Ƿ��������
            if (sum( sum(temp_Y,2) <= 15 ) == J) && (sum(P_real >= P) == J)
                s(k,temp) = temp_j-M(temp)-18;
                pop_X{k} = temp_X;
                pop_Y{k} = temp_Y;
            else
                continue;
            end
        end
    end
%}
