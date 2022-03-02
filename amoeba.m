function [Q,P]= amoeba(L,C,A,a,OL)

n = length(L{1}(:,1));
D = 0.5*ones(n);
D(L{1} == 0) = 0;
L = cello2inf(L);
% sourceNode=1;
% sinkNode=n;
temp_Q = ones(n);
Q = zeros(n);
numite = 0;
while max(max(abs(temp_Q-Q)))>0.1 && numite < 200
    temp_Q=Q;
    B=consdivcell(D, L);
    B = cell4B(B,n);
    % Sole the pressure in each points
    P = cell4P(B, A, n);

%     tempP=repmat(P,1,n)-repmat(P',n,1);
%     Q=(D./L).*tempP;
    Q = zeros(n);
    [i,j] = find(triu(D));
    for ite = 1: length(i)
        fuz1 = P(i(ite), :) ./ [L{3}(i(ite), j(ite)), L{2}(i(ite), j(ite)), L{1}(i(ite), j(ite))];
        fuz2 = P(j(ite), :) ./ [L{3}(i(ite), j(ite)), L{2}(i(ite), j(ite)), L{1}(i(ite), j(ite))];
        temp = sum((fuz1 - fuz2).^2) + (fuz1(2) - fuz2(2))^2;
        temp = temp + (fuz1(1) - fuz2(1))*(fuz1(2) - fuz2(2)) + (fuz1(3) - fuz2(3))*(fuz1(2) - fuz2(2));
        Dis = sqrt(abs(temp/6));
        Q(i(ite), j(ite)) = D(i(ite), j(ite)) * Dis;
    end
    Q = Q + Q';
    D=(Q+D)./2;
    L = fuzzplusdiv(UE(OL,Q,C,a) , L);
    L = cello2inf(L);
    numite = numite +1;
    if rem(numite, 10) == 0
        disp(numite)
        disp(max(max(abs(temp_Q-Q))))
    end
end
P = 1/6* (P(:,1) + 4 * P(:,2) + P(:,3));
end

% subfunctions
function L = UE(OL,Q,C,a)
    L1 = (1-a)*0.15*(Q./C).^(4*(1-a));
    L2 = 0.15*(Q./C).^4;
    L3 = (1+a)*0.15*(Q./C).^(4*(1+a));
    L = {OL+L1, OL+L2, OL+L3};
end

function L = cello2inf(L)
    for ite = 1: length(L)
        L{ite}(L{ite} == 0) = inf; 
    end
end

function B = consdivcell(D, L)
    B = cell(size(L));
    for ite = 1 : length(L)
        B{ite} = D./L{ite};
    end
    temp = B{1};
    B{1} = B{3};
    B{3} = temp;
end

function B = cell4B(B,n)
    for ite = 1: length(B)
        B{ite} = B{ite}-diag(sum(B{ite}));
        B{ite}(:, n) = [];
    end
end

function P = cell4P(B, A,n)
    P = cell(size(B));
   for ite =  1 : length(B)
       P{ite} = B{ite} \ A;
   end
   P = [P{3},P{2},P{1}];
   P(n, :) = 0; 
end

function C = fuzzplusdiv(A, B)
    C = {(A{1}+B{1})/2, (A{2} + B{2})/2, (A{3} + B{3})/2};
end

function A = disfuzz(B)
    A = 1/6*(B{1} + 4*B{2} + B{3});
end