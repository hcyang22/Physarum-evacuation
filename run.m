%% main script
clear; clc
load("./data/chicago.mat")
% net: link info
% tag: tag info of net.mat
% node: node id and cordinate
% exnode: exit node id vector
links = net(:,[1,2]);
n = length(node(:,1)) + 1;
Q = zeros(n);
OL = zeros(n);
C = zeros(n);
for ite = 1 : length(links(:,1))
    Q(links(ite,1),links(ite,2)) = 0.01;
    OL(links(ite,1),links(ite,2)) = net(ite,4);
    C(links(ite,1),links(ite,2)) = net(ite,3);
end
Q(n, exnode) = 0.01;
C(n, exnode) = inf;
OL(n, exnode) = 0.01;
Q(exnode, n) = 0.01;
C(exnode, n) = inf;
OL(exnode, n) = 0.01;
% Q = Q + Q';
% OL = OL + OL';
% C = C + C';
a = 0.2;
OL(OL == 0) = inf;
C = C;
C(C==0) = inf;
L = UE(OL,Q,C,a);
A = - ones(n,1);
A(n) = n - 1;

[Q,P] = amoeba(L, C, A, a, OL);

%% plot script
%  Note: this script use licensed matlab toolbox: davinci. If you encounter
%  license problem, please install that toolbox and get a license. You may
%  also use other methods to visualize the result. 
[x,y] = find(triu(Q));
links = [x,y];
expan0 = 2000;
maxexpan = 4000;
hold on;
plot(node(:, 2), node(:, 3),'.r');
maxl = 0;
for ite = 1 : 11
    plot(node(exnode(ite), 2), node(exnode(ite), 3),'*r');
end

for ite = 1 : length(links(:,1))
    if Q(links(ite,1), links(ite,2)) > 0.1
        if max(links(ite,:)) == 934
            % skip the virtual sink
            continue
        end
        dis = sqrt((node(links(ite,1),2)-node(links(ite,2),2))^2 + (node(links(ite,1),3)-node(links(ite,2),3))^2);
        if dis < 4*expan0
            expan = min([floor(dis/4), maxexpan]);
        else
            expan = expan0;
        end
%         disp(expan)
        if P(links(ite,1)) > P(links(ite,2))
            davinci('arrow','X',[node(links(ite,1),2), node(links(ite,2),2)],'Y',[node(links(ite,1),3), node(links(ite,2),3)],'Head.Length', 3*expan, 'Head.Sweep', 2*expan, 'Head.Width',4*expan)
        else
            davinci('arrow','X',[node(links(ite,2),2), node(links(ite,1),2)],'Y',[node(links(ite,2),3), node(links(ite,1),3)],'Head.Length', 3*expan, 'Head.Sweep', 2*expan, 'Head.Width',4*expan)
        end
    end
end
plot(1e5 * [6.5, 7], 1e6 * 2.1 * [1, 1], 'k-', "LineWidth",2)


% axis setting
axis equal
axis off
% plot(node(918, 2), node(918, 3),'*r');
hold off;

%% subfunctions
function L = UE(OL,Q,C,a)
    L1 = (1-a)*0.15*(Q./C).^(4*(1-a));
    L2 = 0.15*(Q./C).^4;
    L3 = (1+a)*0.15*(Q./C).^(4*(1+a));
    L = {OL+L1, OL+L2, OL+L3};
end