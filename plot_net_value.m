function plot_net_value(Loc, A, value)

stem3(Loc(:, 1), Loc(:, 2), value);
hold on;
plot_net(Loc, A);
box off
grid off
axis off