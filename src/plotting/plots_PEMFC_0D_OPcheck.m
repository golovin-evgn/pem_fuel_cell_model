function plots_PEMFC_0D_OPcheck(x_vec,t,u_traj_pp,voli,p)
%%%
%%% plots of PEMFC model %%%%%%%%%%%%%%
%%%


% style of lines
line_style = '-';
line_style2 = '--';

% --------------------------------------------------------------------------------------  
% from vector to structure
x = p.state2struct(x_vec.');
u = p.inputs2struct(u_traj_pp.data.');

% --------------------------------------------------------------------------------------   
% variables
% relative humidities
a_H2O_c = (p.R * x.T_c .* x.c_H2O_c)./p.p_sat(x.T_c);
a_H2O_a = (p.R * x.T_a .* x.c_H2O_a)./p.p_sat(x.T_a);
% concentrations of water on sides of membrane
c_H2O_ca = x.xi_H2O_ca .* x.p_a ./ (p.R * x.T_s);
c_H2O_cc = x.xi_H2O_cc .* x.p_c ./ (p.R * x.T_s);
% relative humidities on the sides of membrane
a_H2O_ca = (p.R * x.T_s .* c_H2O_ca)./p.p_sat(x.T_s);
a_H2O_cc = (p.R * x.T_s .* c_H2O_cc)./p.p_sat(x.T_s); 
% conductivity of membrane
kappa = p.kappa(x.lambda_m, x.T_s);
% cell current
delta_z = p.L_z / p.N;
I_cell = p.L_x * delta_z * sum(x.i_m,1);

% --------------------------------------------------------------------------------------   
% plots
ps = plot_presets();

% --------------------------------------------------------------------------------------   
% plots

% anode and cathode pressure
h1 = figure(1);
set(gcf, 'unit', 'normalized', 'position', ...
   [0.0740 0.2093 0.6469 0.5852]);

subplot(2,3,1);
plot(t, I_cell / (p.L_x*p.L_z)*1e-4,'Color', ps.color_rbmap(1,:), 'linestyle', line_style, 'linewidth', 1);
grid on; hold on;
xlabel('$t ~/~ \mathrm{s}$','fontsize', 12, 'Interpreter','latex');
ylabel('$ I ~/~ \mathrm{A ~cm^{-2}}$','fontsize',12, 'Interpreter','latex');
% set(gca,'FontSize',14);
axis([t(1), t(end), 0, 1.1*max(I_cell / (p.L_x*p.L_z)*1e-4)]);

subplot(2,3,2);
plot(t, x.U_cell,'Color', ps.color_rbmap(1,:), 'linestyle', line_style, 'linewidth', 1);
grid on; hold on;
xlabel('$t ~/~ \mathrm{s}$','fontsize', 12, 'Interpreter','latex');
ylabel('$U ~/~ \mathrm{V}$','fontsize',12, 'Interpreter','latex');
line([t(1), t(end)], [0.4,0.4], 'linestyle', '-.', 'color','k','linewidth', 1);
% line([t(1), t(end)], [1.1846,1.1846], 'linestyle', '-.', 'color','k','linewidth', 1);
axis([t(1), t(end), 0.3, 1]);


subplot(2,3,3);
plot(t, a_H2O_a(voli,:),'Color', ps.color_rblmap(1,:), 'linestyle', line_style, 'linewidth', 1);
grid on; hold on;
plot(t, a_H2O_c(voli,:),'Color', ps.color_rblmap(2,:), 'linestyle', line_style, 'linewidth', 1);
plot(t, a_H2O_ca(voli,:),'Color', ps.color_rbmap(1,:), 'linestyle', line_style, 'linewidth', 1);
plot(t, a_H2O_cc(voli,:),'Color', ps.color_rbmap(2,:), 'linestyle', line_style, 'linewidth', 1);
line([t(1), t(end)], [0.95,0.95], 'linestyle', '-.', 'color','k','linewidth', 1);
line([t(1), t(end)], [0.05,0.05], 'linestyle', '-.', 'color','k','linewidth', 1);
axis([t(1), t(end), -0.05 1.05]);
xlabel('$t ~/~ \mathrm{s}$','fontsize', 12, 'Interpreter','latex');
ylabel('$RH ~/~ \mathrm{1}$','fontsize',12,'Interpreter','latex');
legend('a_A','a_C','Location','southeast');


subplot(2,3,4);
plot(t, x.c_H2_a(voli,:),'Color', ps.color_rbmap(1,:), 'linestyle', line_style, 'linewidth', 1);
grid on; hold on;
plot(t, x.c_O2_c(voli,:),'Color', ps.color_rbmap(2,:), 'linestyle', line_style, 'linewidth', 1);
line([t(1), t(end)], [0,0], 'linestyle', '-.', 'color','k','linewidth', 1);
axis([t(1), t(end), -1, 120]);
xlabel('$t ~/~ \mathrm{s}$','fontsize', 12, 'Interpreter','latex');
ylabel('$c_{\mathrm{r} ~/~ \mathrm{mol~m^{-3}}}$','fontsize',12,'Interpreter','latex');
legend('c_{H2,a}','c_{O2,c}','Location','northeast');

subplot(2,3,5);
plot(t, x.p_a(voli,:),'Color', ps.color_rbmap(1,:), 'linestyle', line_style, 'linewidth', 1);
grid on; hold on;
plot(t, x.p_c(voli,:),'Color', ps.color_rbmap(2,:), 'linestyle', line_style, 'linewidth', 1);
xlabel('$t ~/~ s$','fontsize', 12, 'Interpreter','latex');
ylabel('$p ~/~ \mathrm{Pa}$','fontsize', 12, 'Interpreter','latex');
line([t(1), t(end)], [1*1e5,1*1e5], 'linestyle', '-.', 'color','k','linewidth', 1);
line([t(1), t(end)], [3.22*1e5,3.22*1e5], 'linestyle', '-.', 'color','k','linewidth', 1);
plot(t, u.p_a_out, 'Color', ps.color_rblmap(1,:),'linestyle', line_style2, 'linewidth', 1);
plot(t, u.p_c_out,'Color', ps.color_rblmap(2,:), 'linestyle', line_style2, 'linewidth', 1);
axis([t(1), t(end), 0.95*1e5, 3.3*1e5]);
legend('p_A','p_C','Location','southeast');


% temperatures of anode gas
subplot(2,3,6);
plot(t, x.T_a(voli,:),'Color', ps.color_rbmap(1,:), 'linestyle', line_style, 'linewidth', 1);
grid on; hold on;
plot(t, x.T_c(voli,:),'Color', ps.color_rbmap(2,:), 'linestyle', line_style, 'linewidth', 1);
plot(t, x.T_s(voli,:),'Color', ps.color_rbmap(4,:), 'linestyle', line_style, 'linewidth', 1);
xlabel('$t ~/~ \mathrm{s}$','fontsize', 12, 'Interpreter','latex');
ylabel('$T ~/~ \mathrm{K}$','fontsize',12,'Interpreter','latex');
line([t(1), t(end)], [273.15+30,273.15+30], 'linestyle', '-.', 'color','k','linewidth', 1);
line([t(1), t(end)], [273.15+90,273.15+90], 'linestyle', '-.', 'color','k','linewidth', 1);
plot(t, u.T_a_in,'Color', ps.color_rblmap(1,:), 'linestyle', line_style2, 'linewidth', 1);
plot(t, u.T_c_in,'Color', ps.color_rblmap(2,:), 'linestyle', line_style2, 'linewidth', 1);
plot(t, u.T_cool,'Color', ps.color_rblmap(4,:), 'linestyle', line_style2, 'linewidth', 1);
axis([t(1), t(end), 300, 370]);
legend('T_A','T_C','T_S','Location','southeast');

% exportgraphics(h1,['results_doe_01.emf'],'ContentType','vector');
% print -djpeg -r500 res\fig_sim_all_10op.jpg

