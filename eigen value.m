%% =========================================================
%  SIMULATION 3: EIGENVALUE SENSITIVITY ANALYSIS
%  Perturb the stiffness matrix and track how natural
%  frequencies (eigenvalues) shift — structural robustness.
%% =========================================================
clc; clear; close all;

fprintf('============================================================\n');
fprintf('  SIMULATION 3: Eigenvalue Sensitivity Analysis\n');
fprintf('============================================================\n\n');

%% --- STRUCTURAL SYSTEM ---
% 3-DOF spring-mass system stiffness matrix
K = [6 -2  0;
    -2  5 -3;
     0 -3  4];

M = eye(3);   % unit mass matrix (normalized)

% Nominal eigenvalues (natural frequencies squared)
[V_nom, D_nom] = eig(K, M);
lambda_nom = sort(diag(D_nom));
omega_nom  = sqrt(lambda_nom);   % natural frequencies (rad/s)

fprintf('Nominal eigenvalues (λ):        %.4f   %.4f   %.4f\n', lambda_nom(1), lambda_nom(2), lambda_nom(3));
fprintf('Nominal nat. frequencies (ω):   %.4f   %.4f   %.4f\n\n', omega_nom(1), omega_nom(2), omega_nom(3));

%% --- MONTE CARLO: PERTURB K ---
N          = 5000;
noise_pct  = 0.05;   % ±5% noise on stiffness entries
rng(42);

lambda_mc = zeros(3, N);
omega_mc  = zeros(3, N);

for i = 1:N
    dK = noise_pct * (rand(3,3) - 0.5) .* abs(K);
    dK = (dK + dK') / 2;            % keep symmetric
    K_p = K + dK;
    K_p = (K_p + K_p') / 2;        % enforce symmetry
    [~, D_p]      = eig(K_p, M);
    lam_p         = sort(diag(D_p));
    lambda_mc(:,i) = lam_p;
    omega_mc(:,i)  = sqrt(abs(lam_p));
end

%% --- STATISTICS PER MODE ---
fprintf('%-8s  %-10s  %-10s  %-10s  %-10s\n', 'Mode','λ nominal','λ mean','λ std','CV (%)');
for k = 1:3
    lm = mean(lambda_mc(k,:));
    ls = std(lambda_mc(k,:));
    fprintf('  λ%d       %.4f      %.4f      %.4f      %.2f\n', k, lambda_nom(k), lm, ls, 100*ls/lm);
end
fprintf('\n');

%% --- SENSITIVITY INDEX (std/nominal) ---
sens_idx = std(lambda_mc, 0, 2) ./ lambda_nom;
fprintf('Sensitivity index (σ/λ_nom):    λ1=%.4f  λ2=%.4f  λ3=%.4f\n\n', sens_idx(1), sens_idx(2), sens_idx(3));

%% --- PLOTTING ---
figure('Name','Eigenvalue Sensitivity Analysis','Color','k','Position',[50 50 1200 900]);
colors3 = {[0.3 0.9 0.5],[1 0.85 0.2],[0.5 0.4 1]};

% 1-3: Histograms for each mode
for k = 1:3
    subplot(3,3,k);
    histogram(lambda_mc(k,:), 60, 'FaceColor', colors3{k}, 'EdgeColor','none','FaceAlpha',0.85);
    hold on;
    xline(lambda_nom(k),'w--','LineWidth',2,'Label',sprintf('λ_{%d}^{nom}',k),'LabelColor','w','FontSize',9);
    xline(mean(lambda_mc(k,:)),'r:','LineWidth',1.5);
    set(gca,'Color','k','XColor','w','YColor','w');
    xlabel(sprintf('\\lambda_%d',k),'Color','w');
    ylabel('Count','Color','w');
    title(sprintf('Mode %d eigenvalue distribution',k),'Color','w','FontWeight','bold');
    grid on; set(gca,'GridColor',[0.3 0.3 0.3]);
end

% 4: Box plot of all three modes
subplot(3,3,4);
boxplot(lambda_mc', 'Labels',{'\lambda_1','\lambda_2','\lambda_3'}, ...
    'Colors','y','MedianStyle','target','Widths',0.4);
hold on;
plot(1:3, lambda_nom, 'wp', 'MarkerSize',10,'MarkerFaceColor','w');
set(gca,'Color','k','XColor','w','YColor','w');
title('Box plots — all modes','Color','w','FontWeight','bold');
ylabel('\lambda','Color','w');
grid on; set(gca,'GridColor',[0.3 0.3 0.3]);

% 5: Sensitivity index bar
subplot(3,3,5);
b5 = bar(sens_idx*100, 0.5, 'FaceColor','flat');
b5.CData = cell2mat(colors3');
set(gca,'Color','k','XColor','w','YColor','w','XTickLabel',{'\lambda_1','\lambda_2','\lambda_3'});
ylabel('Sensitivity index  σ/λ_{nom} (%)','Color','w');
title('Mode sensitivity ranking','Color','w','FontWeight','bold');
grid on; set(gca,'GridColor',[0.3 0.3 0.3]);

% 6: Natural frequency distribution (mode 1)
subplot(3,3,6);
histogram(omega_mc(1,:), 60,'FaceColor',colors3{1},'EdgeColor','none','FaceAlpha',0.85);
hold on;
xline(omega_nom(1),'w--','LineWidth',2);
set(gca,'Color','k','XColor','w','YColor','w');
xlabel('\omega_1  (rad/s)','Color','w'); ylabel('Count','Color','w');
title('Natural frequency \omega_1','Color','w','FontWeight','bold');
grid on; set(gca,'GridColor',[0.3 0.3 0.3]);

% 7: λ1 vs λ2 scatter
subplot(3,3,7);
scatter(lambda_mc(1,:), lambda_mc(2,:), 3, [0.4 0.7 1],'filled','MarkerFaceAlpha',0.3);
hold on;
plot(lambda_nom(1), lambda_nom(2),'rp','MarkerSize',12,'MarkerFaceColor','r');
set(gca,'Color','k','XColor','w','YColor','w');
xlabel('\lambda_1','Color','w'); ylabel('\lambda_2','Color','w');
title('Mode 1 vs Mode 2 scatter','Color','w','FontWeight','bold');
legend({'Perturbed','Nominal'},'TextColor','w','Color','k','EdgeColor','w','FontSize',8);
grid on; set(gca,'GridColor',[0.3 0.3 0.3]);

% 8: λ2 vs λ3 scatter
subplot(3,3,8);
scatter(lambda_mc(2,:), lambda_mc(3,:), 3, [1 0.7 0.2],'filled','MarkerFaceAlpha',0.3);
hold on;
plot(lambda_nom(2), lambda_nom(3),'rp','MarkerSize',12,'MarkerFaceColor','r');
set(gca,'Color','k','XColor','w','YColor','w');
xlabel('\lambda_2','Color','w'); ylabel('\lambda_3','Color','w');
title('Mode 2 vs Mode 3 scatter','Color','w','FontWeight','bold');
legend({'Perturbed','Nominal'},'TextColor','w','Color','k','EdgeColor','w','FontSize',8);
grid on; set(gca,'GridColor',[0.3 0.3 0.3]);

% 9: Running std of λ1 — convergence check
subplot(3,3,9);
run_std = zeros(1,N);
for i = 2:N, run_std(i) = std(lambda_mc(1,1:i)); end
plot(1:N, run_std,'Color',colors3{1},'LineWidth',1.5);
yline(std(lambda_mc(1,:)),'w--','LineWidth',1.5,'Label','Final σ','LabelColor','w','FontSize',9);
set(gca,'Color','k','XColor','w','YColor','w');
xlabel('Trials','Color','w'); ylabel('\sigma(\lambda_1)','Color','w');
title('Convergence of mode 1 σ','Color','w','FontWeight','bold');
grid on; set(gca,'GridColor',[0.3 0.3 0.3]);

sgtitle('Simulation 3: Eigenvalue Sensitivity Analysis','Color','w','FontSize',13,'FontWeight','bold');

fprintf('Simulation 3 complete.\n');
