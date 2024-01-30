%% Plot Model & Data
% Plot models and data associated with successful stop trials
% (Manuscript Figure 1-1)
% 1) 0% change
clear;
figure('Color','w'); % specify size!
width=600;
height=700;
set(gcf,'position',[100, -30 , width, height])

% Written by Darcy Diesburg (2023), Brown University
% figure settings
dip_colors = {[136 204 238], [51 34 136]; [76 145 65], [225 193 110]};
outspike_colors = {[152 152 51], [167 68 150], [204 102 119], [136 34 85]};
dip_ylims = [-4 9];
spike_ylims = [-10 280];
hist_ylims = [0 50];
xlims = [0 500];
sim_trials = 50;

% data
trial_spikes = load(strcat(pwd,'/HNN_sims/SSmodel/',...
    'spk_0.txt'));
% from param.txt, identities of drives and spikes - these can be found in
% the first 8 lines of the param .txt file
L2_pyr = [  35,  134];
L2_basket = [   0,   34];
L5_pyr = [ 170,  269];
L5_basket = [ 135,  169];
prox1 = [ 540,  809];
dist = [ 270,  539];
prox2 = [ 810, 1079];

% 1) Plot the raster of spikes in model pyr and basket cells
% load spk data from exemplar trial
subplot(2,1,1);
outspikes = trial_spikes(trial_spikes(:,2) <= L5_pyr(2),:);
% L2 basket
L2b = outspikes(outspikes(:,2) >= L2_basket(1) & outspikes(:,2) <= L2_basket(2),:);
scatter(L2b(:,1),L2b(:,2),30,outspike_colors{1}/255,'filled'); hold on;
% L2 pyr
L2p = outspikes(outspikes(:,2) >= L2_pyr(1) & outspikes(:,2) <= L2_pyr(2),:);
scatter(L2p(:,1),L2p(:,2),30,outspike_colors{2}/255,'filled');
% L5 basket
L5b = outspikes(outspikes(:,2) >= L5_basket(1) & outspikes(:,2) <= L5_basket(2),:);
scatter(L5b(:,1),L5b(:,2),30,outspike_colors{3}/255,'filled');
% L5 pyr
L5p = outspikes(outspikes(:,2) >= L5_pyr(1) & outspikes(:,2) <= L5_pyr(2),:);
scatter(L5p(:,1),L5p(:,2),30,outspike_colors{4}/255,'filled');
set(gca, 'YDir','reverse'); 
ylim(spike_ylims); xlim(xlims); set(gca,'XTick',[]); set(gca,'YTick',[])
hold off; box off;
%pos = get(gca, 'Position');
%set(gca, 'Position',pos+[0 0.05 0 0.03])

subplot(2,1,2);
% plot individual trial    
% load simulation
trial_simulation = load(strcat(pwd,'/HNN_sims/SSmodel/','dpl_0.txt'));
simulation_time = trial_simulation(:,1); %time
trial_simulation_agg = trial_simulation(:,2);%aggregate dipole

% % plot aggregate dipoles
hold on;
plot(simulation_time, trial_simulation_agg,'Color',...
    dip_colors{1,1}./255,'Linewidth',3); hold on;

%% 100% change
figure('Color','w'); % specify size!
width=600;
height=700;
set(gcf,'position',[100, -30 , width, height])

trial_spikes = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_100perc/',...
    'spk.txt'));

% 1) Plot the raster of spikes in model pyr and basket cells
% load spk data from exemplar trial
subplot(2,1,1);
outspikes = trial_spikes(trial_spikes(:,2) <= L5_pyr(2),:);
% L2 basket
L2b = outspikes(outspikes(:,2) >= L2_basket(1) & outspikes(:,2) <= L2_basket(2),:);
scatter(L2b(:,1),L2b(:,2),30,outspike_colors{1}/255,'filled'); hold on;
% L2 pyr
L2p = outspikes(outspikes(:,2) >= L2_pyr(1) & outspikes(:,2) <= L2_pyr(2),:);
scatter(L2p(:,1),L2p(:,2),30,outspike_colors{2}/255,'filled');
% L5 basket
L5b = outspikes(outspikes(:,2) >= L5_basket(1) & outspikes(:,2) <= L5_basket(2),:);
scatter(L5b(:,1),L5b(:,2),30,outspike_colors{3}/255,'filled');
% L5 pyr
L5p = outspikes(outspikes(:,2) >= L5_pyr(1) & outspikes(:,2) <= L5_pyr(2),:);
scatter(L5p(:,1),L5p(:,2),30,outspike_colors{4}/255,'filled');
set(gca, 'YDir','reverse'); 
ylim(spike_ylims); xlim(xlims); set(gca,'XTick',[]); set(gca,'YTick',[])
hold off; box off;
%pos = get(gca, 'Position');
%set(gca, 'Position',pos+[0 0.05 0 0.03])

subplot(2,1,2);
% plot individual trial    
% load simulation
trial_simulation = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_100perc/',...
    'dpl.txt'));
simulation_time = trial_simulation(:,1); %time
trial_simulation_agg_100 = trial_simulation(:,2);%aggregate dipole

% % plot aggregate dipoles
hold on;
plot(simulation_time, trial_simulation_agg_100,'Color',...
    dip_colors{1,1}./255,'Linewidth',3); hold on;
plot(simulation_time, trial_simulation_agg,'Color',...
    dip_colors{1,2}./255,'Linewidth',3);
%% 200% change
figure('Color','w'); % specify size!
width=600;
height=700;
set(gcf,'position',[100, -30 , width, height])

trial_spikes = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_200perc/',...
    'spk.txt'));

% 1) Plot the raster of spikes in model pyr and basket cells
% load spk data from exemplar trial
subplot(2,1,1);
outspikes = trial_spikes(trial_spikes(:,2) <= L5_pyr(2),:);
% L2 basket
L2b = outspikes(outspikes(:,2) >= L2_basket(1) & outspikes(:,2) <= L2_basket(2),:);
scatter(L2b(:,1),L2b(:,2),30,outspike_colors{1}/255,'filled'); hold on;
% L2 pyr
L2p = outspikes(outspikes(:,2) >= L2_pyr(1) & outspikes(:,2) <= L2_pyr(2),:);
scatter(L2p(:,1),L2p(:,2),30,outspike_colors{2}/255,'filled');
% L5 basket
L5b = outspikes(outspikes(:,2) >= L5_basket(1) & outspikes(:,2) <= L5_basket(2),:);
scatter(L5b(:,1),L5b(:,2),30,outspike_colors{3}/255,'filled');
% L5 pyr
L5p = outspikes(outspikes(:,2) >= L5_pyr(1) & outspikes(:,2) <= L5_pyr(2),:);
scatter(L5p(:,1),L5p(:,2),30,outspike_colors{4}/255,'filled');
set(gca, 'YDir','reverse'); 
ylim(spike_ylims); xlim(xlims); set(gca,'XTick',[]); set(gca,'YTick',[])
hold off; box off;
%pos = get(gca, 'Position');
%set(gca, 'Position',pos+[0 0.05 0 0.03])

subplot(2,1,2);
% plot individual trial    
% load simulation
trial_simulation = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_200perc/',...
    'dpl.txt'));
simulation_time = trial_simulation(:,1); %time
trial_simulation_agg_200 = trial_simulation(:,2);%aggregate dipole

% % plot aggregate dipoles
hold on;
plot(simulation_time, trial_simulation_agg_200,'Color',...
    dip_colors{1,1}./255,'Linewidth',3); hold on;
plot(simulation_time, trial_simulation_agg_100,'Color',...
    dip_colors{1,2}./255,'Linewidth',3);
%% 300% change
figure('Color','w'); % specify size!
width=600;
height=700;
set(gcf,'position',[100, -30 , width, height])

trial_spikes = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_300perc/',...
    'spk.txt'));

% 1) Plot the raster of spikes in model pyr and basket cells
% load spk data from exemplar trial
subplot(2,1,1);
outspikes = trial_spikes(trial_spikes(:,2) <= L5_pyr(2),:);
% L2 basket
L2b = outspikes(outspikes(:,2) >= L2_basket(1) & outspikes(:,2) <= L2_basket(2),:);
scatter(L2b(:,1),L2b(:,2),30,outspike_colors{1}/255,'filled'); hold on;
% L2 pyr
L2p = outspikes(outspikes(:,2) >= L2_pyr(1) & outspikes(:,2) <= L2_pyr(2),:);
scatter(L2p(:,1),L2p(:,2),30,outspike_colors{2}/255,'filled');
% L5 basket
L5b = outspikes(outspikes(:,2) >= L5_basket(1) & outspikes(:,2) <= L5_basket(2),:);
scatter(L5b(:,1),L5b(:,2),30,outspike_colors{3}/255,'filled');
% L5 pyr
L5p = outspikes(outspikes(:,2) >= L5_pyr(1) & outspikes(:,2) <= L5_pyr(2),:);
scatter(L5p(:,1),L5p(:,2),30,outspike_colors{4}/255,'filled');
set(gca, 'YDir','reverse'); 
ylim(spike_ylims); xlim(xlims); set(gca,'XTick',[]); set(gca,'YTick',[])
hold off; box off;
%pos = get(gca, 'Position');
%set(gca, 'Position',pos+[0 0.05 0 0.03])

subplot(2,1,2);
% plot individual trial    
% load simulation
trial_simulation = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_300perc/',...
    'dpl.txt'));
simulation_time = trial_simulation(:,1); %time
trial_simulation_agg_300 = trial_simulation(:,2);%aggregate dipole

% % plot aggregate dipoles
hold on;
plot(simulation_time, trial_simulation_agg_300,'Color',...
    dip_colors{1,1}./255,'Linewidth',3); hold on;
plot(simulation_time, trial_simulation_agg_200,'Color',...
    dip_colors{1,2}./255,'Linewidth',3);
%% 400% change
figure('Color','w'); % specify size!
width=600;
height=700;
set(gcf,'position',[100, -30 , width, height])

trial_spikes = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_400perc/',...
    'spk.txt'));

% 1) Plot the raster of spikes in model pyr and basket cells
% load spk data from exemplar trial
subplot(2,1,1);
outspikes = trial_spikes(trial_spikes(:,2) <= L5_pyr(2),:);
% L2 basket
L2b = outspikes(outspikes(:,2) >= L2_basket(1) & outspikes(:,2) <= L2_basket(2),:);
scatter(L2b(:,1),L2b(:,2),30,outspike_colors{1}/255,'filled'); hold on;
% L2 pyr
L2p = outspikes(outspikes(:,2) >= L2_pyr(1) & outspikes(:,2) <= L2_pyr(2),:);
scatter(L2p(:,1),L2p(:,2),30,outspike_colors{2}/255,'filled');
% L5 basket
L5b = outspikes(outspikes(:,2) >= L5_basket(1) & outspikes(:,2) <= L5_basket(2),:);
scatter(L5b(:,1),L5b(:,2),30,outspike_colors{3}/255,'filled');
% L5 pyr
L5p = outspikes(outspikes(:,2) >= L5_pyr(1) & outspikes(:,2) <= L5_pyr(2),:);
scatter(L5p(:,1),L5p(:,2),30,outspike_colors{4}/255,'filled');
set(gca, 'YDir','reverse'); 
ylim(spike_ylims); xlim(xlims); set(gca,'XTick',[]); set(gca,'YTick',[])
hold off; box off;
%pos = get(gca, 'Position');
%set(gca, 'Position',pos+[0 0.05 0 0.03])

subplot(2,1,2);
% plot individual trial    
% load simulation
trial_simulation = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_400perc/',...
    'dpl.txt'));
simulation_time = trial_simulation(:,1); %time
trial_simulation_agg_400 = trial_simulation(:,2);%aggregate dipole

% % plot aggregate dipoles
hold on;
plot(simulation_time, trial_simulation_agg_400,'Color',...
    dip_colors{1,1}./255,'Linewidth',3); hold on;
plot(simulation_time, trial_simulation_agg_300,'Color',...
    dip_colors{1,2}./255,'Linewidth',3);
%% 500% change
figure('Color','w'); % specify size!
width=600;
height=700;
set(gcf,'position',[100, -30 , width, height])

trial_spikes = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_500perc/',...
    'spk.txt'));

% 1) Plot the raster of spikes in model pyr and basket cells
% load spk data from exemplar trial
subplot(2,1,1);
outspikes = trial_spikes(trial_spikes(:,2) <= L5_pyr(2),:);
% L2 basket
L2b = outspikes(outspikes(:,2) >= L2_basket(1) & outspikes(:,2) <= L2_basket(2),:);
scatter(L2b(:,1),L2b(:,2),30,outspike_colors{1}/255,'filled'); hold on;
% L2 pyr
L2p = outspikes(outspikes(:,2) >= L2_pyr(1) & outspikes(:,2) <= L2_pyr(2),:);
scatter(L2p(:,1),L2p(:,2),30,outspike_colors{2}/255,'filled');
% L5 basket
L5b = outspikes(outspikes(:,2) >= L5_basket(1) & outspikes(:,2) <= L5_basket(2),:);
scatter(L5b(:,1),L5b(:,2),30,outspike_colors{3}/255,'filled');
% L5 pyr
L5p = outspikes(outspikes(:,2) >= L5_pyr(1) & outspikes(:,2) <= L5_pyr(2),:);
scatter(L5p(:,1),L5p(:,2),30,outspike_colors{4}/255,'filled');
set(gca, 'YDir','reverse'); 
ylim(spike_ylims); xlim(xlims); set(gca,'XTick',[]); set(gca,'YTick',[])
hold off; box off;
%pos = get(gca, 'Position');
%set(gca, 'Position',pos+[0 0.05 0 0.03])

subplot(2,1,2);
% plot individual trial    
% load simulation
trial_simulation = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_500perc/',...
    'dpl.txt'));
simulation_time = trial_simulation(:,1); %time
trial_simulation_agg_500 = trial_simulation(:,2);%aggregate dipole

% % plot aggregate dipoles
hold on;
plot(simulation_time, trial_simulation_agg_500,'Color',...
    dip_colors{1,1}./255,'Linewidth',3); hold on;
plot(simulation_time, trial_simulation_agg_400,'Color',...
    dip_colors{1,2}./255,'Linewidth',3);