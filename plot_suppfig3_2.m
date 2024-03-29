%% Plot Model & Data
% Plot models and data associated with successful stop trials
% (Manuscript Supplementary Figure 3-4)
clear;
figure('Color','w'); % specify size!
width=600;
height=450;
set(gcf,'position',[10, 10 , width, height])

% Written by Darcy Diesburg (2023), Brown University
% figure settings
dip_colors = {[136 204 238], [51 34 136]; [76 145 65], [225 193 110]};
outspike_colors = {[152 152 51], [167 68 150], [204 102 119], [136 34 85]};
inspike_colors = {[170 0 0], [0 128 0]};
dip_ylims = [-7 11];
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

% 1) Plot histogram of incoming prox, dist drive spikes
inspikes = trial_spikes(trial_spikes(:,2) > L5_pyr(2),:);
proxin1 = inspikes(inspikes(:,2) >= prox1(1) & inspikes(:,2) <= prox1(2),1);
proxin2 = inspikes(inspikes(:,2) >= prox2(1) & inspikes(:,2) <= prox2(2),1);
distin = inspikes(inspikes(:,2) >= dist(1) & inspikes(:,2) <= dist(2),1);
edges = linspace(0,500);
subplot(2,1,1);
histogram(distin,edges,'EdgeColor','none','FaceColor',inspike_colors{2}/255,'FaceAlpha',1); hold on;
histogram([proxin1; proxin2],edges,'EdgeColor','none','FaceColor',inspike_colors{1}/255,'FaceAlpha',1); 
ylim(hist_ylims); xlim(xlims); set(gca,'XTick',[])
hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0.2 0 -0.2])

subplot(2,1,2);
% plot individual trials
for trial = 1:sim_trials
    
    % load simulation
    trial_simulation = load(strcat(pwd,'/HNN_sims/SSmodel/',...
        'dpl_',num2str(trial-1),'.txt'));
    simulation_time = trial_simulation(:,1); %time
    trial_simulation_agg = trial_simulation(:,2);%aggregate dipole
    trial_simulation_L2 = trial_simulation(:,3); %layerII/III dipole
    trial_simulation_L5 = trial_simulation(:,4);%layerV dipole
    
    % plot aggregate dipoles
    hold on;
    plot(simulation_time, trial_simulation_agg,...
        'Color',[173, 173, 173]./255)
end

% plot average
% load ERP
data = load(strcat('SUCCSTOP.txt'));
data_time = data(:,1); % time
data = data(:,2); % voltage

% load simulation
simulation = load(strcat(pwd,'/HNN_sims/SSmodel/dpl.txt'));
simulation_time = simulation(:,1); % time
simulation_agg = simulation(:,2); % aggregate dipole
simulation_L2 = simulation(:,3); % layerII/III dipole
simulation_L5 = simulation(:,4); % layerV dipole

%plot data and aggregate dipoles
l1(2) = plot(simulation_time, simulation_agg,'Color',...
    dip_colors{1,1}./255,'Linewidth',3); hold on;
l1(1) = plot(data_time, data ,'Color',...
    dip_colors{1,2}./255,'Linewidth',3); hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0 0 0.2])
%legend(l1, 'Data','Model'); 
xlabel('Time (ms)');
ylabel('Amplitude (mV)');
ylim(dip_ylims);
xlim(xlims); xticks([100 200 300 400 500]);
hold off;

%% Plot models and data associated with SS drives run with default network params
clear;
figure('Color','w'); % specify size!
width=600;
height=450;
set(gcf,'position',[10, 10 , width, height])

% figure settings
dip_colors = {[136 204 238], [51 34 136]; [76 145 65], [225 193 110]};
outspike_colors = {[152 152 51], [167 68 150], [204 102 119], [136 34 85]};
inspike_colors = {[170 0 0], [0 128 0]};
dip_ylims = [-7 11];
spike_ylims = [-10 280];
hist_ylims = [0 50];
xlims = [0 500];
sim_trials = 50;

% data
trial_spikes = load(strcat(pwd,'/HNN_sims/supp figs/default_w_SSdrives/',...
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

% 1) Plot histogram of incoming prox, dist drive spikes
inspikes = trial_spikes(trial_spikes(:,2) > L5_pyr(2),:);
proxin1 = inspikes(inspikes(:,2) >= prox1(1) & inspikes(:,2) <= prox1(2),1);
proxin2 = inspikes(inspikes(:,2) >= prox2(1) & inspikes(:,2) <= prox2(2),1);
distin = inspikes(inspikes(:,2) >= dist(1) & inspikes(:,2) <= dist(2),1);
edges = linspace(0,500);
subplot(2,1,1);
histogram(distin,edges,'EdgeColor','none','FaceColor',inspike_colors{2}/255,'FaceAlpha',1); hold on;
histogram([proxin1; proxin2],edges,'EdgeColor','none','FaceColor',inspike_colors{1}/255,'FaceAlpha',1); 
ylim(hist_ylims); xlim(xlims); set(gca,'XTick',[])
hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0.2 0 -0.2])

subplot(2,1,2);
% plot individual trials
for trial = 1:sim_trials
    
    % load simulation
    trial_simulation = load(strcat(pwd,'/HNN_sims/supp figs/default_w_SSdrives/',...
        'dpl_',num2str(trial-1),'.txt'));
    simulation_time = trial_simulation(:,1); %time
    trial_simulation_agg = trial_simulation(:,2);%aggregate dipole
    trial_simulation_L2 = trial_simulation(:,3); %layerII/III dipole
    trial_simulation_L5 = trial_simulation(:,4);%layerV dipole
    
    % plot aggregate dipoles
    hold on;
    plot(simulation_time, trial_simulation_agg,...
        'Color',[173, 173, 173]./255)
end

% plot average
% load ERP
data = load(strcat('SUCCSTOP.txt'));
data_time = data(:,1); % time
data = data(:,2); % voltage

% load simulation
simulation = load(strcat(pwd,'/HNN_sims/supp figs/default_w_SSdrives/dpl.txt'));
simulation_time = simulation(:,1); % time
simulation_agg = simulation(:,2); % aggregate dipole
simulation_L2 = simulation(:,3); % layerII/III dipole
simulation_L5 = simulation(:,4); % layerV dipole

%plot data and aggregate dipoles
l1(2) = plot(simulation_time, simulation_agg,'Color',...
    dip_colors{1,1}./255,'Linewidth',3); hold on;
l1(1) = plot(data_time, data ,'Color',...
    dip_colors{1,2}./255,'Linewidth',3); hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0 0 0.2])
%legend(l1, 'Data','Model'); 
xlabel('Time (ms)');
ylabel('Amplitude (mV)');
ylim(dip_ylims);
xlim(xlims); xticks([100 200 300 400 500]);
hold off;