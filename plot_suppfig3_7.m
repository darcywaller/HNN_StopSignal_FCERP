%% Plot Model & Data
% (Manuscript Supplementary Figure 3-5)
clear;
figure('Color','w'); % specify size!
width=600;
height=450;
set(gcf,'position',[10, 10 , width, height])

% Written by Darcy Diesburg (2023), Brown University
% figure settings
dip_colors = {[136 204 238], [51 34 136]; [76 145 65], [225 193 110]};
inspike_colors = {[170 0 0], [0 128 0]};
dip_ylims = [-7 11];
hist_ylims = [0 50];
xlims = [0 500];
sim_trials = 50;

% data
trial_spikes = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_p2n2_50/spk_0.txt'));
% from param.txt, identities of drives and spikes - these can be found in
% the first 8 lines of the param .txt file
L5_pyr = [ 170,  269];
prox1 = [ 540,  809];
dist = [ 270,  539];

% 1) Plot histogram of incoming prox, dist drive spikes
inspikes = trial_spikes(trial_spikes(:,2) > L5_pyr(2),:);
proxin1 = inspikes(inspikes(:,2) >= prox1(1) & inspikes(:,2) <= prox1(2),1);
distin = inspikes(inspikes(:,2) >= dist(1) & inspikes(:,2) <= dist(2),1);
edges = linspace(0,500);
subplot(2,3,1);
histogram(distin,edges,'EdgeColor','none','FaceColor',inspike_colors{2}/255,'FaceAlpha',1); hold on;
histogram(proxin1,edges,'EdgeColor','none','FaceColor',inspike_colors{1}/255,'FaceAlpha',1); 
ylim(hist_ylims); xlim(xlims); set(gca,'XTick',[])
hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0.2 0 -0.2])

% p2 n2 only
subplot(2,3,4);
% plot individual trials
for trial = 1:sim_trials
    
    % load simulation
    trial_simulation = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_p2n2_50/',...
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
data = load(strcat('p2n2only.txt'));
data_time = data(:,1); % time
data = data(:,2); % voltage

% load simulation
simulation = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_p2n2_50/dpl.txt'));
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

% p3 only
% data
trial_spikes = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_p3only_50/spk_0.txt'));

% 1) Plot histogram of incoming prox, dist drive spikes
inspikes = trial_spikes(trial_spikes(:,2) > L5_pyr(2),:);
proxin2 = inspikes(inspikes(:,2) >= dist(1) & inspikes(:,2) <= dist(2),1);

subplot(2,3,2);
histogram(proxin2,edges,'EdgeColor','none','FaceColor',inspike_colors{1}/255,'FaceAlpha',1); 
ylim(hist_ylims); xlim(xlims); set(gca,'XTick',[])
hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0.2 0 -0.2])

% p2 n2 only
subplot(2,3,5);
% plot individual trials
for trial = 1:sim_trials
    
    % load simulation
    trial_simulation = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_p3only_50/',...
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
data = load(strcat('p3only.txt'));
data_time = data(:,1); % time
data = data(:,2); % voltage

% load simulation
simulation = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_p3only_50/dpl.txt'));
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

% sum of the dipoles
subplot(2,3,3);
histogram(distin,edges,'EdgeColor','none','FaceColor',inspike_colors{2}/255,'FaceAlpha',1); hold on;
histogram([proxin1; proxin2],edges,'EdgeColor','none','FaceColor',inspike_colors{1}/255,'FaceAlpha',1); 
ylim(hist_ylims); xlim(xlims); set(gca,'XTick',[])
hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0.2 0 -0.2])

data = load(strcat('SUCCSTOP.txt'));
data_time = data(:,1); % time
data = data(:,2); % voltage

simulation = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_p3only_50/dpl.txt'));
simulation_time = simulation(:,1); % time
simulation_agg = simulation(:,2); % aggregate dipole

% SS opt one spike
simulation_p2n2 = load(strcat(pwd,'/HNN_sims/supp figs/SSmodel_p2n2_50/dpl.txt'));
simulation_agg_p2n2 = simulation_p2n2(:,2); % aggregate dipole
simulation_both = simulation_agg+simulation_agg_p2n2;

subplot(2,3,6);
l1(2) = plot(simulation_time, simulation_both,'Color',...
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

%% Plot models and data associated with failed stop trials
clear;
figure('Color','w'); % specify size!
width=600;
height=450;
set(gcf,'position',[10, 10 , width, height])

% figure settings
dip_colors = {[212 85 0], [232 170 128]; [76 145 65], [225 193 110]};
inspike_colors = {[170 0 0], [0 128 0]};
dip_ylims = [-7 11];
hist_ylims = [0 50];
xlims = [0 500];
sim_trials = 50;

% data
trial_spikes = load(strcat(pwd,'/HNN_sims/supp figs/FStimingmodel_p2n2_only/spk_0.txt'));
% from param.txt, identities of drives and spikes - these can be found in
% the first 8 lines of the param .txt file
L5_pyr = [ 170,  269];
prox1 = [ 540,  809];
dist = [ 270,  539];

% 1) Plot histogram of incoming prox, dist drive spikes
inspikes = trial_spikes(trial_spikes(:,2) > L5_pyr(2),:);
proxin1 = inspikes(inspikes(:,2) >= prox1(1) & inspikes(:,2) <= prox1(2),1);
distin = inspikes(inspikes(:,2) >= dist(1) & inspikes(:,2) <= dist(2),1);
edges = linspace(0,500);
subplot(2,3,1);
histogram(distin,edges,'EdgeColor','none','FaceColor',inspike_colors{2}/255,'FaceAlpha',1); hold on;
histogram(proxin1,edges,'EdgeColor','none','FaceColor',inspike_colors{1}/255,'FaceAlpha',1); 
ylim(hist_ylims); xlim(xlims); set(gca,'XTick',[])
hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0.2 0 -0.2])

% p2 n2 only
subplot(2,3,4);
% plot individual trials
for trial = 1:sim_trials
    
    % load simulation
    trial_simulation = load(strcat(pwd,'/HNN_sims/supp figs/FStimingmodel_p2n2_only/',...
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
data = load(strcat('p2n2only_fs.txt'));
data_time = data(:,1); % time
data = data(:,2); % voltage

% load simulation
simulation = load(strcat(pwd,'/HNN_sims/supp figs/FStimingmodel_p2n2_only/dpl.txt'));
simulation_time = simulation(:,1); % time
simulation_agg = simulation(:,2); % aggregate dipole
simulation_L2 = simulation(:,3); % layerII/III dipole
simulation_L5 = simulation(:,4); % layerV dipole

%plot data and aggregate dipoles
l1(2) = plot(simulation_time, simulation_agg,'Color',...
    dip_colors{1,2}./255,'Linewidth',3); hold on;
l1(1) = plot(data_time, data ,'Color',...
    dip_colors{1,1}./255,'Linewidth',3); hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0 0 0.2])
%legend(l1, 'Data','Model'); 
xlabel('Time (ms)');
ylabel('Amplitude (mV)');
ylim(dip_ylims);
xlim(xlims); xticks([100 200 300 400 500]);
hold off;

% p3 only
% data
trial_spikes = load(strcat(pwd,'/HNN_sims/supp figs/FStimingmodel_p3_only/spk_0.txt'));

% 1) Plot histogram of incoming prox, dist drive spikes
inspikes = trial_spikes(trial_spikes(:,2) > L5_pyr(2),:);
proxin2 = inspikes(inspikes(:,2) >= dist(1) & inspikes(:,2) <= dist(2),1);

subplot(2,3,2);
histogram(proxin2,edges,'EdgeColor','none','FaceColor',inspike_colors{1}/255,'FaceAlpha',1); 
ylim(hist_ylims); xlim(xlims); set(gca,'XTick',[])
hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0.2 0 -0.2])

% p2 n2 only
subplot(2,3,5);
% plot individual trials
for trial = 1:sim_trials
    
    % load simulation
    trial_simulation = load(strcat(pwd,'/HNN_sims/supp figs/FStimingmodel_p3_only/',...
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
data = load(strcat('p3only_fs.txt'));
data_time = data(:,1); % time
data = data(:,2); % voltage

% load simulation
simulation = load(strcat(pwd,'/HNN_sims/supp figs/FStimingmodel_p3_only/dpl.txt'));
simulation_time = simulation(:,1); % time
simulation_agg = simulation(:,2); % aggregate dipole
simulation_L2 = simulation(:,3); % layerII/III dipole
simulation_L5 = simulation(:,4); % layerV dipole

%plot data and aggregate dipoles
l1(2) = plot(simulation_time, simulation_agg,'Color',...
    dip_colors{1,2}./255,'Linewidth',3); hold on;
l1(1) = plot(data_time, data ,'Color',...
    dip_colors{1,1}./255,'Linewidth',3); hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0 0 0.2])
%legend(l1, 'Data','Model'); 
xlabel('Time (ms)');
ylabel('Amplitude (mV)');
ylim(dip_ylims);
xlim(xlims); xticks([100 200 300 400 500]);
hold off;

% sum of the dipoles
subplot(2,3,3);
histogram(distin,edges,'EdgeColor','none','FaceColor',inspike_colors{2}/255,'FaceAlpha',1); hold on;
histogram([proxin1; proxin2],edges,'EdgeColor','none','FaceColor',inspike_colors{1}/255,'FaceAlpha',1); 
ylim(hist_ylims); xlim(xlims); set(gca,'XTick',[])
hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0.2 0 -0.2])

data = load(strcat('FAILSTOP.txt'));
data_time = data(:,1); % time
data = data(:,2); % voltage

simulation = load(strcat(pwd,'/HNN_sims/supp figs/FStimingmodel_p3_only/dpl.txt'));
simulation_time = simulation(:,1); % time
simulation_agg = simulation(:,2); % aggregate dipole

% SS opt one spike
simulation_p2n2 = load(strcat(pwd,'/HNN_sims/supp figs/FStimingmodel_p2n2_only/dpl.txt'));
simulation_agg_p2n2 = simulation_p2n2(:,2); % aggregate dipole
simulation_both = simulation_agg+simulation_agg_p2n2;

subplot(2,3,6);
l1(2) = plot(simulation_time, simulation_both,'Color',...
    dip_colors{1,2}./255,'Linewidth',3); hold on;
l1(1) = plot(data_time, data ,'Color',...
    dip_colors{1,1}./255,'Linewidth',3); hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0 0 0.2])
%legend(l1, 'Data','Model'); 
xlabel('Time (ms)');
ylabel('Amplitude (mV)');
ylim(dip_ylims);
xlim(xlims); xticks([100 200 300 400 500]);
hold off;