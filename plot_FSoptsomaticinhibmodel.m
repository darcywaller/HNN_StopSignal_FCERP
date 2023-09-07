%% Plot Model & Data
% Plot models and data associated with failed stop trials
% (Manuscript Figure 5B)
clear;
figure('Color','w'); % specify size!
width=600;
height=800;
set(gcf,'position',[10, 10 , width, height])

% Written by Darcy Diesburg (2023), Brown University
% figure settings
dip_colors = {[212 85 0], [232 170 128]; [76 145 65], [225 193 110]};
outspike_colors = {[152 152 51], [167 68 150], [204 102 119], [136 34 85]};
inspike_colors = {[170 0 0], [0 128 0]};
dip_ylims = [-4 9];
spike_ylims = [-10 280];
hist_ylims = [0 50];
xlims = [0 500];
sim_trials = 50;

% data
trial_spikes = load(strcat(pwd,'/HNN_sims/FSoptsomaticinhibmodel/',...
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
subplot(3,1,1);
histogram(distin,edges,'EdgeColor','none','FaceColor',inspike_colors{2}/255,'FaceAlpha',1); hold on;
histogram([proxin1; proxin2],edges,'EdgeColor','none','FaceColor',inspike_colors{1}/255,'FaceAlpha',1); 
ylim(hist_ylims); xlim(xlims); set(gca,'XTick',[])
hold off; box off;
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0.12 0 -0.12])

% 2) Plot the raster of spikes in model pyr and basket cells
% load spk data from exemplar trial
subplot(3,1,2);
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
pos = get(gca, 'Position');
set(gca, 'Position',pos+[0 0.09 0 0.05])

subplot(3,1,3);
% plot individual trials
for trial = 1:sim_trials
    
    % load simulation
    trial_simulation = load(strcat(pwd,'/HNN_sims/FSoptsomaticinhibmodel/',...
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
data = load(strcat('FAILSTOP.txt'));
data_time = data(:,1); % time
data = data(:,2); % voltage

% load simulation
simulation = load(strcat(pwd,'/HNN_sims/FSoptsomaticinhibmodel/dpl.txt'));
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
set(gca, 'Position',pos+[0 0 0 0.1])
%legend(l1, 'Data','Model'); 
xlabel('Time (ms)');
ylabel('Amplitude (mV)');
ylim(dip_ylims);
xlim(xlims); xticks([100 200 300 400 500]);

% plot layer-specific dipoles in inset
axes('Position',[.58 .13 .18 .1]);
box on;

l2(1) = plot(simulation_time, simulation_L2, 'Color',...
    dip_colors{2,2}./255,'Linewidth',2); hold on;
l2(2) = plot(simulation_time, simulation_L5, 'Color',...
    dip_colors{2,1}./255,'Linewidth',2);
ylim(dip_ylims); xlim(xlims); xticks([100 200 300 400 500]);
grid on; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); hold off; 
