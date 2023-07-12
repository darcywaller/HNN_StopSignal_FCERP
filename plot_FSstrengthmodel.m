%% Plot Model & Data
% Plot models and data associated with successful stop trials
% (Manuscript Figure 4A)

% Written by Darcy Diesburg (2023), Brown University

clear;
% Figure Settings
colours = {[255 140 0], [255 201 157]; [76 145 65], [225 193 110]};
ylims = [-4 8];
xlims = [0 500];
sim_trials = 50;

% plot individual trials
for trial = 1:sim_trials
    
    % load simulation
    trial_simulation = load(strcat(pwd,'/HNN_sims/FSstrengthmodel/',...
        'dpl_',num2str(trial-1),'.txt'));
    simulation_time = trial_simulation(:,1); %time
    trial_simulation_agg = trial_simulation(:,2);%aggregate dipole
    trial_simulation_L2 = trial_simulation(:,3); %layerII/III dipole
    trial_simulation_L5 = trial_simulation(:,4);%layerV dipole
    
    % plot aggregate dipoles
    subplot(1,2,1)
    hold on
    plot(simulation_time, trial_simulation_agg,...
        'Color', [127, 127, 127]./255)
    
    % plot layer-specific dipoles
    subplot(1,2,2)
    hold on
    plot(simulation_time, trial_simulation_L5,...
        'Color', [127, 127, 127]./255)
    plot(simulation_time, trial_simulation_L2,...
        'Color', [127, 127, 127]./255)
end

% plot average
% load ERP
data = load(strcat('FAILSTOP.txt'));
data_time = data(:,1); % time
data = data(:,2); % voltage

% load simulation
simulation = load(strcat(pwd,'/HNN_sims/FSstrengthmodel/dpl.txt'));
simulation_time = simulation(:,1); % time
simulation_agg = simulation(:,2); % aggregate dipole
simulation_L2 = simulation(:,3); % layerII/III dipole
simulation_L5 = simulation(:,4); % layerV dipole

%plot data and aggregate dipoles
subplot(1,2,1)
l1(1) = plot(data_time, data ,'Color',...
    colours{1,1}./255,'Linewidth',3);
l1(2) = plot(simulation_time, simulation_agg,'Color',...
    colours{1,2}./255,'Linewidth',3);
xlabel('Time (ms)')
ylabel('Amplitude (mV)')
ylim(ylims)
xlim(xlims)
title(strcat('Failed stops'))
legend(l1, 'Data','Model')

% plot layer-specific dipoles
subplot(1,2,2)
l2(1) = plot(simulation_time, simulation_L5, 'Color',...
    colours{2,1}./255,'Linewidth',3);
l2(2) = plot(simulation_time, simulation_L2, 'Color',...
    colours{2,2}./255,'Linewidth',3);
xlabel('Time (ms)')
ylabel('Amplitude (mV)')
ylim(ylims)
xlim(xlims)
title('Layer-specific Dipoles')
legend(l2, 'Layer V','Layer II/III')
