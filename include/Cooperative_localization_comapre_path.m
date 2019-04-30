%% clear
clc;
clear;
close all;
%% ---- Load map ---- %%
convert_map;
map_image = imread('map/map.pgm');
map_imageCropped = map_image(750:1250,950:1400);
map_imageBW = map_imageCropped < 100;
map = robotics.OccupancyGrid(map_imageBW,20);
map.GridLocationInWorld = 0.5*[-25 -25];
inflate(map,0.01);
%map.show();

clear map_image map_imageCropped map_imageBW x_size y_size vars threashold resulution rawData1 newData1 name i ii jj map_raw;
%show(map);

%% ---- Load simulation data ---- %%
load('simulation/simulation_data.mat');

%% ---- Networks and estimators parameters ---- %%
Nf = 3; % Number of robots

sigma_L = 0.5; % Likelihood param
sigma_U = 0.5; % UCLT fusion param
Np = 200; % Number of particles
sigma_x_step = 0.5; % process noise
sigma_y_step = 0.5; % process noise
sigma_thets_step = 0.2; % process noise
Cnoise = [0.01 0 0;0 0.01 0; 0 0 0.001]; % process noise
fusion_reg = 0.01; % fusion reglarization param
max_laser_samples = 20;
total_time_steps = 100;
start_t = 200;
MonteCarloRuns = 1; % number of Monte Carlo simulations




currentScan = cell(Nf,1);

pose_history_PI = zeros(total_time_steps,6);
pose_history_No_fusion = zeros(total_time_steps,6);
pose_history_UCLT = zeros(total_time_steps,6);

%% initial filters

Xt = cell(Nf);
PF_PI.X = cell(Nf,1); % Particles intersection fusion network
PF_LPI.X = cell(Nf,1); % Linear Particles intersection fusion network
PF_S.X = cell(Nf,1); % No fusion
PF_UC.X = cell(Nf,1);% UCLT - centralized PF
PF_PI.W = cell(Nf,1);
PF_LPI.W = cell(Nf,1);
PF_S.W = cell(Nf,1);
PF_UC.W = cell(Nf,1);
for ii = 1:Nf
    PF_PI.X{ii} = zeros(Np,3);
    PF_LPI.X{ii} = zeros(Np,3);
    PF_S.X{ii} = zeros(Np,3);
    PF_UC.X{ii} = zeros(Np,3);
    PF_PI.W{ii} = ones(Np,1);
    PF_LPI.W{ii} = ones(Np,1);
    PF_S.W{ii} = ones(Np,1);
    PF_UC.W{ii} = ones(Np,1);
end

X1_pose = zeros(total_time_steps,3);
X2_pose = zeros(total_time_steps,3);
X3_pose = zeros(total_time_steps,3);
for ii = 1:total_time_steps
    X1_pose(ii,:) = X1{start_t+ii};
    X2_pose(ii,:) = X2{start_t+ii};
    X3_pose(ii,:) = X3{start_t+ii};
end



for Nmc = 1:MonteCarloRuns
    %% ---- Generate fusion network ---- %%
    A = GenerateFusionNetwork(Nf,0.2);
    
    %% ---- Initial starting position of the PFs ---- %%
    init_x(1,:) = X1{start_t};
    init_x(2,:) = X2{start_t};
    init_x(3,:) = X3{start_t};
    
    for ii = 1:Nf
        for jj = 1:Np
            PF_PI.X{ii}(jj,:) = init_x(ii,:) + mvnrnd(zeros(1,3),eye(3));
            PF_LPI.X{ii}(jj,:) = init_x(ii,:) + mvnrnd(zeros(1,3),eye(3));
            PF_S.X{ii}(jj,:) = init_x(ii,:) + mvnrnd(zeros(1,3),eye(3));
            PF_UC.X{ii}(jj,:) = init_x(ii,:) + mvnrnd(zeros(1,3),eye(3));
        end
    end
    counter = 1;
    t = start_t;
    end_t = start_t+total_time_steps;
    while t<end_t
        
        
        Xt{1,2} =  X1{t}-X2{t};
        Xt{2,1} =  X2{t}-X1{t};
        Xt{1,3} =  X1{t}-X3{t};
        Xt{3,1} =  X3{t}-X1{t};
        Xt{2,3} =  X2{t}-X3{t};
        Xt{3,2} =  X3{t}-X2{t};
        Xt{1,1} =  X1{t}-X1{t};
        Xt{2,2} =  X1{t}-X1{t};
        Xt{3,3} =  X1{t}-X1{t};
        
        scandata1 = scans1{t};
        
        V(1) = omometry_data1{t-1}(1);
        
        omega(1) = omometry_data1{t-1}(2);
        n1 = length(scandata1.Ranges);
        currentScan{1}  = lidarScan(scans1{t});
        currentScan{1} = removeInvalidData(currentScan{1},'RangeLimits',[0.1 12]);
        
        scandata2 = scans2{t};
        V(2) = omometry_data2{t-1}(1);
        omega(2) = omometry_data2{t-1}(2);
        n2 = length(scandata2.Ranges);
        currentScan{2}  = lidarScan(scans2{t});
        currentScan{2} = removeInvalidData(currentScan{2},'RangeLimits',[0.1 12]);
        
        scandata3 = scans3{t};
        V(3) = omometry_data3{t-1}(1);
        omega(3) = omometry_data3{t-1}(2);
        n3 = length(scandata3.Ranges);
        currentScan{3}  = lidarScan(scans3{t});
        currentScan{3} = removeInvalidData(currentScan{3},'RangeLimits',[0.1 12]);
        
        %% model new step
        for ii = 1:Nf
            for jj = 1:Np
                PF_PI.X{ii}(jj,:) = PF_PI.X{ii}(jj,:) + [(V(ii)*normrnd(1,sigma_x_step))*cos( PF_PI.X{ii}(jj,3)),(V(ii)*normrnd(1,sigma_y_step))*sin( PF_PI.X{ii}(jj,3)),omega(ii)*normrnd(1,sigma_thets_step)].*dt{t}+mvnrnd(zeros(1,3),Cnoise);
                PF_LPI.X{ii}(jj,:) = PF_LPI.X{ii}(jj,:) + [(V(ii)*normrnd(1,sigma_x_step))*cos( PF_LPI.X{ii}(jj,3)),(V(ii)*normrnd(1,sigma_y_step))*sin( PF_LPI.X{ii}(jj,3)),omega(ii)*normrnd(1,sigma_thets_step)].*dt{t}+mvnrnd(zeros(1,3),Cnoise);
                PF_S.X{ii}(jj,:) = PF_S.X{ii}(jj,:) + [(V(ii)*normrnd(1,sigma_x_step))*cos( PF_S.X{ii}(jj,3)),(V(ii)*normrnd(1,sigma_y_step))*sin( PF_S.X{ii}(jj,3)),omega(ii)*normrnd(1,sigma_thets_step)].*dt{t}+mvnrnd(zeros(1,3),Cnoise);
                PF_UC.X{ii}(jj,:) = PF_UC.X{ii}(jj,:) + [(V(ii)*normrnd(1,sigma_x_step))*cos( PF_UC.X{ii}(jj,3)),(V(ii)*normrnd(1,sigma_y_step))*sin( PF_UC.X{ii}(jj,3)),omega(ii)*normrnd(1,sigma_thets_step)].*dt{t}+mvnrnd(zeros(1,3),Cnoise);
            end
        end
        
        %% Likelihood
        for ii = 1:Nf
            for qq = 1:Np
                    Relitive_scan = transformScan(currentScan{ii},PF_PI.X{ii}(qq,:));
                    estimated_obsticles = Relitive_scan.Cartesian;
                    PF_PI.W{ii}(qq) = Likelihood(estimated_obsticles,max_laser_samples,obsticle_vector,sigma_L);
            end
            if sum(isnan(PF_PI.W{ii}))>=1
                PF_PI.W{ii} = ones(Np,1);
            end
            PF_PI.W{ii} = PF_PI.W{ii}./sum(PF_PI.W{ii});
            
            for qq = 1:Np
                    Relitive_scan = transformScan(currentScan{ii},PF_LPI.X{ii}(qq,:));
                    estimated_obsticles = Relitive_scan.Cartesian;
                    PF_LPI.W{ii}(qq) = Likelihood(estimated_obsticles,max_laser_samples,obsticle_vector,sigma_L);
            end
            if sum(isnan(PF_LPI.W{ii}))>=1
                PF_LPI.W{ii} = ones(Np,1);
            end
            PF_LPI.W{ii} = PF_LPI.W{ii}./sum(PF_LPI.W{ii});
            
            for qq = 1:Np
                    Relitive_scan = transformScan(currentScan{ii},PF_S.X{ii}(qq,:));
                    estimated_obsticles = Relitive_scan.Cartesian;
                    PF_S.W{ii}(qq) = Likelihood(estimated_obsticles,max_laser_samples,obsticle_vector,sigma_L);
            end
            if sum(isnan(PF_S.W{ii}))>=1
                PF_S.W{ii} = ones(Np,1);
            end
            PF_S.W{ii} = PF_S.W{ii}./sum(PF_S.W{ii});
            
            
            
            for qq = 1:Np
                    Relitive_scan = transformScan(currentScan{ii},PF_UC.X{ii}(qq,:));
                    estimated_obsticles = Relitive_scan.Cartesian;
                    PF_UC.W{ii}(qq) = Likelihood(estimated_obsticles,max_laser_samples,obsticle_vector,sigma_L);
            end
            if sum(isnan(PF_UC.W{ii}))>=1
                PF_UC.W{ii} = ones(Np,1);
            end
            PF_UC.W{ii} = PF_UC.W{ii}./sum(PF_UC.W{ii});
            
        end
        
        
        
        
        %%--------UPF--------%%
        for qq = 1:Np
            
            D12 = (PF_UC.X{1}(qq,1:2) - PF_UC.X{2}(qq,1:2))- Xt{1,2}(1:2);
            D13 = (PF_UC.X{1}(qq,1:2) - PF_UC.X{3}(qq,1:2))- Xt{1,3}(1:2);
            D23 = (PF_UC.X{2}(qq,1:2) - PF_UC.X{3}(qq,1:2))- Xt{2,3}(1:2);
            L12 = exp(-sigma_U*norm(D12));
            L13 = exp(-sigma_U*norm(D13));
            L23 = exp(-sigma_U*norm(D23));
            Uw(qq)=PF_UC.W{1}(qq)*PF_UC.W{2}(qq)*PF_UC.W{3}(qq)*L12*L13*L23;
        end
        if sum(isnan(Uw))>=1
            Uw = ones(Np,1);
        end
        Uw = Uw./sum(Uw);
        
        for ii = 1:Nf
            PF_UC.X{ii} = PF_UC.X{ii}(randsample(Np,Np,true,Uw),:);
        end
        %% fusion
        A(1,2) = det_localization( PF_LPI.X{1},PF_LPI.W{1},PF_LPI.X{2},PF_LPI.W{2},Xt{1,2}(1:2) );
        A(2,1) = A(1,2);
        A(1,3) = det_localization( PF_LPI.X{1},PF_LPI.W{1},PF_LPI.X{3},PF_LPI.W{3},Xt{1,3}(1:2) );
        A(3,1) = A(1,3);
        A(3,2) = det_localization( PF_LPI.X{3},PF_LPI.W{3},PF_LPI.X{2},PF_LPI.W{2},Xt{3,2}(1:2) );
        A(2,3) = A(3,2);
        A(1,1) = 1- A(1,2)-A(1,3);
        A(2,2) = 1- A(2,1)-A(2,3);
        A(3,3) = 1- A(3,2)-A(3,1);
        for ii = 1:Nf
            
            PF(:,:,ii) = PF_LPI.X{ii};
            W(:,ii) = PF_LPI.W{ii};
            
        end
        
        [ PF_out ] = LCPI( PF,W,A,Xt );
        
        for ii = 1:Nf
            PF_LPI.X{ii} = PF_out(:,:,ii);
        end
        %% resmpling
        for ii = 1:Nf
            PF_PI.X{ii} = PF_PI.X{ii}(randsample(Np,Np,true,PF_PI.W{ii}),:);
            PF_S.X{ii} = PF_S.X{ii}(randsample(Np,Np,true,PF_S.W{ii}),:);
            PF_PI.W{ii} = ones(Np,1)./Np;
            PF_S.W{ii} = ones(Np,1)./Np;
        end
        
      
        
        
        % find alpha using det
        fu = det_localization( PF_PI.X{1},PF_PI.W{1},PF_PI.X{2},PF_PI.W{2},Xt{1,2}(1:2) );
        PF_PI.W{1} = ParticlesIntersectionK_localization( PF_PI.X{1}, PF_PI.W{1},PF_PI.X{2}, PF_PI.W{2},fu,0.08 ,Xt{1,2}(1:2));
        PF_PI.X{1} = PF_PI.X{1}(randsample(Np,Np,true,PF_PI.W{1}),:)+ normrnd(0,fusion_reg,size(PF_PI.X{3}));
        PF_PI.W{1} = ones(Np,1)./Np;
        
        fu = det_localization( PF_PI.X{2},PF_PI.W{2},PF_PI.X{3},PF_PI.W{3},Xt{2,3}(1:2) );
        PF_PI.W{2} = ParticlesIntersectionK_localization( PF_PI.X{2}, PF_PI.W{2},PF_PI.X{3}, PF_PI.W{3},fu,0.08 ,Xt{2,3}(1:2));
        PF_PI.X{2} = PF_PI.X{2}(randsample(Np,Np,true,PF_PI.W{2}),:)+ normrnd(0,fusion_reg,size(PF_PI.X{3}));
        PF_PI.W{2} = ones(Np,1)./Np;
        
        fu = det_localization( PF_PI.X{3},PF_PI.W{3},PF_PI.X{1},PF_PI.W{1},Xt{3,1}(1:2) );
        PF_PI.W{3} = ParticlesIntersectionK_localization( PF_PI.X{3}, PF_PI.W{3},PF_PI.X{1}, PF_PI.W{1},fu,0.08 ,Xt{3,1}(1:2));
        PF_PI.X{3} = PF_PI.X{3}(randsample(Np,Np,true,PF_PI.W{3}),:)+ normrnd(0,fusion_reg,size(PF_PI.X{3}));
        PF_PI.W{3} = ones(Np,1)./Np;
        
        
        
        
        
        
        
        
        
        % ---- save history ---- %
        for ii = 1:Nf
            pose_history_PI(counter,2*ii-1) = mean(PF_PI.X{ii}(:,1));
            pose_history_PI(counter,2*ii) = mean(PF_PI.X{ii}(:,2));
            
            pose_history_No_fusion(counter,2*ii-1) = mean(PF_S.X{ii}(:,1));
            pose_history_No_fusion(counter,2*ii) = mean(PF_S.X{ii}(:,2));
            
            pose_history_UCLT(counter,2*ii-1) = mean(PF_UC.X{ii}(:,1));
            pose_history_UCLT(counter,2*ii) = mean(PF_UC.X{ii}(:,2));
        end
        
        %% PLOT visualization
        figure(1)
        show(map);
%         
%         hold on;
%         for ii = 1:Nf
%             scatter(PF_PI.X{ii}(:,1),PF_PI.X{ii}(:,2),'r.');
%             
%         end
%         
%         
%         
%         plot(X1{t}(1),X1{t}(2),'y+','MarkerSize',20);
%         plot(X2{t}(1),X2{t}(2),'y+','MarkerSize',20);
%         plot(X3{t}(1),X3{t}(2),'y+','MarkerSize',20);
%         %         scatter(obsticle_vector(:,1),obsticle_vector(:,2),'k.');
%         %         scatter(estimated_obsticles(:,1),estimated_obsticles(:,2),'r.');
%         hold off;
        
        hold on;
        for ii = 1:Nf
            plot(pose_history_PI(1:counter,2*ii-1),pose_history_PI(1:counter,2*ii),'b-');
            plot(pose_history_No_fusion(1:counter,2*ii-1),pose_history_No_fusion(1:counter,2*ii),'g-');
            plot(pose_history_UCLT(1:counter,2*ii-1),pose_history_UCLT(1:counter,2*ii),'w-');
        end
        
        
        
        plot(X1_pose(1:counter,1),X1_pose(1:counter,2),'k--');
        plot(X2_pose(1:counter,1),X2_pose(1:counter,2),'k--');
        plot(X3_pose(1:counter,1),X3_pose(1:counter,2),'k--');
      
        hold off;
        
        
        xlabel('X[m]');
        ylabel('Y[m]');
        xlim([-11.5 10.5]);
        ylim([-10.5 10.5]);
        
        pause(0.0001);
        
        
        
        t = t+1;
        counter = counter + 1;
    end
end








