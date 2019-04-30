function [ w ] = Likelihood( estimated_obsticles,max_laser_samples,obsticle_vector,sigma )

idxs = fix(linspace(1,length(estimated_obsticles),max_laser_samples));
idxM = knnsearch(obsticle_vector,estimated_obsticles(idxs,:));
Mu = obsticle_vector(idxM,:);
for ii = 1:max_laser_samples
    C(:,:,ii) = sigma*eye(2);
end
GM = gmdistribution(Mu,C);
w = prod(pdf(GM,estimated_obsticles(idxs,:)));
if w == 0
    w = eps;
end
end

% plot(obsticle_vector(:,1),obsticle_vector(:,2),'k.');
% hold on;
% plot(estimated_obsticles(:,1),estimated_obsticles(:,2),'r.');
