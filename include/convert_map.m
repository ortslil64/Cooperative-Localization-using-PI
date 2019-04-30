%% import map 
rawData1 = importdata('map/map.pgm');
[~,name] = fileparts('map/map.pgm');
newData1.(genvarname(name)) = rawData1; 
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));

end
%% map parameters from yaml 
map_raw = map;
threashold = 205;
obsticle_vector = [];
resulution = 0.05;

%% convert map to vector 
[y_size,x_size] = size(map_raw.cdata);
for ii = 1:x_size
    for jj = 1:y_size
        if map_raw.cdata(jj,ii) < threashold
            obsticle_vector = [resulution*[ii,-jj]+[-58.8654998, -52.684998]+[0,resulution*y_size];obsticle_vector];
        end
    end
end
