CompRndApod = zeros(500,4,1024);
for set = 1 : 500
    list = [1 : 1024];
    for i = 1 : 4
        tmp = randperm(length(list),256);
        CompRndApod(set,i,list(tmp)) = 1;
        list(tmp) = [];
    end
end

%% Test
tmp = CompRndApod(500,1,:) + CompRndApod(500,2,:) + ...
    CompRndApod(500,3,:) + CompRndApod(500,4,:);
figure; imagesc(reshape(tmp,[32,32]));

%% Save
save('saved500CompRndApod','CompRndApod');