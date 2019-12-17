function txt = myupdatefcn(~,event_obj, posArray, posArray2, stations,jumpmag)
% Customizes text of data tips
pos = get(event_obj,'Position');
% I = get(event_obj, 'DataIndex');
I2 = find(posArray(:, 1) == pos(1) & posArray(:, 2) == pos(2));
if isempty(I2)
    I2 = find(posArray2(:, 1) == pos(1) & posArray2(:, 2) == pos(2));
end
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Jump: ',num2str(jumpmag(I2))], ...
       ['Station: ', stations{I2}]};