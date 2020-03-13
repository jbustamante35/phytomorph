function [file] = getTmpWriteTicket()
    

    remoteExchange='/iplant/home/nmiller/publicData/pTES/';
    [~,tm] = system('date +%s');
    tm(end) = [];
    
    [~,q] = system(['/mnt/snapper/nate/phytoKeyPool/computeWindow ' tm ' 300']);
    
    tm = str2double(tm);
    windowTime = round(tm / 300);
    windowTime = windowTime* 300;
    [~,hash] = system(['echo -n ' num2str(windowTime) '|sha256sum ']);
    hash=hash(1:end-4);
    
   
    file = [remoteExchange hash '_' num2str(randi(3,1))];
    
end