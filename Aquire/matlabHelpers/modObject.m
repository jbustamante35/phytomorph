function [] = modObject(baseObject,op,parameters)
    baseObject = fileread(baseObject);
    parameter = fileread(parameters);
    
    parameter = projectPtr(parameter);
    
    cmd = ['/mnt/scratch1/phytomorph_dev/Aquire/phytoMeta/arrayOp remove ''' baseObject ''' usbPorts ''' parameter ''''];
    remove
end