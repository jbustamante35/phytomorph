function [] = writeToDisk(oPath,rise_data_out,recover_data_out,rise_affineSequence,recover_affineSequence,p)
    matFile = [oPath '{P_[' num2str(p) ']}.mat'];
    save(matFile,'rise_data_out','recover_data_out','rise_affineSequence','recover_data_out','p')
end
            