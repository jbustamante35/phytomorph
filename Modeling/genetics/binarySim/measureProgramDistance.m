function [d] = measureProgramDistance(hPROGRAM,PROGRAM)
    if issparse(hPROGRAM)
        hPROGRAM = uint8(full(hPROGRAM));
    end
    if issparse(PROGRAM)
        PROGRAM = uint8(full(PROGRAM));
    end
    d = bitwise_hamming(hPROGRAM,PROGRAM);
end