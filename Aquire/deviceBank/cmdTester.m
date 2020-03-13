function [] = cmdTester(queue,cmd,n)
    if cmd == 1
        tic;
        pause(n);
        send(queue,[num2str(toc) ':done\n']);
    else
        send(queue,'not started\n')
    end
end