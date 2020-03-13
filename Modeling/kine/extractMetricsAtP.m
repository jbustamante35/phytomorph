function [f,j] = extractMetricsAtP(metricFunc,p,h)
    f = metricFunc(p);
    j = myJacobian(metricFunc,p,h);
    %tmpFunc = @(x)metricFunc(x)';
    %close all
    %[j,err] = jacobianest(metricFunc,p,.005,2.001);
end