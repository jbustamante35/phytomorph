function [traceC] = traceVec_ver2(toMetricSpace,traceC,dv,N,alpha,haltF)
    for e = 1:N
        %{
        hmag = linspace(.01,1,100);
        for hb = 1:numel(hmag)
            [f,jT] = extractMetricsAtP(toMetricSpace,traceC(end,:),hmag(hb)*ones(1,size(traceC,2)),typicalX);
            dc(hb,:) = (pinv(jT')*dv)';
            
        end
        %}
        [f,jT] = extractMetricsAtP(toMetricSpace,traceC(end,:),.0001*ones(1,size(traceC,2)));
         
        
        
        
        
        %jT = j';
        
        %{
        % reported as working
        dc = (pinv(jT)*dv)';
        dc = alpha*(dc / norm(dc)); % was this line
        %}
        
        
        %{
        jT(:,end) = [];
        dc = (pinv(jT)*dv)';
        dc = alpha*(dc / norm(dc));
        dc = [dc 0];
        %}
        
        %[U,S,V] = svd(T);
        %{
        ijT = pinv(jT);
        ijT(end,:) = [];
        dc = (ijT*dv)';
        dc = alpha*(dc / norm(dc));
        dc = [dc 0];
        %}
        
        %jT = j;
        % experiment
        % reported as working
        %dc = (inv(jT')*dv)';
        dc = (inv(jT')*dv)';
        %dc(end) = 0;
        dc = alpha*(dc / norm(dc)); % was this line
        
        
        if any(isnan(dc) | isinf(dc))
            here = 1;
        end
        
        % focus on the mapping - bttom is the source top is the target
        % [d{target}/d{source}] = Tij = i,target j,source
        %dc = (pinv(jT)*dv)';
        %dc = dc / dc(end);
        %dc(1:3) = alpha*(dc(1:3)/norm(dc(1:3)));
        %dc = dc / dc(end);
        %dc = alpha*(dc / norm(dc)); % was this line
        traceC(e+1,:)  = traceC(e,:) + dc;
    end
end