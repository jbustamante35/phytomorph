FilePath = '/mnt/tetra/nate/maizeSeedlingsW/';
FileList_J = {};
FileExt = {'j2son'};
FileList_J = gdig(FilePath,FileList_J,FileExt,1);
%%
angleStack = [];
CURVE_SIZE = 40;
cnt = 1;
TIM = [];
for doc = 1:100
    close all
    
    h1 = figure;
    h2 = figure;
    fid = fopen(FileList_J{doc});
    A = fread(fid,'*char')';
    fclose(fid);
    seedlingDoc = jsondecode(A);
    for e = 1:numel(seedlingDoc(1))
        for p = 1:numel(seedlingDoc(1).seedlingDoc(e).phenotypes.paths)
            pth = seedlingDoc(1).seedlingDoc(e).phenotypes.paths(p).d(:,1:2);
            pth = flip(pth,1);
            
            [pth,~,aLEN] = arcLength(pth,'arcLen');
            
            
            pth = pth(1:4:end,:);
            
            dpth = diff(pth,1,1);
            ang = atan2(-dpth(:,1),dpth(:,2));
            ang = unwrap(ang);
            %{
            figure(h1)
            plot(ang)
            hold on
            pause(.3)
            figure(h2)
            %}
            dang = diff(ang,1,1);
            
            if size(dang,1) > 200
                TIM(:,cnt) = dang(1:200);
                cnt = cnt + 1;
            end
            
            dang = ang;
            if size(dang,1) > CURVE_SIZE
                aLEN = im2colF(aLEN(1:size(dang,1))',[CURVE_SIZE 1],[1 1]);
                data = im2colF(dang,[CURVE_SIZE 1],[1 1]);
                data = cat(3,data,aLEN);
                data = permute(data,[1 3 2]);
                angleStack = cat(3,angleStack,data);
                
                
                
                
            end
            plot(pth(:,2),-pth(:,1),'r')
            hold on
            
        end
        drawnow
    end
end
%%
net = narnet(1:10,10);
for i = 1:size(TIM,1)
    for j = 1:size(TIM,2)
        TIMY{i,j} = TIM(i,j);
    end
end
%TIM = mat2cell(TIM',size(TIM,2),size(TIM,1));
[Xs,Xi,Ai,Ts] = preparets(net,{},{},TIMY(:,1:2));
%%
net = train(net,Xs,Ts,Xi,Ai);
[Y,Xf,Af] = net(Xs,Xi,Ai);
[netc,Xic,Aic] = closeloop(net,Xf,Af);
y2 = netc(cell(0,20),Xic,Aic);
%%
G = angleStack(:,1,:);
TH = linspace(min(G(:)),max(G(:)),10);
YD = discretize(angleStack(:,1,:),TH);
Y = TH(YD);
angleStackTrain = angleStack;
angleStackTrain(:,1,:) = Y;
%%
Y = squeeze(angleStackTrain(end,1,:));
X = angleStack(1:end-1,:,:);
X = reshape(X,[size(X,1) size(X,2) 1 size(X,3)]);
layers = [imageInputLayer([size(X,1) size(X,2) size(X,3)],'Normalization','none');
          convolution2dLayer([5,2],5);
          reluLayer();
          maxPooling2dLayer([1 1],'Stride',2);
          convolution2dLayer([5,1],5);
          reluLayer();
          maxPooling2dLayer([3 1],'Stride',2);
          fullyConnectedLayer(numel(TH)-1);
          softmaxLayer();
          classificationLayer()];
       options = trainingOptions('sgdm','minibatch',128,'ValidationPatience',30,'Shuffle','every-epoch','Plots','training-progress','MaxEpochs',100,'InitialLearnRate',.01,'ExecutionEnvironment','parallel');

    curveNetwork = trainNetwork(X,categorical(double(Y)),layers,options);
      
    %%
    X = angleStack(1:end-1,:,:);
X = reshape(X,[size(X,1) size(X,2) 1 size(X,3)]);
    LAB = ind2vec(squeeze(YD(end,1,:))');
    curveNetwork = patternnet(10);
    curveNetwork = train(curveNetwork,squeeze(X(:,1,:,:)),full(LAB),'useParallel','yes');
    %%
curveNetwork = narnet(1:10,5);

[Xs,Xi,Ai,Ts] = preparets(net,{},{},T);
net = train(net,Xs,Ts,Xi,Ai);
      %% try drawing
      close all
      
      initTH = X(:,:,:,30);
      for d = 1:500
          %probs = curveNetwork.predict(initTH(end-(CURVE_SIZE-2):end,:,:));
          probs = curveNetwork(initTH(end-(CURVE_SIZE-2):end,1,:));
          T = gendist(probs',1,1);
          initTH = cat(1,initTH,[TH(T) initTH(end,2)+1]);
          dTH = diff(initTH(:,1));
          CUR = 4*[cos(dTH) sin(dTH)];
          CUR = cumsum(CUR,1);
          plot(CUR(:,2),CUR(:,1));
          axis equal
          drawnow
          
      end
      