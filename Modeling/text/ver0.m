location = '/mnt/snapper/nate/abstractText/abstractText/';
cdir = dir(location)
cdir(1:2) = [];
for e = 1:numel(cdir)
    inFile = [location cdir(e).name];
    sentence = extractFileText(inFile);
    abstracts{e} = sentence;
    e
end
%%
T_abstracts = tokenizedDocument(abstracts);
T_abstracts = docfun(@lower,T_abstracts);

%%
NT_abstracts = removeStopWords(T_abstracts);
%%
NT_abstracts = normalizeWords(NT_abstracts,'Style','lemma');
NT_abstracts = erasePunctuation(NT_abstracts);
NT_abstracts = removeEmptyDocuments(NT_abstracts);
%%
tbl = context(NT_abstracts,'grps','contextLength',50);
head(tbl)
%%
T_abstracts = docfun(@lower,T_abstracts);
PT_abstracts = addPartOfSpeechDetails(T_abstracts);
SPT_abstracts = normalizeWords(PT_abstracts,'Style','lemma');
SSPT_abstracts = removeStopWords(SPT_abstracts);
MSSPT_abstracts = erasePunctuation(SSPT_abstracts);
MSSPT_abstracts = removeEmptyDocuments(MSSPT_abstracts);
%%
BOW = bagOfWords(NT_abstracts);
cleanedBag = removeInfrequentWords(BOW,2);
%%
numTopics = 4;
InitialTopicConcentration = .6;
wordConcentration = 30;
mdl = fitlda(cleanedBag,numTopics,'Solver','cvb0','FitTopicConcentration',0,'Verbose',1,'WordConcentration',wordConcentration,'InitialTopicConcentration',InitialTopicConcentration);
%%
figure;
for topicIdx = 1:numTopics
    subplot(2,4,topicIdx)
    wordcloud(mdl,topicIdx);
    title("Topic " + topicIdx)
end
%%
topicMixtures = transform(mdl,cleanedBag);
%%
for doc = 1:numel(MSSPT_abstracts)    
    [tIDX,topicScore] = mdl.predict(MSSPT_abstracts(1,doc))
end
%%
[S C U E L ERR LAM] = PCA_FIT_FULL(topicMixtures,3);
%%
%emb = trainWordEmbedding(MSSPT_abstracts,'Dimension',20,'NumEpochs',10,'Verbose',1,'NumNegativeSamples',20);
emb = trainWordEmbedding(NT_abstracts,'LossFunction','hs','Window',5,'Model','skipgram','Dimension',30,'NumEpochs',10,'Verbose',1,'NumNegativeSamples',20);
%%
words = emb.Vocabulary;
V1 = word2vec(emb,'phenotype');
V2 = word2vec(emb,'hormone');
V3 = word2vec(emb,'gravitropism');
V4 = word2vec(emb,'auxin');
V5 = word2vec(emb,'hydrotropism');
d1 = V4 - V2
d2 = V4 - V3
delta1 = V2 - V1;
delta2 = V3 - V1;
W1 = V2 + delta2;
W2 = V3 + delta1;
WT = .5*(W1 + W2);
t1 = vec2word(emb,delta1-d2,5,'Distance','euclidean')
t2 = vec2word(emb,delta2-d1,5,'Distance','euclidean')
word1 = vec2word(emb,W1,5,'Distance','euclidean');
word2 = vec2word(emb,W2,5,'Distance','euclidean');
word3 = vec2word(emb,WT,5,'Distance','euclidean');
%%
V1 = word2vec(emb,'root');
V2 = word2vec(emb,'hormone');
V4 = vec2word(emb,V2+V1,5,'Distance','euclidean')
V4 = vec2word(emb,V2+V1,5,'Distance','cosine')
%%
V1 = word2vec(emb,'man');
V2 = word2vec(emb,'king');
V3 = word2vec(emb,'woman');
V4 = vec2word(emb,V2-V1+V3)
%%
king = word2vec(emb,"king");
man = word2vec(emb,"man");
woman = word2vec(emb,"woman");
queen = word2vec(emb,"queen");
CEN = mean(cat(1,man,woman,king,queen));
king = king - CEN;
man = man - CEN;
woman = woman - CEN;
word = vec2word(emb,king - man + woman)
%%
V1 = word2vec(emb,'treatment');
V2 = word2vec(emb,'signal');
V3 = word2vec(emb,'cold');
V4 = vec2word(emb,V2-V1+V3)
%%
V1 = word2vec(emb,'treatment');
V2 = word2vec(emb,'protein');
V3 = word2vec(emb,'cold');
V4 = vec2word(emb,V2-V1+V3)
%%
V1 = word2vec(emb,'root');
V2 = word2vec(emb,'control');
V3 = word2vec(emb,'gravitropism');
V4 = vec2word(emb,V2-V1+V3)
%%
V1 = word2vec(emb,'hormone');
V2 = word2vec(emb,'phenotype');
V3 = word2vec(emb,'gravitropism');
V4 = vec2word(emb,(V2-V1)+(V3-V1),5,'Distance','euclidean')
%%
italy = word2vec(emb,"Italy");
rome = word2vec(emb,"Rome");
paris = word2vec(emb,"Paris");
word = vec2word(emb,italy - rome + paris)
%%
V1 = word2vec(emb,'man');
V2 = word2vec(emb,'king');
V3 = word2vec(emb,'woman');


V4 = vec2word(emb,V1 + (V2-V1)+(V3-V1),5,'Distance','euclidean')
%%
V1 = word2vec(emb,'phenotype');
V2 = word2vec(emb,'hormone');
V3 = word2vec(emb,'gravitropism');
V4 = vec2word(emb,V1 + (V2-V1)+(V3-V1),5,'Distance','euclidean')
%%
V1 = word2vec(emb,'treatment');
V2 = word2vec(emb,'hormone');
delta = word2vec(emb,'plant');
V3 = word2vec(emb,'heat');
V4 = vec2word(emb,V1 + (delta+V2-V1)+(V3-V1),15,'Distance','euclidean')
%%
V1 = word2vec(emb,'gravitropism');
V2 = word2vec(emb,'hormone');
V3 = word2vec(emb,'root');
V4 = vec2word(emb,(V1-V2),5,'Distance','euclidean')
V4 = vec2word(emb,(V1-V2),5,'Distance','cosine')
%%
S1 = string(normalizeWords(tokenizedDocument('hormone')))
S2 = string(normalizeWords(tokenizedDocument('gravitropism')))
S3 = string(normalizeWords(tokenizedDocument('root')))

V1 = word2vec(emb,S1);
V2 = word2vec(emb,S2);
V3 = word2vec(emb,S3);
V4 = vec2word(emb,V2-V1+V3)



