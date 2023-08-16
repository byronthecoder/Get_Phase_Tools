
function [AMsigs,tAmp,FSamp]=getAmpSigs(fileName,inTxtFile,spGroup,chunksXfile)

addpath([userpath,'\TBXS\leong_SAMPH']);%MFB_coeffs,sample_advance
addpath([userpath,'\utils']);
% fileName='FR_F_AB4_2016_10_24_ModuleDiadoco_bababa';
% inTxtFile='textgrid_info_FR_GE.csv';
% spGroup=1;
% chunksXfile=1;
outRoot='C:/Users/lancia/Desktop/reading_corpora/DDK-all-ctrl-wav/';

root=[outRoot,'DDK-all-ctrl-wav/'];


df = readtable([outRoot,inTxtFile]);

FSsig=16000;
FSamp=1000 ;%sampling frequency of the AM signals
startMeasureT=0;

i=find(strcmp(fileName,df.filename));


wav_file_path=[root, df.filename{i},'.wav'];
txt_file_path=[root, df.filename{i},'.TextGrid'];
if exist(wav_file_path,'file')>0
    if spGroup==2
        condition=df.condition{i};  
        [tierNames,txtData]=readTextGrid(txt_file_path);
        if isempty(txtData)
            return
        end
        badguys=(cellfun(@isempty,txtData{1}(:,3)) | contains(txtData{1}(:,3),'#'));
        syllStartTs=[txtData{1}{~badguys,1}];
        syllEndTs=[txtData{1}{~badguys,2}];
        if length(syllStartTs)<2
            return;
        end
        realDur=syllEndTs(end)-syllStartTs(1);
        chunkDur=realDur/3;
        ChunkStarts=[syllStartTs(1),syllStartTs(1)+chunkDur,syllStartTs(1)+2*chunkDur];
        ChunkEnds=[ChunkStarts(2),ChunkStarts(3),syllEndTs(end)];

        [ChunkStartsT,chunkStartsSyllIdx]=get_closest_Vals_Idxs(ChunkStarts,syllStartTs);
        [ChunkEndsT,ChunkEndsSyllIdx]=get_closest_Vals_Idxs(ChunkEnds,syllEndTs);

        ChunkEndsSyllIdx=ChunkEndsSyllIdx+1;

        chunkSyllNs=ChunkEndsSyllIdx-chunkStartsSyllIdx;
        chunkSyllRates=  chunkSyllNs./(ChunkEndsT-ChunkStartsT);  
        syllRate=mean(chunkSyllRates); %syllable rate
    else
        [tierNames,txtData]=readTextGrid(txt_file_path);
        if isempty(txtData)
            return
        end
        badguys=(cellfun(@isempty,txtData{1}(:,3)) | contains(txtData{1}(:,3),'#'));
        syllStartTs=[txtData{1}{~badguys,1}];
        syllEndTs=[txtData{1}{~badguys,2}];
        ChunkStarts=df.startT(i);
        ChunkEnds=df.endT(i);
        ChunkStartsT=ChunkStarts;
        ChunkEndsT=ChunkEnds;
        chunkDur=ChunkEndsT-ChunkStartsT;
        chunkSyllRates=df.rate(i);
        chunkSyllNs=chunkDur*chunkSyllRates;
        ChunkStartsT=ChunkStartsT+startMeasureT;
        syllRate=mean(chunkSyllRates); %syllable rate
        condition='';
    end

     if (min(chunkSyllNs)<2) || (min(chunkSyllNs)==inf)
        return
     end
    %derive AM filter edge frequencies from sentences' durations
    segRate=(syllRate)*2; %segment rate

    stressRate=(syllRate)/3; %accent rate

    [audioIn,audioSr]=audioread(wav_file_path);
    if audioSr~=FSsig
        x=resample(audioIn,FSsig,audioSr);
    else
        x=audioIn;
    end
    
    AMsig=abs(hilbert(zscore(x)));
    
    lpfilt =fir1(500,[1/(FSsig/2) (0.75*segRate)/(FSsig/2)]);

    AMsig=filter(lpfilt,1,AMsig);
    AMsig=sample_advance(AMsig, floor(500/2), 1e-7);

    % newCentFreq = findOptOScFreq2(AMsig,pars.FSamp,[],[pars.startT,pars.endT]);

    AMsig=resample(AMsig,FSamp,FSsig);


    tAmp=[1:length(AMsig)]./FSamp;

    CF0 = [stressRate/2;stressRate;syllRate;segRate;max(segRate,syllRate*3)];
    edges=CF0(1:end-1)+diff(CF0)/2;

    [CF, bpfs] = MFB_coeffs(edges,FSamp,0);

    
    AMsigs=zeros(size(AMsig,1),3);
     for n = 1:length(CF)
        l_fil = bpfs(1,n);
        fil_n = bpfs(2:1+l_fil,n);
        nshift = floor(l_fil/2);
        AMfil = filter(fil_n, 1, AMsig);
        AMfil = sample_advance(AMfil, nshift, 1e-7);
        AMsigs(:,n) =AMfil;
     end

end