T_deoxy = table();
T_oxy   = table();
T_eeg = table();
T_eegoxy   = table();
T_eegdeoxy = table();
T_eegdeoxyoxy   = table();
T_oxydeoxy = table();

%预分配减少运算量
C = cell(33,5,29);   %这个地方的 5对应的折数需要修改为实际的数量
slopeall = cell(1,29);

EegMyDataDir = 'F:\2016_Data\EEG_01-29';
NirsMyDataDir = 'F:\2016_Data\NIRS_01-29';
WorkingDir = 'F:\2016_Data';

%subject名称
% parameters
subdir_list = {'\subject 01','\subject 02','\subject 03','\subject 04','\subject 05','\subject 06',...
    '\subject 07','\subject 08','\subject 09','\subject 10','\subject 11','\subject 12','\subject 13'...
    ,'\subject 14','\subject 15','\subject 16','\subject 17','\subject 18','\subject 19','\subject 20'...
    ,'\subject 21','\subject 22','\subject 23','\subject 24','\subject 25','\subject 26','\subject 27'...
    ,'\subject 28','\subject 29'}; % subject

subNum = 29;
for id = 1:subNum
    
    startup_bbci_toolbox('DataDir',EegMyDataDir,'TmpDir','F:\2016_Data\TMP');
    BTB.History = 0; % to aviod error for merging cnt
    
    % load nirs data
    loadDir = fullfile(NirsMyDataDir, subdir_list{id});
    cd(loadDir);
    load cnt; load mrk;load mnt; % load continous eeg signal (cnt), marker (mrk)
    cd(WorkingDir)
    
    % MBLL and save post-MBLL data
    for idx = 1 : 6
        filename{idx} = fullfile(subdir_list{id}, ['session',num2str(idx)]);
        cntHb{idx} = proc_BeerLambert(cnt{idx}, 'Opdist', 3, 'DPF', [5.98 7.15], 'Citation', 1);
        file_saveNIRSMatlab(filename{idx}, cntHb{idx}, mrk{idx}, mnt);
    end
    
    clear cnt cntHb mrk;
    % load deoxy- and oxy-hemoglobin data
    for idx = 1 : 6
        
        [cnt_temp.deoxy{idx}, mrk_temp{idx}, mnt] = file_loadNIRSMatlab(filename{idx}, 'Signal','deoxy');
        [cnt_temp.oxy{idx}  , ~, ~]   = file_loadNIRSMatlab(filename{idx}, 'Signal','oxy');
    end
    % merge cnts in each session
    % for motor imagery: imag
    % for mental arithmetic: ment
    
    [cnt.imag.deoxy, mrk.imag.deoxy] = proc_appendCnt({cnt_temp.deoxy{1}, cnt_temp.deoxy{3}, cnt_temp.deoxy{5}}, {mrk_temp{1}, mrk_temp{3}, mrk_temp{5}}); % merged motor imagery cnts
    [cnt.imag.oxy, mrk.imag.oxy]     = proc_appendCnt({cnt_temp.oxy{1}, cnt_temp.oxy{3}, cnt_temp.oxy{5}}, {mrk_temp{1}, mrk_temp{3}, mrk_temp{5}}); % merged motor imagery cnts
    clear cnt_temp;
    % band-pass filtering
    band_freq = [0.01 0.1]; % Hz
    ord.nirs = 3;
    [z,p,k] = cheby2(ord.nirs, 30, band_freq/10*2, 'bandpass');
    [SOS,G] = zp2sos(z,p,k);
    
    cnt.imag.deoxy = proc_filtfilt(cnt.imag.deoxy, SOS, G);
    cnt.imag.oxy   = proc_filtfilt(cnt.imag.oxy, SOS, G);
    
    % channel selection
    FrontalChannel = {'AF7Fp1','AF3Fp1','AF3AFz','FpzFp1','FpzAFz','FpzFp2','AF4AFz','AF4Fp2','AF8Fp2'};
    MotorChannel = {'C5CP5','C5FC5','C5C3','FC3FC5','FC3C3','FC3FC1','CP3CP5','CP3C3','CP3CP1','C1C3','C1FC1','C1CP1','C2FC2','C2CP2','C2C4','FC4FC2','FC4C4','FC4FC6','CP4CP6','CP4CP2','CP4C4','C6CP6','C6C4','C6FC6'};
    OccipitalChannel = {'OzPOz','OzO1','OzO2'};
    ParientalChannel = {'Pz','P3','P4','P7','P8'};
    
    cnt_org.imag.deoxy = cnt.imag.deoxy; % backup
    cnt_org.imag.oxy   = cnt.imag.oxy; % backup
    
    cnt.imag.deoxy = proc_selectChannels(cnt.imag.deoxy, [MotorChannel,ParientalChannel]);
    cnt.imag.oxy   = proc_selectChannels(cnt.imag.oxy, [MotorChannel,ParientalChannel]);
    
    % segmentation (epoching)
    ival_epo  = [-10 25]*1000; % epoch range (unit: msec)
    ival_base = [-3 0]*1000; % baseline correction range (unit: msec)
    
    epo.imag.deoxy = proc_segmentation(cnt.imag.deoxy, mrk.imag.deoxy, ival_epo);
    epo.imag.oxy   = proc_segmentation(cnt.imag.oxy, mrk.imag.oxy, ival_epo);
    epo.imag.deoxy = proc_baseline(epo.imag.deoxy,ival_base);
    epo.imag.oxy   = proc_baseline(epo.imag.oxy,ival_base);
    
    % 获取排序后的索引
    [~, deoxy_sort_indices] = sort(epo.imag.deoxy.y(1, :));
    [~, oxy_sort_indices] = sort(epo.imag.oxy.y(1, :));
    
    % 使用排序后的索引对矩阵和标签矩阵进行重新排列
    epo.imag.deoxy.x = epo.imag.deoxy.x(:, :, deoxy_sort_indices);
    epo.imag.deoxy.y = epo.imag.deoxy.y(:, deoxy_sort_indices);
    epo.imag.oxy.x = epo.imag.oxy.x(:, :, oxy_sort_indices);
    epo.imag.oxy.y = epo.imag.oxy.y(:, oxy_sort_indices);
    
    % 读取EEG数据
    loadDir = fullfile(EegMyDataDir, subdir_list{id},'with occular artifact');
    cd(loadDir);
    load cnt; load mrk, load mnt; % load continous eeg signal (cnt), marker (mrk) and montage (mnt)
    cd(WorkingDir)
    
    % merge cnts in each session
    % for motor imagery: imag
    % for mental arithmetic: ment
    cnt_temp = cnt; mrk_temp = mrk; % save data temporarily
    clear cnt mrk mnt;
    [cnt.imag.eeg, mrk.imag.eeg] = proc_appendCnt({cnt_temp{1}, cnt_temp{3}, cnt_temp{5}}, {mrk_temp{1}, mrk_temp{3}, mrk_temp{5}}); % merged motor imagery cnts
    clear cnt_temp;
    % Select EEG channels only (excluding EOG channels) for classification
    cnt.imag.eeg = proc_selectChannels(cnt.imag.eeg,'not','*EOG'); % remove EOG channels (VEOG, HEOG)
    
    % common average reference
    cnt.imag.eeg = proc_commonAverageReference(cnt.imag.eeg);
    
    % segmentation (epoching)
    ival_epo  = [-10 25]*1000; % epoch range (unit: msec)
    ival_base = [-3 0]*1000; % baseline correction range (unit: msec)
    
    epo.imag.eeg = proc_segmentation(cnt.imag.eeg, mrk.imag.eeg, ival_epo);
    epo.imag.eeg = proc_baseline(epo.imag.eeg,ival_base);
    
    % frequency band selection for common spatial pattern (CSP)
    MotorChannel = {'CFC5','CFC6','CFC3','CFC4','Cz,','CCP5','CCP6','CCP3','CCP4'};
    ParientalChannel = {'Pz','P3','P4','P7','P8'};
    FrontalChannel = {'F7','FAF5','F3','AFp1','FAF1','AFp2','FAF2','FAF6','F4','F8'};
    OccipitalChannel = {'PPO1','OPO1','OPO2','PPO2'};
    
    % channel selection
    cnt_org.imag.eeg = cnt.imag.eeg; % backup
    cnt.imag.eeg = proc_selectChannels(cnt.imag.eeg, [MotorChannel,ParientalChannel]);
    
    % narrow frequency band selection for CSP
    band_csp.imag = select_bandnarrow(cnt.imag.eeg, mrk.imag.eeg, [0 10]*1000); % band selection using 0~10 sec epoch for motor imagery
    
    % Cheby2 bandpass filter with a passband of band_csp, with at most Rp dB of passband ripple and at least Rs dB attenuation in the stopbands that are 3 Hz wide on both sides of the passband
    Wp.imag = band_csp.imag/epo.imag.eeg.fs*2;
    Ws.imag = [band_csp.imag(1)-3, band_csp.imag(end)+3]/epo.imag.eeg.fs*2;
    Rp.imag = 3; % in dB
    Rs.imag = 30; % in dB
    [ord.imag, Ws.imag] = cheb2ord(Wp.imag, Ws.imag, Rp.imag, Rs.imag);
    [filt_b.imag,filt_a.imag] = cheby2(ord.imag, Rs.imag, Ws.imag);
    
    epo.imag.eeg = proc_filtfilt(epo.imag.eeg, filt_b.imag, filt_a.imag);
    % 获取排序后的索引
    [~, sort_indices] = sort(epo.imag.eeg.y(1, :));
    
    % 使用排序后的索引对矩阵和标签矩阵进行重新排列
    epo.imag.eeg.x = epo.imag.eeg.x(:, :, sort_indices);
    epo.imag.eeg.y = epo.imag.eeg.y(:, sort_indices);
    
    %% classification by using moving time windows
    StepSize = 1*1000; % msec
    WindowSize = 3*1000; % msec
    ival_start = (ival_epo(1):StepSize:ival_epo(end)-WindowSize)';
    ival_end = ival_start+WindowSize;
    ival = [ival_start, ival_end];
    nStep = length(ival);%nStep = 33
    
    for stepIdx = 1:nStep
        %clear slope;
        
        %NIRS特征提取
        slope.imag.deoxy{stepIdx} = proc_slopeAcrossTime(epo.imag.deoxy, ival(stepIdx,:));%提取全部的NIRS数据（未做处理）
        slope.imag.oxy{stepIdx}   = proc_slopeAcrossTime(epo.imag.oxy,   ival(stepIdx,:));
        slope.imag.deoxy{stepIdx}.x = (f_FilteredDataMatrixtoFeatures_NIRS(epo.imag.deoxy.x,10,1,[stepIdx-1 stepIdx+2],1,1,1,1,0,0,0))';
        slope.imag.oxy{stepIdx}.x = (f_FilteredDataMatrixtoFeatures_NIRS(epo.imag.oxy.x,10,1,[stepIdx-1 stepIdx+2],1,1,1,1,0,0,0))';
        %                                                                  Mean:求均值;  Variance: 求方差;     Skewness: 偏度;    Kurtosis: 峰度
        
        
        %EEG特征提取
        segment.imag{stepIdx} = proc_selectIval(epo.imag.eeg, ival(stepIdx,:));%提取全部的EEG数据（未做处理）
        segment.imag{stepIdx} = rmfield( segment.imag{stepIdx},'clab');
        segment.imag{stepIdx} = rmfield( segment.imag{stepIdx},'className');
        segment.imag{stepIdx} = rmfield( segment.imag{stepIdx},'t');
        segment.imag{stepIdx} = rmfield( segment.imag{stepIdx},'event');
        segment.imag{stepIdx} = rmfield( segment.imag{stepIdx},'mrk_info');
        segment.imag{stepIdx} = rmfield( segment.imag{stepIdx},'cnt_info');
        segment.imag{stepIdx} = rmfield( segment.imag{stepIdx},'refIval');
        segment.imag{stepIdx}.y = segment.imag{stepIdx}.y(1,:);
        segment.imag{stepIdx}.CSPMatrix = learnCSP(segment.imag{stepIdx});%求协方差矩阵
        segment.imag{stepIdx}.x = extractCSPFeatures(segment.imag{stepIdx}, segment.imag{stepIdx}.CSPMatrix,2)';%特征提取
        
        
        
        %数据混合
        slope.imag.eegoxy{stepIdx}.x = [slope.imag.oxy{stepIdx}.x;segment.imag{stepIdx}.x];%EEG与Oxy混合
        slope.imag.eegoxy{stepIdx}.y = slope.imag.oxy{stepIdx}.y;
        slope.imag.eegdeoxy{stepIdx}.x = [slope.imag.deoxy{stepIdx}.x;segment.imag{stepIdx}.x];%EEG与Deoxy混合
        slope.imag.eegdeoxy{stepIdx}.y = slope.imag.oxy{stepIdx}.y;
        slope.imag.eegdeoxyoxy{stepIdx}.x = [slope.imag.eegdeoxy{stepIdx}.x;slope.imag.oxy{stepIdx}.x];%EEG与Deoxy与Oxy混合
        slope.imag.eegdeoxyoxy{stepIdx}.y = slope.imag.oxy{stepIdx}.y;
        slope.imag.oxydeoxy{stepIdx}.x = [slope.imag.oxy{stepIdx}.x;slope.imag.deoxy{stepIdx}.x];%Deoxy与Oxy混合
        slope.imag.oxydeoxy{stepIdx}.y = slope.imag.oxy{stepIdx}.y;
        slope.imag.eeg{stepIdx}.x = segment.imag{stepIdx}.x;%把eeg数据放到slope里面
        slope.imag.eeg{stepIdx}.y = slope.imag.oxy{stepIdx}.y;
        
        
        % 最大—最小值归一化
        slope.imag.deoxy{stepIdx}.x = mystand(slope.imag.deoxy{stepIdx}.x);
        slope.imag.oxy{stepIdx}.x = mystand(slope.imag.oxy{stepIdx}.x);
        slope.imag.eeg{stepIdx}.x = mystand(slope.imag.eeg{stepIdx}.x);
        slope.imag.eegoxy{stepIdx}.x = mystand(slope.imag.eegoxy{stepIdx}.x);%EEG与Oxy混合
        slope.imag.eegdeoxy{stepIdx}.x = mystand(slope.imag.eegdeoxy{stepIdx}.x);%EEG与Deoxy混合
        slope.imag.eegdeoxyoxy{stepIdx}.x = mystand(slope.imag.eegdeoxyoxy{stepIdx}.x);%EEG与Deoxy与Oxy混合
        slope.imag.oxydeoxy{stepIdx}.x = mystand(slope.imag.oxydeoxy{stepIdx}.x);%Deoxy与Oxy混合
        
        slopebackup = slope;
        
        %------------------------------------------------------------------------------------------------------------------------
        %deoxy训练集逐步回归分析筛选
        [~,~,~,slope.imag.deoxy{stepIdx}.in,~,~,~] = stepwisefit((slope.imag.deoxy{stepIdx}.x)',(slope.imag.deoxy{stepIdx}.y(1,:))','maxiter',500);
        slope.imag.deoxy{stepIdx}.in = slope.imag.deoxy{stepIdx}.in';
        slope.imag.deoxy{stepIdx}.x =slope.imag.deoxy{stepIdx}.x(slope.imag.deoxy{stepIdx}.in,:);
        %deoxy测试集参考deoxy进行筛选
        %slope.imag.deoxy{stepIdx}.x_test = slope.imag.deoxy{stepIdx}.x_test(slope.imag.deoxy{stepIdx}.in,:);
        
        %oxy训练集逐步回归分析筛选
        [~,~,~,slope.imag.oxy{stepIdx}.in,~,~,~] = stepwisefit((slope.imag.oxy{stepIdx}.x)',(slope.imag.oxy{stepIdx}.y(1,:))','maxiter',500);
        slope.imag.oxy{stepIdx}.in = slope.imag.oxy{stepIdx}.in';
        slope.imag.oxy{stepIdx}.x =slope.imag.oxy{stepIdx}.x(slope.imag.oxy{stepIdx}.in,:);
        %oxy测试集参考oxy进行筛选
        %slope.imag.oxy{stepIdx}.x_test = slope.imag.oxy{stepIdx}.x_test(slope.imag.oxy{stepIdx}.in,:);
        
        %eeg训练集逐步回归分析筛选
        [~,~,~,slope.imag.eeg{stepIdx}.in,~,~,~] = stepwisefit((slope.imag.eeg{stepIdx}.x)',(slope.imag.eeg{stepIdx}.y(1,:))','maxiter',500);
        slope.imag.eeg{stepIdx}.in = slope.imag.eeg{stepIdx}.in';
        slope.imag.eeg{stepIdx}.x =slope.imag.eeg{stepIdx}.x(slope.imag.eeg{stepIdx}.in,:);
        %eeg测试集参考eeg进行筛选
        %slope.imag.eeg{stepIdx}.x_test = slope.imag.eeg{stepIdx}.x_test(slope.imag.eeg{stepIdx}.in,:);
        
        %eegoxy训练集逐步回归分析筛选
        [~,~,~,slope.imag.eegoxy{stepIdx}.in,~,~,~] = stepwisefit((slope.imag.eegoxy{stepIdx}.x)',(slope.imag.eegoxy{stepIdx}.y(1,:))','maxiter',500);
        slope.imag.eegoxy{stepIdx}.in = slope.imag.eegoxy{stepIdx}.in';
        slope.imag.eegoxy{stepIdx}.x =slope.imag.eegoxy{stepIdx}.x(slope.imag.eegoxy{stepIdx}.in,:);
        %eegoxy测试集参考eegoxy进行筛选
        %slope.imag.eegoxy{stepIdx}.x_test = slope.imag.eegoxy{stepIdx}.x_test(slope.imag.eegoxy{stepIdx}.in,:);
        
        %eegdeoxy训练集逐步回归分析筛选
        [~,~,~,slope.imag.eegdeoxy{stepIdx}.in,~,~,~] = stepwisefit((slope.imag.eegdeoxy{stepIdx}.x)',(slope.imag.eegdeoxy{stepIdx}.y(1,:))','maxiter',500);
        slope.imag.eegdeoxy{stepIdx}.in = slope.imag.eegdeoxy{stepIdx}.in';
        slope.imag.eegdeoxy{stepIdx}.x =slope.imag.eegdeoxy{stepIdx}.x(slope.imag.eegdeoxy{stepIdx}.in,:);
        %eegdeoxy测试集参考eegdeoxy进行筛选
        %slope.imag.eegdeoxy{stepIdx}.x_test = slope.imag.eegdeoxy{stepIdx}.x_test(slope.imag.eegdeoxy{stepIdx}.in,:);
        
        %eegdeoxyoxy训练集逐步回归分析筛选
        [~,~,~,slope.imag.eegdeoxyoxy{stepIdx}.in,~,~,~] = stepwisefit((slope.imag.eegdeoxyoxy{stepIdx}.x)',(slope.imag.eegdeoxyoxy{stepIdx}.y(1,:))','maxiter',500);
        slope.imag.eegdeoxyoxy{stepIdx}.in = slope.imag.eegdeoxyoxy{stepIdx}.in';
        slope.imag.eegdeoxyoxy{stepIdx}.x =slope.imag.eegdeoxyoxy{stepIdx}.x(slope.imag.eegdeoxyoxy{stepIdx}.in,:);
        %eegdeoxyoxy测试集参考eegdeoxyoxy进行筛选
        %slope.imag.eegdeoxyoxy{stepIdx}.x_test = slope.imag.eegdeoxyoxy{stepIdx}.x_test(slope.imag.eegdeoxyoxy{stepIdx}.in,:);
        
        %oxydeoxy训练集逐步回归分析筛选
        [~,~,~,slope.imag.oxydeoxy{stepIdx}.in,~,~,~] = stepwisefit((slope.imag.oxydeoxy{stepIdx}.x)',(slope.imag.oxydeoxy{stepIdx}.y(1,:))','maxiter',500);
        slope.imag.oxydeoxy{stepIdx}.in = slope.imag.oxydeoxy{stepIdx}.in';
        slope.imag.oxydeoxy{stepIdx}.x =slope.imag.oxydeoxy{stepIdx}.x(slope.imag.oxydeoxy{stepIdx}.in,:);
        %oxydeoxy测试集参考oxydeoxy进行筛选
        %slope.imag.oxydeoxy{stepIdx}.x_test = slope.imag.oxydeoxy{stepIdx}.x_test(slope.imag.oxydeoxy{stepIdx}.in,:);
        
        if(isempty(slope.imag.deoxy{stepIdx}.x) )%|| isempty(slope.imag.deoxy{stepIdx}.x_test))
            slope.imag.deoxy{stepIdx}.x = slopebackup.imag.deoxy{stepIdx}.x;%保留全部训练集
            %slope.imag.deoxy{stepIdx}.x_test = B1;%保留全部测试集
        end
        if(isempty(slope.imag.oxy{stepIdx}.x) )%|| isempty(slope.imag.oxy{stepIdx}.x_test))
            slope.imag.oxy{stepIdx}.x = slopebackup.imag.oxy{stepIdx}.x;%保留全部训练集
            %slope.imag.oxy{stepIdx}.x_test = B2;%保留全部测试集
        end
        if(isempty(slope.imag.eeg{stepIdx}.x) )%|| isempty(slope.imag.eeg{stepIdx}.x_test))
            slope.imag.eeg{stepIdx}.x = slopebackup.imag.eeg{stepIdx}.x;%保留全部训练集
            %slope.imag.eeg{stepIdx}.x_test = B3;%保留全部测试集
        end
        if(isempty(slope.imag.eegoxy{stepIdx}.x))% || isempty(slope.imag.eegoxy{stepIdx}.x_test))
            slope.imag.eegoxy{stepIdx}.x = slopebackup.imag.eegoxy{stepIdx}.x;%保留全部训练集
            %slope.imag.eegoxy{stepIdx}.x_test = B4;%保留全部测试集
        end
        if(isempty(slope.imag.eegdeoxy{stepIdx}.x) )%|| isempty(slope.imag.eegdeoxy{stepIdx}.x_test))
            slope.imag.eegdeoxy{stepIdx}.x = slopebackup.imag.eegdeoxy{stepIdx}.x;%保留全部训练集
            %slope.imag.eegdeoxy{stepIdx}.x_test = B5;%保留全部测试集
        end
        if(isempty(slope.imag.eegdeoxyoxy{stepIdx}.x))% || isempty(slope.imag.eegdeoxyoxy{stepIdx}.x_test))
            slope.imag.eegdeoxyoxy{stepIdx}.x = slopebackup.imag.eegdeoxyoxy{stepIdx}.x;%保留全部训练集
            %slope.imag.eegdeoxyoxy{stepIdx}.x_test = B6;%保留全部测试集
        end
        if(isempty(slope.imag.oxydeoxy{stepIdx}.x) )%|| isempty(slope.imag.oxydeoxy{stepIdx}.x_test))
            slope.imag.oxydeoxy{stepIdx}.x = slopebackup.imag.oxydeoxy{stepIdx}.x;%保留全部训练集
            %slope.imag.oxydeoxy{stepIdx}.x_test = B7;%保留全部测试集
        end
        
        %fprintf('Motor imagery, Step: %d/%d，ID：%d/%d\n', stepIdx, nStep,id,subNum);
        
        
        %每组提取10%的盲源信号
        bilv = 1/6;
        
        % deoxy随机取10%作为最后的训练集
        [M1,N1]=size(slope.imag.deoxy{stepIdx}.x'); %读取矩阵行列数;
        A1 = slope.imag.deoxy{stepIdx}.x';
        Y1 = slope.imag.deoxy{stepIdx}.y';
        num1 = round(M1*bilv); % 取A的1/3行作为盲源特征，round是让里面的数字为四舍五入取整;
        [~,idx]=sort(rand(M1,1));%随机排列生成index;
        B1=(A1(idx(1:num1),:))';%测试集B;
        Y1_1=(Y1(idx(1:num1),:))';%测试集B标签
        C1=A1(idx(num1+1:M1),:);%训练集C;
        Y1_2=(Y1(idx(num1+1:M1),:))';%训练集C标签
        slope.imag.deoxy{stepIdx}.x = C1';%保存训练数据
        slope.imag.deoxy{stepIdx}.y = Y1_2;
        slope.imag.deoxy{stepIdx}.x_test = B1;
        slope.imag.deoxy{stepIdx}.y_test = Y1_1;
        
        % oxy随机取10%作为最后的训练集
        [M2,N2]=size(slope.imag.oxy{stepIdx}.x'); %读取矩阵行列数;
        A2 = slope.imag.oxy{stepIdx}.x';
        num2 = round(M2*bilv); % 取A的1/3行作为训练集，round为四舍五入取整;
        [~,idx]=sort(rand(M2,1));%随机排列生成index;
        B2=(A2(idx(1:num2),:))';%测试集B2;
        Y2_1=(Y1(idx(1:num2),:))';%测试集B2标签
        C2=A2(idx(num2+1:M2),:);%训练集C2;
        Y2_2=(Y1(idx(num2+1:M2),:))';%训练集C2标签
        slope.imag.oxy{stepIdx}.x = C2';%保存训练数据
        slope.imag.oxy{stepIdx}.y = Y2_2;
        slope.imag.oxy{stepIdx}.x_test = B2;
        slope.imag.oxy{stepIdx}.y_test = Y2_1;
        
        % eeg随机取10%作为最后的训练集
        [M3,N3]=size(slope.imag.eeg{stepIdx}.x'); %读取矩阵行列数;
        A3 = slope.imag.eeg{stepIdx}.x';
        num3 = round(M3*bilv); % 取A的1/3行作为训练集，round为四舍五入取整;
        [~,idx]=sort(rand(M3,1));%随机排列生成index;
        B3=(A3(idx(1:num3),:))';%测试集B2;
        Y3_1=(Y1(idx(1:num3),:))';%测试集B2标签
        C3=A3(idx(num3+1:M3),:);%训练集C2;
        Y3_2=(Y1(idx(num3+1:M3),:))';%训练集C2标签
        slope.imag.eeg{stepIdx}.x = C3';%保存训练数据
        slope.imag.eeg{stepIdx}.y = Y3_2;
        slope.imag.eeg{stepIdx}.x_test = B3;
        slope.imag.eeg{stepIdx}.y_test = Y3_1;
        
        % eegoxy随机取10%作为最后的训练集
        [M4,N4]=size(slope.imag.eegoxy{stepIdx}.x'); %读取矩阵行列数;
        A4 = slope.imag.eegoxy{stepIdx}.x';
        num4 = round(M4*bilv); % 取A的1/3行作为训练集，round为四舍五入取整;
        [~,idx]=sort(rand(M4,1));%随机排列生成index;
        B4=(A4(idx(1:num4),:))';%测试集B2;
        Y4_1=(Y1(idx(1:num4),:))';%测试集B2标签
        C4=A4(idx(num4+1:M4),:);%训练集C2;
        Y4_2=(Y1(idx(num4+1:M4),:))';%训练集C2标签
        slope.imag.eegoxy{stepIdx}.x = C4';%保存训练数据
        slope.imag.eegoxy{stepIdx}.y = Y4_2;
        slope.imag.eegoxy{stepIdx}.x_test = B4;
        slope.imag.eegoxy{stepIdx}.y_test = Y4_1;
        
        % eegdeoxy随机取10%作为最后的训练集
        [M5,N5]=size(slope.imag.eegdeoxy{stepIdx}.x'); %读取矩阵行列数;
        A5 = slope.imag.eegdeoxy{stepIdx}.x';
        num5 = round(M5*bilv); % 取A的1/3行作为训练集，round是让里面的数字为四舍五入取整;
        [~,idx]=sort(rand(M5,1));%随机排列生成index;
        B5=(A5(idx(1:num5),:))';%测试集B;
        Y5_1=(Y1(idx(1:num5),:))';%测试集B标签
        C5=A5(idx(num5+1:M5),:);%训练集C;
        Y5_2=(Y1(idx(num5+1:M5),:))';%训练集C标签
        slope.imag.eegdeoxy{stepIdx}.x = C5';%保存训练数据
        slope.imag.eegdeoxy{stepIdx}.y = Y5_2;
        slope.imag.eegdeoxy{stepIdx}.x_test = B5;
        slope.imag.eegdeoxy{stepIdx}.y_test = Y5_1;
        
        % eegdeoxyoxy随机取10%作为最后的训练集
        [M6,N6]=size(slope.imag.eegdeoxyoxy{stepIdx}.x'); %读取矩阵行列数;
        A6 = slope.imag.eegdeoxyoxy{stepIdx}.x';
        num6 = round(M6*bilv); % 取A的1/3行作为训练集，round是让里面的数字为四舍五入取整;
        [~,idx]=sort(rand(M6,1));%随机排列生成index;
        B6=(A6(idx(1:num6),:))';%测试集B;
        Y6_1=(Y1(idx(1:num6),:))';%测试集B标签
        C6=A6(idx(num6+1:M6),:);%训练集C;
        Y6_2=(Y1(idx(num6+1:M6),:))';%训练集C标签
        slope.imag.eegdeoxyoxy{stepIdx}.x = C6';%保存训练数据
        slope.imag.eegdeoxyoxy{stepIdx}.y = Y6_2;
        slope.imag.eegdeoxyoxy{stepIdx}.x_test = B6;
        slope.imag.eegdeoxyoxy{stepIdx}.y_test = Y6_1;
        
        % oxydeoxy随机取10%作为最后的训练集
        [M7,N7]=size(slope.imag.oxydeoxy{stepIdx}.x'); %读取矩阵行列数;
        A7 = slope.imag.oxydeoxy{stepIdx}.x';
        num7 = round(M7*bilv); % 取A的1/3行作为训练集，round是让里面的数字为四舍五入取整;
        [~,idx]=sort(rand(M7,1));%随机排列生成index;
        B7=(A7(idx(1:num7),:))';%测试集B;
        Y7_1=(Y1(idx(1:num7),:))';%测试集B标签
        C7=A7(idx(num7+1:M7),:);%训练集C;
        Y7_2=(Y1(idx(num7+1:M7),:))';%训练集C标签
        slope.imag.oxydeoxy{stepIdx}.x = C7';%保存训练数据
        slope.imag.oxydeoxy{stepIdx}.y = Y7_2;
        slope.imag.oxydeoxy{stepIdx}.x_test = B7;
        slope.imag.oxydeoxy{stepIdx}.y_test = Y7_1;
    end
    
    %slopebackup = slope;
    
    %% 交叉验证参数，Shift为重复次数，Fold为几折
    nShift = 1;
    nFold = 5;
    datasize = 60 - 60*bilv;
    %group.imag = slope.imag.deoxy{1,1}.y; % epo.imag.deoxy.y == epo.imag.oxy.y == epo.imag.eeg.y
    epobackup = epo;
    clear stepIdx;
    
    
    
    % 混合部分
    for shiftIdx = 1:nShift
        indices.imag{shiftIdx} = crossvalind('Kfold',datasize,nFold);
        for stepIdx = 1:nStep
            for foldIdx = 1:nFold
                fprintf('Motor imagery, Repeat: %d/%d, Step: %d/%d，Fold: %d/%d, ID：%d/%d\n',shiftIdx, nShift, stepIdx, nStep, foldIdx, nFold, id, subNum);
                
                test = (indices.imag{shiftIdx} == foldIdx);
                train = ~test;
                
                
                %十折交叉数据划分
                
                [hang,lie] = size(slope.imag.eeg{stepIdx}.x);
                if hang > lie
                    slope.imag.eeg{stepIdx}.x = slope.imag.eeg{stepIdx}.x';
                end
                
                fv_train.deoxy.x = (slope.imag.deoxy{stepIdx}.x(:,train))';
                fv_train.oxy.x = (slope.imag.oxy{stepIdx}.x(:,train))';
                fv_train.eeg.x = (slope.imag.eeg{stepIdx}.x(:,train))';
                fv_train.eegoxy.x = (slope.imag.eegoxy{stepIdx}.x(:,train))';
                fv_train.eegdeoxy.x = (slope.imag.eegdeoxy{stepIdx}.x(:,train))';
                fv_train.eegdeoxyoxy.x = (slope.imag.eegdeoxyoxy{stepIdx}.x(:,train))';
                fv_train.oxydeoxy.x = (slope.imag.oxydeoxy{stepIdx}.x(:,train))';
                
                fv_test.deoxy.x = (slope.imag.deoxy{stepIdx}.x(:,test))';
                fv_test.oxy.x = (slope.imag.oxy{stepIdx}.x(:,test))';
                fv_test.eeg.x = (slope.imag.eeg{stepIdx}.x(:,test))';
                fv_test.eegoxy.x = (slope.imag.eegoxy{stepIdx}.x(:,test))';
                fv_test.eegdeoxy.x = (slope.imag.eegdeoxy{stepIdx}.x(:,test))';
                fv_test.eegdeoxyoxy.x = (slope.imag.eegdeoxyoxy{stepIdx}.x(:,test))';
                fv_test.oxydeoxy.x = (slope.imag.oxydeoxy{stepIdx}.x(:,test))';
                
                y_train1 = (slope.imag.deoxy{stepIdx}.y(1,train))';
                y_train2 = (slope.imag.oxy{stepIdx}.y(1,train))';
                y_train3 = (slope.imag.eeg{stepIdx}.y(1,train))';
                y_train4 = (slope.imag.eegoxy{stepIdx}.y(1,train))';
                y_train5 = (slope.imag.eegdeoxy{stepIdx}.y(1,train))';
                y_train6 = (slope.imag.eegdeoxyoxy{stepIdx}.y(1,train))';
                y_train7 = (slope.imag.oxydeoxy{stepIdx}.y(1,train))';
                
                y_test1  = (slope.imag.deoxy{stepIdx}.y(1,test))';
                y_test2  = (slope.imag.oxy{stepIdx}.y(1,test))';
                y_test3 = (slope.imag.eeg{stepIdx}.y(1,test))';
                y_test4 = (slope.imag.eegoxy{stepIdx}.y(1,test))';
                y_test5 = (slope.imag.eegdeoxy{stepIdx}.y(1,test))';
                y_test6 = (slope.imag.eegdeoxyoxy{stepIdx}.y(1,test))';
                y_test7 = (slope.imag.oxydeoxy{stepIdx}.y(1,test))';
                
                C{stepIdx,foldIdx,id}.deoxy = fitcdiscr(fv_train.deoxy.x,y_train1,'DiscrimType','pseudoLinear');
                C{stepIdx,foldIdx,id}.oxy   = fitcdiscr(fv_train.oxy.x,y_train2,'DiscrimType','pseudoLinear');
                C{stepIdx,foldIdx,id}.eeg   = fitcdiscr(fv_train.eeg.x,y_train3,'DiscrimType','pseudoLinear');
                C{stepIdx,foldIdx,id}.eegoxy   = fitcdiscr(fv_train.eegoxy.x,y_train4,'DiscrimType','pseudoLinear');
                C{stepIdx,foldIdx,id}.eegdeoxy   = fitcdiscr(fv_train.eegdeoxy.x,y_train5,'DiscrimType','pseudoLinear');
                C{stepIdx,foldIdx,id}.eegdeoxyoxy   = fitcdiscr(fv_train.eegdeoxyoxy.x,y_train6,'DiscrimType','pseudoLinear');
                C{stepIdx,foldIdx,id}.oxydeoxy   = fitcdiscr(fv_train.oxydeoxy.x,y_train7,'DiscrimType','pseudoLinear');
                
                % classification
                grouphat.deoxy(foldIdx,:)  = predict(C{stepIdx,foldIdx,id}.deoxy,fv_test.deoxy.x);
                grouphat.oxy(foldIdx,:)    = predict(C{stepIdx,foldIdx,id}.oxy,  fv_test.oxy.x);
                grouphat.eeg(foldIdx,:)    = predict(C{stepIdx,foldIdx,id}.eeg,  fv_test.eeg.x);
                grouphat.eegoxy(foldIdx,:)  = predict(C{stepIdx,foldIdx,id}.eegoxy, fv_test.eegoxy.x);
                grouphat.eegdeoxy(foldIdx,:)  = predict(C{stepIdx,foldIdx,id}.eegdeoxy, fv_test.eegdeoxy.x);
                grouphat.eegdeoxyoxy(foldIdx,:)  = predict(C{stepIdx,foldIdx,id}.eegdeoxyoxy, fv_test.eegdeoxyoxy.x);
                grouphat.oxydeoxy(foldIdx,:)  = predict(C{stepIdx,foldIdx,id}.oxydeoxy, fv_test.oxydeoxy.x);
                
                cmat.deoxy(:,:,foldIdx)  = confusionmat(y_test1, grouphat.deoxy(foldIdx,:));
                cmat.oxy(:,:,foldIdx)    = confusionmat(y_test2, grouphat.oxy(foldIdx,:));
                cmat.eeg(:,:,foldIdx)    = confusionmat(y_test3, grouphat.eeg(foldIdx,:));
                cmat.eegoxy(:,:,foldIdx)  = confusionmat(y_test4, grouphat.eegoxy(foldIdx,:));
                cmat.eegdeoxy(:,:,foldIdx)  = confusionmat(y_test5, grouphat.eegdeoxy(foldIdx,:));
                cmat.eegdeoxyoxy(:,:,foldIdx)  = confusionmat(y_test6, grouphat.eegdeoxyoxy(foldIdx,:));
                cmat.oxydeoxy(:,:,foldIdx)  = confusionmat(y_test7, grouphat.oxydeoxy(foldIdx,:));
                
                
                acc.imag.deoxy(stepIdx,foldIdx,id)    = trace((sum(cmat.deoxy,3))) / sum(sum(sum(cmat.deoxy,3),2),1);
                acc.imag.oxy(stepIdx,foldIdx,id)      = trace((sum(cmat.oxy,3)))   / sum(sum(sum(cmat.oxy,3),2),1);
                acc.imag.eeg(stepIdx,foldIdx,id)      = trace((sum(cmat.eeg,3)))   / sum(sum(sum(cmat.eeg,3),2),1);
                acc.imag.eegoxy(stepIdx,foldIdx,id)    = trace((sum(cmat.eegoxy,3)))  / sum(sum(sum(cmat.eegoxy,3),2),1);
                acc.imag.eegdeoxy(stepIdx,foldIdx,id)    = trace((sum(cmat.eegdeoxy,3)))  / sum(sum(sum(cmat.eegdeoxy,3),2),1);
                acc.imag.eegdeoxyoxy(stepIdx,foldIdx,id)    = trace((sum(cmat.eegdeoxyoxy,3)))  / sum(sum(sum(cmat.eegdeoxyoxy,3),2),1);
                acc.imag.oxydeoxy(stepIdx,foldIdx,id)    = trace((sum(cmat.oxydeoxy,3)))  / sum(sum(sum(cmat.oxydeoxy,3),2),1);
                %slopebackup = slope;
            end
        end
    end
    
    
    
    mean_acc.imag.deoxy(:,id)  = mean(squeeze(acc.imag.deoxy(:,:,id))',1);
    colName = sprintf('Subject %d', id);
    T_deoxy.(colName) = mean_acc.imag.deoxy(:, id);
    
    mean_acc.imag.oxy(:,id)    = mean(squeeze(acc.imag.oxy(:,:,id))',1);
    colName = sprintf('Subject %d', id);
    T_oxy.(colName) = mean_acc.imag.oxy(:, id);
    
    mean_acc.imag.eeg(:,id)    = mean(squeeze(acc.imag.eeg(:,:,id))',1);
    colName = sprintf('Subject %d', id);
    T_eeg.(colName) = mean_acc.imag.eeg(:, id);
    
    mean_acc.imag.eegoxy(:,id)  = mean(squeeze(acc.imag.eegoxy(:,:,id))',1);
    colName = sprintf('Subject %d', id);
    T_eegoxy.(colName) = mean_acc.imag.eegoxy(:, id);
    
    mean_acc.imag.eegdeoxy(:,id)  = mean(squeeze(acc.imag.eegdeoxy(:,:,id))',1);
    colName = sprintf('Subject %d', id);
    T_eegdeoxy.(colName) = mean_acc.imag.eegdeoxy(:, id);
    
    mean_acc.imag.eegdeoxyoxy(:,id)  = mean(squeeze(acc.imag.eegdeoxyoxy(:,:,id))',1);
    colName = sprintf('Subject %d', id);
    T_eegdeoxyoxy.(colName) = mean_acc.imag.eegdeoxyoxy(:, id);
    
    mean_acc.imag.oxydeoxy(:,id)  = mean(squeeze(acc.imag.oxydeoxy(:,:,id))',1);
    colName = sprintf('Subject %d', id);
    T_oxydeoxy.(colName) = mean_acc.imag.oxydeoxy(:, id);
    slopeall{1,id} = slope;
    clear slope;
end
T_deoxy = table2array (T_deoxy);
T_oxy = table2array (T_oxy);
T_eeg = table2array (T_eeg);
T_eegoxy   = table2array (T_eegoxy);
T_eegdeoxy = table2array (T_eegdeoxy);
T_eegdeoxyoxy   = table2array (T_eegdeoxyoxy);
T_oxydeoxy = table2array (T_oxydeoxy);

[Y_row_T_deoxy,Ind_col_T_deoxy]=max(T_deoxy); %每列的最大值及列号
[Y_row_T_oxy,Ind_col_T_oxy]=max(T_oxy); %每列的最大值及列号
[Y_row_T_eeg,Ind_col_T_eeg]=max(T_eeg); %每列的最大值及列号
[Y_row_T_eegoxy,Ind_col_T_eegoxy]=max(T_eegoxy); %每列的最大值及列号
[Y_row_T_eegdeoxy,Ind_col_T_eegdeoxy]=max(T_eegdeoxy); %每列的最大值及列号
[Y_row_T_eegdeoxyoxy,Ind_col_T_eegdeoxyoxy]=max(T_eegdeoxyoxy); %每列的最大值及列号
[Y_row_T_oxydeoxy,Ind_col_T_oxydeoxy]=max(T_oxydeoxy); %每列的最大值及列号

mean(Y_row_T_deoxy)
mean(Y_row_T_oxy)
mean(Y_row_T_eeg)
mean(Y_row_T_eegoxy)
mean(Y_row_T_eegdeoxy)
mean(Y_row_T_eegdeoxyoxy)
mean(Y_row_T_oxydeoxy)



%寻找每个个体的7种模态集合的对应最优窗口后，找到N折所对应的准确度最大值，和其所对应的分类器，将盲源数据通过分类器分类预测
%每个个体的33个窗口都要用最优折的分类器去预测盲源数据
for id = 1:29
    
    %部分结构初始化
    grouphat_my = cell(1,33);
    cmat_my     = cell(1,33);
    
    for winNum = 1:33
        
        %EEG部分
        %winID.eeg = Ind_col_T_eeg(1,id);                    %提取对应的最优窗口索引到winID结构体的eeg
        fold_acc.eeg = acc.imag.eeg(winNum,:,id);        %根据winNum找到准确度集合并存到fold_acc的eeg
        [~,foldID.eeg]=max(fold_acc.eeg);                   %取出N折准确度集合的最大值索引
        C_Best.eeg = C{winNum,foldID.eeg,id}.eeg;        %取出N折准确度集合最大值对应的分类器并存放到C_Best.eeg
        %盲源数据预测
        grouphat_my{id,winNum}.eeg    = predict(C_Best.eeg,(slopeall{1,id}.imag.eeg{1,winNum}.x_test)');
        cmat_my{id,winNum}.eeg        = confusionmat((slopeall{1,id }.imag.eeg{1,winNum}.y_test(1,:))', grouphat_my{id,winNum}.eeg);
        acc_my.eeg{id,winNum}         = trace((sum(cmat_my{id,winNum}.eeg,3)))/sum(sum(sum(cmat_my{id,winNum}.eeg,3),2),1);
        
        
        %OXY部分
        fold_acc.oxy = acc.imag.oxy(winNum,:,id);        %根据winNum找到准确度集合并存到fold_acc的eeg
        [~,foldID.oxy]=max(fold_acc.oxy);                   %取出N折准确度集合的最大值索引
        C_Best.oxy = C{winNum,foldID.oxy,id}.oxy;        %取出N折准确度集合最大值对应的分类器并存放到C_Best.eeg
        %盲源数据预测
        grouphat_my{id,winNum}.oxy    = predict(C_Best.oxy,(slopeall{1,id}.imag.oxy{1,winNum}.x_test)');
        cmat_my{id,winNum}.oxy        = confusionmat((slopeall{1,id }.imag.oxy{1,winNum}.y_test(1,:))', grouphat_my{id,winNum}.oxy);
        acc_my.oxy{id,winNum}         = trace((sum(cmat_my{id,winNum}.oxy,3)))/sum(sum(sum(cmat_my{id,winNum}.oxy,3),2),1);
        
        %DEOXY部分
        fold_acc.deoxy = acc.imag.deoxy(winNum,:,id);        %根据winNum找到准确度集合并存到fold_acc的eeg
        [~,foldID.deoxy]=max(fold_acc.deoxy);                   %取出N折准确度集合的最大值索引
        C_Best.deoxy = C{winNum,foldID.deoxy,id}.deoxy;        %取出N折准确度集合最大值对应的分类器并存放到C_Best.eeg
        %盲源数据预测
        grouphat_my{id,winNum}.deoxy    = predict(C_Best.deoxy,(slopeall{1,id}.imag.deoxy{1,winNum}.x_test)');
        cmat_my{id,winNum}.deoxy        = confusionmat((slopeall{1,id }.imag.deoxy{1,winNum}.y_test(1,:))', grouphat_my{id,winNum}.deoxy);
        acc_my.deoxy{id,winNum}         = trace((sum(cmat_my{id,winNum}.deoxy,3)))/sum(sum(sum(cmat_my{id,winNum}.deoxy,3),2),1);
        
        %EEGOXY部分
        fold_acc.eegoxy = acc.imag.eegoxy(winNum,:,id);        %根据winNum找到准确度集合并存到fold_acc的eeg
        [~,foldID.eegoxy]=max(fold_acc.eegoxy);                   %取出N折准确度集合的最大值索引
        C_Best.eegoxy = C{winNum,foldID.eegoxy,id}.eegoxy;        %取出N折准确度集合最大值对应的分类器并存放到C_Best.eeg
        %盲源数据预测
        grouphat_my{id,winNum}.eegoxy    = predict(C_Best.eegoxy,(slopeall{1,id}.imag.eegoxy{1,winNum}.x_test)');
        cmat_my{id,winNum}.eegoxy        = confusionmat((slopeall{1,id }.imag.eegoxy{1,winNum}.y_test(1,:))', grouphat_my{id,winNum}.eegoxy);
        acc_my.eegoxy{id,winNum}         = trace((sum(cmat_my{id,winNum}.eegoxy,3)))/sum(sum(sum(cmat_my{id,winNum}.eegoxy,3),2),1);        
        
        %EEGDEOXY部分
        fold_acc.eegdeoxy = acc.imag.eegdeoxy(winNum,:,id);        %根据winNum找到准确度集合并存到fold_acc的eeg
        [~,foldID.eegdeoxy]=max(fold_acc.eegdeoxy);                   %取出N折准确度集合的最大值索引
        C_Best.eegdeoxy = C{winNum,foldID.eegdeoxy,id}.eegdeoxy;        %取出N折准确度集合最大值对应的分类器并存放到C_Best.eeg
        %盲源数据预测
        grouphat_my{id,winNum}.eegdeoxy    = predict(C_Best.eegdeoxy,(slopeall{1,id}.imag.eegdeoxy{1,winNum}.x_test)');
        cmat_my{id,winNum}.eegdeoxy        = confusionmat((slopeall{1,id }.imag.eegdeoxy{1,winNum}.y_test(1,:))', grouphat_my{id,winNum}.eegdeoxy);
        acc_my.eegdeoxy{id,winNum}         = trace((sum(cmat_my{id,winNum}.eegdeoxy,3)))/sum(sum(sum(cmat_my{id,winNum}.eegdeoxy,3),2),1);        

        %EEGDEOXYOXY部分
        fold_acc.eegdeoxyoxy = acc.imag.eegdeoxyoxy(winNum,:,id);        %根据winNum找到准确度集合并存到fold_acc的eeg
        [~,foldID.eegdeoxyoxy]=max(fold_acc.eegdeoxyoxy);                   %取出N折准确度集合的最大值索引
        C_Best.eegdeoxyoxy = C{winNum,foldID.eegdeoxyoxy,id}.eegdeoxyoxy;        %取出N折准确度集合最大值对应的分类器并存放到C_Best.eeg
        %盲源数据预测
        grouphat_my{id,winNum}.eegdeoxyoxy    = predict(C_Best.eegdeoxyoxy,(slopeall{1,id}.imag.eegdeoxyoxy{1,winNum}.x_test)');
        cmat_my{id,winNum}.eegdeoxyoxy        = confusionmat((slopeall{1,id }.imag.eegdeoxyoxy{1,winNum}.y_test(1,:))', grouphat_my{id,winNum}.eegdeoxyoxy);
        acc_my.eegdeoxyoxy{id,winNum}         = trace((sum(cmat_my{id,winNum}.eegdeoxyoxy,3)))/sum(sum(sum(cmat_my{id,winNum}.eegdeoxyoxy,3),2),1);        
 
        %DEOXYOXY部分
        fold_acc.oxydeoxy = acc.imag.oxydeoxy(winNum,:,id);        %根据winNum找到准确度集合并存到fold_acc的eeg
        [~,foldID.oxydeoxy]=max(fold_acc.oxydeoxy);                   %取出N折准确度集合的最大值索引
        C_Best.oxydeoxy = C{winNum,foldID.oxydeoxy,id}.oxydeoxy;        %取出N折准确度集合最大值对应的分类器并存放到C_Best.eeg
        %盲源数据预测
        grouphat_my{id,winNum}.oxydeoxy    = predict(C_Best.oxydeoxy,(slopeall{1,id}.imag.oxydeoxy{1,winNum}.x_test)');
        cmat_my{id,winNum}.oxydeoxy        = confusionmat((slopeall{1,id }.imag.oxydeoxy{1,winNum}.y_test(1,:))', grouphat_my{id,winNum}.oxydeoxy);
        acc_my.oxydeoxy{id,winNum}         = trace((sum(cmat_my{id,winNum}.oxydeoxy,3)))/sum(sum(sum(cmat_my{id,winNum}.oxydeoxy,3),2),1);        

        clear fold_acc foldID C_Best
    end
    
end



% a = size(slope.imag.deoxy{1,1}.x,1)+size(slope.imag.deoxy{1,2}.x,1)+size(slope.imag.deoxy{1,3}.x,1)+size(slope.imag.deoxy{1,4}.x,1)+size(slope.imag.deoxy{1,5}.x,1)+size(slope.imag.deoxy{1,6}.x,1)+size(slope.imag.deoxy{1,7}.x,1)+size(slope.imag.deoxy{1,8}.x,1)+size(slope.imag.deoxy{1,9}.x,1)+size(slope.imag.deoxy{1,10}.x,1)+size(slope.imag.deoxy{1,11}.x,1)+size(slope.imag.deoxy{1,12}.x,1)+size(slope.imag.deoxy{1,13}.x,1)+size(slope.imag.deoxy{1,14}.x,1)+size(slope.imag.deoxy{1,15}.x,1)+size(slope.imag.deoxy{1,16}.x,1)+size(slope.imag.deoxy{1,17}.x,1)+size(slope.imag.deoxy{1,18}.x,1)+size(slope.imag.deoxy{1,19}.x,1)+size(slope.imag.deoxy{1,20}.x,1)+size(slope.imag.deoxy{1,21}.x,1)+size(slope.imag.deoxy{1,22}.x,1)+size(slope.imag.deoxy{1,23}.x,1)+size(slope.imag.deoxy{1,24}.x,1)+size(slope.imag.deoxy{1,25}.x,1)+size(slope.imag.deoxy{1,26}.x,1)+size(slope.imag.deoxy{1,27}.x,1)+size(slope.imag.deoxy{1,28}.x,1)+size(slope.imag.deoxy{1,29}.x,1)+size(slope.imag.deoxy{1,30}.x,1)+size(slope.imag.deoxy{1,31}.x,1)+size(slope.imag.deoxy{1,32}.x,1)+size(slope.imag.deoxy{1,33}.x,1)
% b = size(slope.imag.oxy{1,1}.x,1)+size(slope.imag.oxy{1,2}.x,1)+size(slope.imag.oxy{1,3}.x,1)+size(slope.imag.oxy{1,4}.x,1)+size(slope.imag.oxy{1,5}.x,1)+size(slope.imag.oxy{1,6}.x,1)+size(slope.imag.oxy{1,7}.x,1)+size(slope.imag.oxy{1,8}.x,1)+size(slope.imag.oxy{1,9}.x,1)+size(slope.imag.oxy{1,10}.x,1)+size(slope.imag.oxy{1,11}.x,1)+size(slope.imag.oxy{1,12}.x,1)+size(slope.imag.oxy{1,13}.x,1)+size(slope.imag.oxy{1,14}.x,1)+size(slope.imag.oxy{1,15}.x,1)+size(slope.imag.oxy{1,16}.x,1)+size(slope.imag.oxy{1,17}.x,1)+size(slope.imag.oxy{1,18}.x,1)+size(slope.imag.oxy{1,19}.x,1)+size(slope.imag.oxy{1,20}.x,1)+size(slope.imag.oxy{1,21}.x,1)+size(slope.imag.oxy{1,22}.x,1)+size(slope.imag.oxy{1,23}.x,1)+size(slope.imag.oxy{1,24}.x,1)+size(slope.imag.oxy{1,25}.x,1)+size(slope.imag.oxy{1,26}.x,1)+size(slope.imag.oxy{1,27}.x,1)+size(slope.imag.oxy{1,28}.x,1)+size(slope.imag.oxy{1,29}.x,1)+size(slope.imag.oxy{1,30}.x,1)+size(slope.imag.oxy{1,31}.x,1)+size(slope.imag.oxy{1,32}.x,1)+size(slope.imag.oxy{1,33}.x,1)

