close all
clear integratedGrayValues

[files,paths] = uigetfile('*.czi','SELECT MASTER IMAGES','SELECT MASTER IMAGES','MultiSelect','on');
disp('Master Images Selected')
DESTINATION_FOLDER = uigetdir('C:\','SELECT DESTINATION FOLDER');
disp('Destination folder selected')
filescell = fullfile(paths,files);
%%%%Cell Channel Stitching
for i=1:length(filescell)
Datai = ReadImage6D(char(filescell(i)));
image6d = Datai{1};
squeezei = squeeze(image6d(1,1,:,1,:,:));

sumz = [];
for k=1:size(squeezei,1)
    sumz = [sumz,sum(squeezei(k,:,:),'all')];
end
    
    
MidAi = find(sumz == max(sumz));    
Ai = squeezei(MidAi,:,:);
varname=['A',num2str(i)];
assignin('base',varname,squeezei(MidAi,:,:)); 
%% 

Abz = sum(StitchedvolumeAi(:,:,MidAi),3);

varname=['Abz',num2str(i)];
assignin('base',varname,Abz);

%%%%Nuclei Channel Stitching
Datai = ReadImage6D(char(filescell(i)));
image6d = Datai{1};
squeezei = squeeze(image6d(1,1,:,2,:,:));

sumzb = [];
for k=1:size(squeezei,1)
    sumzb = [sumzb,sum(squeezei,'all')]
end

MidBi = find(sumzb == max(sumzb));
Bi = StitchedvolumeBi(:,:,MidBi);

BAbz = sum(StitchedvolumeBi(:,:,MidBi),3);

varname=['BAbz',num2str(i)];
assignin('base',varname,BAbz);
%%%%%%%%%%
disp('Generating original master image')
figure, imshow(Ai), colorbar,caxis([0 5000]),title(['Original Master Image',num2str(i)])
saveas(figure(2*i-1),fullfile(DESTINATION_FOLDER,['1_Original_Master_Image',num2str(i),'.png']))
%%%%%%%%%%
disp('Generating original master image')
figure, imshow(Bi), colorbar,caxis([0 25000]),title(['Original Nuclei Image',num2str(i)])
saveas(figure(2*i),fullfile(DESTINATION_FOLDER,['2_Nuclei_Master_Image',num2str(i),'.png']))
end

%% 
repcomb = cat(3,A1,A2,A3);
repcombAbz = cat(3,Abz1,Abz2,Abz3);
repcombBAbz = cat(3,BAbz1,BAbz2,BAbz3);
%Option for filters
for i=1:length(filescell)
prompt = 'Do you want to apply filters? Y/N: ';
x = input(prompt,'s');
%%%%
fla = 'N';
while fla == 'N'

%Condition for filters or not
if x == 'Y'
    prompt = {'Threshold Adjusting Factor:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'18000'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    
  
    Adjustfactori = str2num(answer{1});

    MedAi = medfilt2(repcombBAbz(:,:,i),[4 4]);
    
    MaxAi = nlfilter(MedAi,[4 4],@(x)max(x(:)));
    
    Threshi = MaxAi;
    CustomMedi = median(MaxAi,'all')+Adjustfactori;
    Threshi(Threshi<CustomMedi) = 0;
    Threshi(Threshi>=CustomMedi) = 1;

    flamsi = imfill(Threshi,'holes');
    seDAi = strel('diamond',1);
    BWfinalBi = imerode(flamsi, seDAi);
    BWfinalBi = imerode(BWfinalBi, seDAi);
    
varname=['BWfinalB',num2str(i)];
assignin('base',varname,BWfinalBi);
    
disp('Generating segmented master image')
figure, imshow(BWfinalBi), colorbar,caxis([0 1]),title(['Segmented Master Image',num2str(i)])
saveas(figure(i+6),fullfile(DESTINATION_FOLDER,['3_Segmented_Master_Image',num2str(i),'.png']))

elseif x == 'N'
    prompt = {'Area Threshold:','Edge Sensitivity:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'4000','0.15'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    Areathresh = str2num(answer{1});
    [~,thresholdBi] = edge(repcombBAbz(:,:,i), 'sobel');
    
    fudgeFactorBi = str2num(answer{2});
    
    Threshi = edge(repcombBAbz(:,:,i),'sobel',thresholdBi*fudgeFactorBi);
    
    se90Ai = strel('line',3,90);
    se0Ai = strel('line',3,0);
    BWsdilBi = imdilate(Threshi,[se90Ai,se0Ai]);
    BWdfillBi = imfill(BWsdilBi, 'holes');
    
    seDAi = strel('diamond',1);
    BWfinalBi = imerode(BWdfillBi, seDAi);
    BWfinalBi = imerode(BWfinalBi, seDAi);
else
    disp('Invalid input: Please run the code again')
    return
end

prompt = 'Are you happy? Y/N: ';
fla = input(prompt,'s');
if fla == 'N'
   close
end
end
end
%% 
segcomb = cat(3,BWfinalB1,BWfinalB2,BWfinalB3);
%Option for connected sample
for i = 1:length(filescell)

segfla = 'N';
while segfla == 'N'

    prompt = 'Are your cells connected after segmentation? Y/N: ';
    y = input(prompt,'s');

if y == 'Y'
    prompt = {'Segmentation Factor:','Area Threshold:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'5','4000'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    
    Segmentationfactor = str2num(answer{1});
    Areathresh = str2num(answer{2});

   
    Sobelsegmentation = segcomb(:,:,i);
    segmentopen = ~bwareaopen(~Sobelsegmentation,10);
    %%%
    D = -bwdist(~Sobelsegmentation);
    %%%
    disp('First Watershed')
    FirstWS = watershed(D);
    %%%
    segmentopen = Sobelsegmentation;
    segmentopen(FirstWS == 0) = 0;
    %%%
    mask = imextendedmin(D,Segmentationfactor);
    %%%
    disp('Second Watershed')
    D2 = imimposemin(D,mask);
    Ld2 = watershed(D2);
    bw3 = Sobelsegmentation;
    bw3(Ld2 == 0) = 0;
    BWfinalBi = bw3;
elseif y == 'N'
  BWfinalBi = segcomb(:,:,i);
  prompt = {'Area Threshold:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'4000'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
   
    Areathresh = str2num(answer{1});
end
%%%%%%%% Finding objects in the segment and filtering using area

CC = bwconncomp(BWfinalBi);
L = labelmatrix(CC);
stats = regionprops('table',L,'Area');

varname=['L',num2str(i)];
assignin('base',varname,L);


numstats = table2array(stats);
label = find(numstats(:)>Areathresh);
f = transpose(1:length(label));
AMI = numstats(label,:);
combine = [f,AMI];

filareas = ismember(L,label);

varname=['label',num2str(i)];
assignin('base',varname,label);

varname=['areas',num2str(i)];
assignin('base',varname,filareas);

disp('Generating Filtered Segments')
figure,imshow(filareas),colorbar,caxis([0 1]),title(['Filtered Segments',num2str(i)])
saveas(figure(i+9),fullfile(DESTINATION_FOLDER,['4_Filtered_Segments',num2str(i),'.png']))

prompt = 'Are you happy? Y/N: ';
segfla = input(prompt,'s');
if segfla == 'N'
   close
end
end
end

%% 
%%%%% Puncta Detection
Lcomb = cat(3,L1,L2,L3);
filareascomb = cat(3,areas1,areas2,areas3);
close all
for i = 1:length(filescell)
    
hist = histogram(repcomb(:,:,i),256);
vals = hist.Values;
close
vals(vals==0) = [];
vals = transpose(vals);
vals = log10(vals);
lenv = length(vals);
lindex = transpose([1:lenv]);
dvals = [lindex,vals];


x = lindex;
y = vals;


p = polyfit(x,y,1);



y1 = polyval(p,x);
figure
plot(x,y,'o')
hold on
plot(x,y1)

flango = p(1).*x(:,1)+p(2);
difz = flango-y;
normie = max(difz);
binpoint = find(difz==normie);


Punk = 'N';

while Punk == 'N'
koolaid = repcomb(:,:,i);
PunctaThresh = binpoint*255.996;
prompt = {'PunctaThresh Adjustment:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'0.1'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    
    PTAdjust = str2num(answer{1});
%%%%%%%%%%%%%%%%%%%%%
koolaid(koolaid<PunctaThresh+PunctaThresh*(PTAdjust)) = 0;
close

figure,imshow(koolaid),colorbar,caxis([0 1])
saveas(figure(1),fullfile(DESTINATION_FOLDER,['5_Detected_Puncta',num2str(i),'.png']))


prompt = 'Are you happy? Y/N: ';
Punk = input(prompt,'s');

if Punk == 'N' 
    close
end

end
%%%%%%%%%%%%Generate Punctate Mask 

Dink = 'N';

while Dink == 'N'
prompt = {'Punctate Area Threshold:','Punctate Maximum Size'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'10','40'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    
    Areathresh2 = str2num(answer{1});
    Areathresh3 = str2num(answer{2});

Pseglog = logical(koolaid);
Pseg = bwconncomp(Pseglog);
LabelPseg = labelmatrix(Pseg);

Punctastats = regionprops('table',LabelPseg,repcomb(:,:,i),'Area','MeanIntensity');

Pnumstats = table2array(Punctastats);
plabel = find(Pnumstats(:,1)>Areathresh2 & Pnumstats(:,1)<Areathresh3);
o = transpose(1:length(plabel));
pAMI = Pnumstats(plabel,:);
pcombine = [o,pAMI];

pareas1 = ismember(LabelPseg,plabel);
varname=['pareas1',num2str(i)];
assignin('base',varname,pareas1);

figure,imshow(pareas1),colorbar,caxis([0 1])
saveas(figure(2),fullfile(DESTINATION_FOLDER,['6_Filtered_Detected_Puncta',num2str(i),'.png']))

prompt = 'Are you happy? Y/N: ';
Dink = input(prompt,'s');

if Dink == 'N' 
    close
end

end
close all
end

%% 
pareacomb = cat(3,pareas11,pareas12,pareas13);

for i=1:length(filescell)
pfinalstats = regionprops('table',pareacomb(:,:,i),repcomb(:,:,i),'Centroid','Area','MeanIntensity');
pfinalstatsarr = table2array(pfinalstats);

varname=['pfinalstatsarr',num2str(i)];
assignin('base',varname,pfinalstatsarr);



FINALTABLE = table(transpose(1:length(pfinalstatsarr)),pfinalstatsarr(:,1),pfinalstatsarr(:,4),'VariableNames',{'Label','Area(Px)','MeanIntensity'});
writetable(FINALTABLE,fullfile(DESTINATION_FOLDER,['Puncta_Statistics',num2str(i),'.csv']));

varname=['FinalTable',num2str(i)];
assignin('base',varname,FINALTABLE);
end

%% 
areascomb = cat(3,areas1,areas2,areas3);

for i = 1:length(filescell)
structarea = bwconncomp(areascomb(:,:,i));
labmat = labelmatrix(structarea);
    

structarea = bwconncomp(areascomb(:,:,i));
Cellz = [1:structarea.NumObjects];

varname=['DataPuncta',num2str(i)];
assignin('base',varname,Cellz);

DataPunctaTab = table(Cellz(:),'VariableNames',{'Cell Labels'});
writetable(DataPunctaTab,fullfile(DESTINATION_FOLDER,['Cell Labels',num2str(i),'.csv']));

end
%% 
