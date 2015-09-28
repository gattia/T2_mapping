 
clc 
clear

cd('/Volumes/Anthony.Gatti_MacintoshHD/Users/Gatti/Desktop/T2_map_testing')

file=[dir('exam*')];

for e=1:length(file);
    files=file(char(e));
    exam(e,1)=str2num(files.name(6:9));
end
echos=8

for t=1:length(exam);
    
    f=exam(t)
    eval(['T2orig_' num2str(f) '=nrrdread(''' 'T2MAP_' num2str(f) '.nrrd' '''); ']);
    eval(['cd exam_' num2str(f) ])
    Ser=dir('Ser*00');

    if Ser.name =='Ser300'
        series=3
    elseif Ser.name =='Ser400'
        series=4
    elseif Ser.name =='Ser500'
        series=5
    elseif Ser.name =='Ser600'
        series=6
    elseif Ser.name =='Ser700'
        series=7        
    elseif Ser.name =='Ser800'
        series=8
    elseif Ser.name =='Ser900'
        series=9
    end
    
    
    eval(['cd Ser' num2str(series)]);
    
    images=dir('E*');

    numsl=length(images)/echos;




imref=[1:numsl];
index=0; 
    for h=1:numsl,
        for te=1:echos,
            index=echos*(h-1)+imref(te)

              if index <= 9,
                 eval(['II_' num2str(f) '(:,:,h,te)=dicomread(''' 'E' num2str(f) 'S' num2str(series) 'I00' num2str(index) '.MR.dcm' '''); ']);
              elseif index <= 99,
                 eval(['II_' num2str(f) '(:,:,h,te)=dicomread(''' 'E' num2str(f) 'S' num2str(series) 'I0' num2str(index) '.MR.dcm' '''); ']);
              else,
                 eval(['II_' num2str(f) '(:,:,h,te)=dicomread(''' 'E' num2str(f) 'S' num2str(series) 'I' num2str(index) '.MR.dcm' '''); ']);
              end;

            if h==1,
                % This is just to read the headers for the correct TE times!
                % I did this because I'm not sure if some oddity could change it
                eval(['info=dicominfo(''' 'E' num2str(f) 'S' num2str(series) 'I00' num2str(index) '.MR.dcm' '''); ']);
                eval(['TE_' num2str(f) '(te)=info.EchoTime;']);
            end;
        end;
    end;
    


    cd ..
    cd ..



%% 

% SO at this point you have 2 matrices- one is vector of TEs the other is a
%     monster huge 4D matrix.
% Now need to do regression analysis one at a time (yup slow!)
% So need a triple loop of X, Y and slice and regress S vs TE
% Where, S is the vector defined by the 4th variable in II (x,y,slnum,:)
% The result will be 2 3D matrices called T2, PD and Rsq with (x,y,slnum)
% T2 is the slope, PD (proton density) is intercept and Rsq error

% Also- 2 ways to speed things up:
% 1.  Don't fit the noise
% 2.  Is there a point to fitting T2 over the whole image?  I don't think
% so.  So, define a region to fit the T2s over.

% Define the region you want to fit, bsaed on image to choose for zooming

figure(1);
for qq=1:numsl,
    subplot(5,6,qq); 
    eval (['imagesc(II_' num2str(f) '(:,:,qq,1));' ]); 
    axis('image'); 
    axis off; 
    colormap('gray');
end
set(gcf,'units','normalized','outerposition',[0 0 1 1])

figch=input('Which image would you like to zoom on (1 is top left, the last one is bottom right)? ');
close(1);

figure(1); 
eval (['imagesc(II_' num2str(f) '(:,:,figch,1));' ]); 
axis('image'); 
colormap('gray');

disp('click mouse LEFT on top left then on bottom right for ROI');
disp(' ');
  disp('CHOOSE A HUGE ROI TO ENVELOPE ALL AREA OF INTEREST !!');
  disp(' ');
  disp('T2s will be calculated over area of interest, on every slice');

[xx1,yy1,button]=ginput(1);
if (button == 1),
  [xx2,yy2,button]=ginput(1);
end;

hold on;
set(figure(1),'DefaultLineColor','red');
vecx=[xx1,xx1,xx2,xx2];
vecy=[yy1,yy2,yy2,yy1];
% fill(vecx,vecy,'r');
line( [vecx,vecx(1)],[vecy,vecy(1)] );
set(line( [vecx,vecx(1)],[vecy,vecy(1)] ), 'Color',[1 0 0], 'lineWidth', [1.0] ); 

hold off;



% the next 4 lines round the coordinates off to a whole number becuase x1
% x2 y1 y2 are used as indices pointing to specific locations in the II
% matrix.  Matrices have discrete dimensions (i.e. no fractions)
eval (['x1_' num2str(f) '=floor(xx1);' ]);
eval (['x2_' num2str(f) '=floor(xx2);' ]);
eval (['y1_' num2str(f) '=floor(yy1);' ]);
eval (['y2_' num2str(f) '=floor(yy2);' ]);


% extract the zoomed in area out of the large 4D matrix called II

    eval (['IIzoom_' num2str(f) '=II_' num2str(f) '(y1_' num2str(f) ':y2_' num2str(f) ',x1_' num2str(f) ':x2_' num2str(f) ',:,:);' ]);

end


%% CREATE MAPS
cd('/Volumes/Anthony.Gatti_MacintoshHD/Users/Gatti/Desktop/T2_map_testing/Analyzed')
for u=1:length(exam); 
    
    g=exam(u)
% allocate the memory for each fitted parameter to help speed things up
    eval(['numsl=length(IIzoom_' num2str(g) '(1,1,:,1));' ]);    
    
    eval(['TE=TE_' num2str(g) ';' ]);
    
    TE=TE';
    
    eval (['T2=zeros(size(IIzoom_' num2str(g) ',1),size(IIzoom_' num2str(g) ',2),numsl);' ]);
    eval (['PD=zeros(size(IIzoom_' num2str(g) ',1),size(IIzoom_' num2str(g) ',2),numsl);' ]);
    eval (['Rsq=zeros(size(IIzoom_' num2str(g) ',1),size(IIzoom_' num2str(g) ',2),numsl);' ]);
    
    disp('T2 pre fitting started') ;
    disp(g);  
    
    eval (['map=IIzoom_' num2str(g) ';' ]);
    
    T2 = zeros(size(map(:,:,:))) ;
    PD = zeros(size(map(:,:,:))) ;
    Rsq = zeros(size(map(:,:,:))) ;

    
    [c,r,s,e] = ind2sub(size(map),find(map(:,:,5,1)>500 & map(:,:,5,1)<2100)); % c= column (y), r= row (x), s=slice, e=echo
    
    tic
    
    for l = 1:length(c); 
         
        S=squeeze(map(c(l),r(l),s(l),:));
        [T2fit,G]=fit(TE,S,'exp1'); warning('OFF');
        T2(c(l),r(l),s(l))=-1/T2fit.b;  % T2 in milliseconds
        PD(c(l),r(l),s(l))=T2fit.a;  % y-intercept is proton density, PD
        Rsq(c(l),r(l),s(l))=G.rsquare;  % R^2 fit quality  
       
    end
    
    toc
    display ('Slice 5 completed')
    
%     for z=1:numsl,
%         tic
%         for y=1:size(map,1),
%            parfor x=1:size(map,2),
%                 S=squeeze(map(y,x,z,:));
%                 if S(1)<350, %
%                     % So if the signal is likely noise don't fit it and set all
%                     % the parameters to zero.  The signal threshold of 99 is just a guess and can be changed. 
%                     T2(y,x,z)=0; PD(y,x,z)=0; Rsq(y,x,z)=0;
%                 elseif S(1)>3750, %get rid of some more noise. 
%                     T2(y,x,z)=0; PD(y,x,z)=0; Rsq(y,x,z)=0;
%                 else
%                     [T2fit,G]=fit(TE,S,'exp1'); warning('OFF');
%                     T2(y,x,z)=-1/T2fit.b;  % T2 in milliseconds
%                     PD(y,x,z)=T2fit.a;  % y-intercept is proton density, PD
%                     Rsq(y,x,z)=G.rsquare;  % R^2 fit quality
%                 end;
%     %             if T2pre(y,x,z)>300 ,% get rid of noise from the finished map. 
%     %                 T2pre(y,x,z)=0;
%     %                 PDpre(y,x,z)=0;
%     %                 Rsqpre(y,x,z)=0;
%     %             end  
%     %    
%     %             if T2pre(y,x,z)<0 ,% get rid of noise from the finished map. 
%     %                 T2pre(y,x,z)=0;
%     %                 PDpre(y,x,z)=0;
%     %                 Rsqpre(y,x,z)=0;
%     %             end  
%             end;
%         end;
%         toc
%         disp('Slice completed = '); z 
%     end;
    
eval (['T2new_' num2str(g) '=zeros(256,256,numsl);' ]);
eval (['PDnew_' num2str(g) '=zeros(256,256,numsl);' ]);
eval (['Rsqnew_' num2str(g) '=zeros(256,256,numsl);' ]);

eval (['T2new_' num2str(g) '(y1_' num2str(g) ':y2_' num2str(g) ',x1_' num2str(g) ':x2_' num2str(g) ',:)=T2;' ]);
eval (['PDnew_' num2str(g) '(y1_' num2str(g) ':y2_' num2str(g) ',x1_' num2str(g) ':x2_' num2str(g) ',:)=PD;' ]);
eval (['Rsqnew_' num2str(g) '(y1_' num2str(g) ':y2_' num2str(g) ',x1_' num2str(g) ':x2_' num2str(g) ',:)=Rsq;' ]);

eval(['save T2maps_' num2str(g) '.mat T2new_' num2str(g) ' PDnew_' num2str(g) ' Rsqnew_' num2str(g) '' ]);

end;

for v=1:length(exam); 
    
    h=exam(v);
    
    disp(h)
    tic
    
    
    eval(['numsl=length(T2new_' num2str(h) '(1,1,:));' ]); 
   
    T2_convert=zeros(256,256,numsl);
    PD_convert=zeros(256,256,numsl);
    Rsq_convert=zeros(256,256,numsl);
    
    eval (['T2new=T2new_' num2str(h) ';' ]); 
    eval (['PDnew=PDnew_' num2str(h) ';' ]); 
    eval (['Rsqnew=Rsqnew_' num2str(h) ';' ]); 

    
    for x=1:256;
        for y=1:256;
            for z=1:numsl;
                T2_convert(x,y,z)=T2new(y,x,z);
            end
        end
    end

    
    for x=1:256;
        for y=1:256;
            for z=1:numsl;
                PD_convert(x,y,z)=PDnew(y,x,z);
            end
        end
    end
    
    
    for x=1:256;
        for y=1:256;
            for z=1:numsl;
                Rsq_convert(x,y,z)=Rsqnew(y,x,z);
            end
        end
    end
    

    eval (['T2nrrd_' num2str(h) '=T2orig_' num2str(h) ';' ]); 
    eval (['T2nrrd_' num2str(h) '.pixelData=T2_convert;' ]);

    eval (['PDnrrd_' num2str(h) '=T2orig_' num2str(h) ';' ]); 
    eval (['PDnrrd_' num2str(h) '.pixelData=PD_convert;' ]);

    eval (['Rsqnrrd_' num2str(h) '=T2orig_' num2str(h) ';' ]); 
    eval (['Rsqnrrd_' num2str(h) '.pixelData=Rsq_convert;' ]);

    toc
    
end
    

for w=1:length(exam);
    k=exam(w);
    eval(['nrrdwrite(' '''' 'T2map_' num2str(k) '.nrrd''' ', T2nrrd_' num2str(k) '); ']);
    eval(['nrrdwrite(' '''' 'PDmap_' num2str(k) '.nrrd''' ', PDnrrd_' num2str(k) '); ']);
    eval(['nrrdwrite(' '''' 'Rsqmap_' num2str(k) '.nrrd''' ', Rsqnrrd_' num2str(k) '); ']);
end

% cd('/Volumes/Anthony.Gatti Macintosh HD/Users/Gatti/Desktop/MRI/MRI/ToAnalyze')
% 
% for x=1:length(exam);
%     l=exam(x);
%     eval(['movefile(''' 'T2MAP_' num2str(l) '.nrrd' ''',''' '/Volumes/Anthony.Gatti Macintosh HD/Users/Gatti/Desktop/MRI/MRI/Analyzed/exam_folder' '''); ']);   
%     eval(['movefile(''' 'exam_' num2str(l) ''',''' '/Volumes/Anthony.Gatti Macintosh HD/Users/Gatti/Desktop/MRI/MRI/Analyzed/exam_folder' '''); ']);
% end
% 
% cd('/Volumes/Anthony.Gatti Macintosh HD/Users/Gatti/Desktop/MRI/MRI/Analyzed')
% 
% for x=1:length(exam);
%     l=exam(x);
%     eval(['copyfile(''' 'T2map_' num2str(l) '.nrrd' ''',''' '/Volumes/Anthony.Gatti Macintosh HD/Users/Gatti/Desktop/MRI/MRI/Analyzed/exam_folder/exam_' num2str(l) '''); ']);   
%     eval(['copyfile(''' 'Rsqmap_' num2str(l) '.nrrd' ''',''' '/Volumes/Anthony.Gatti Macintosh HD/Users/Gatti/Desktop/MRI/MRI/Analyzed/exam_folder/exam_' num2str(l) '''); ']);   
%     eval(['copyfile(''' 'PDmap_' num2str(l) '.nrrd' ''',''' '/Volumes/Anthony.Gatti Macintosh HD/Users/Gatti/Desktop/MRI/MRI/Analyzed/exam_folder/exam_' num2str(l) '''); ']);   
% end